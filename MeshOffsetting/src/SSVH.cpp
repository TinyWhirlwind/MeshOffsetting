#include "SSVH.h"
#include "Tolerance.h"
#include <algorithm>
#include <Eigen/Dense>
namespace meshoffsetting
{
	SSVH::SSVH(const CMeshO& mesh)
	{
		Initialize(mesh);
	}

	void SSVH::Initialize(const CMeshO& mesh)
	{
		_mesh = &mesh;
		Build();
	}

	void SSVH::Build()
	{
		_nodes.clear();
		_faceIndices.clear();
		_triangles.clear();
		_root = -1;

		if (_mesh == nullptr)
		{
			return;
		}
		_faceIndices.reserve(_mesh->face.size());
		_triangles.resize(_mesh->face.size());

		for (std::size_t face_id = 0; face_id < _mesh->face.size(); ++face_id)
		{
			const CFaceO& face = _mesh->face[face_id];
			if (face.IsD())
				continue;
			_faceIndices.push_back(static_cast<int>(face_id));

			CachedTriangle& triangle = _triangles[face_id];
			triangle.a = face.P(0);
			triangle.ab = face.P(1) - triangle.a;
			triangle.ac = face.P(2) - triangle.a;
		}

		if (_faceIndices.empty())
		{
			return;
		}
		std::size_t faceCount = _faceIndices.size();
		_nodes.reserve(faceCount * 2 - 1);
		_nodes.emplace_back();
		_root = 0;
		BuildNode(_root, 0, faceCount);
	}

	ClosestPointResult SSVH::QueryClosestPoint(const Point3m& q) const
	{
		ClosestPointResult result;
		if (!Empty())
		{
			QueryNode(_root, q, result);
		}
		return result;
	}

	void SSVH::QueryNode(int node_index, const Point3m& query_p, ClosestPointResult& result) const
	{
		const SSVHNode& node = _nodes[node_index];
		if (node.isLeaf())
		{
			for (int i = 0; i < node.faceCount; ++i)
			{
				int id = node.firstChild + i;
				int face_id = _faceIndices[id];

				UpdateClosestPointWithTriangle(face_id, query_p, result);
			}
			return;
		}

		int left_child = node.firstChild;
		int right_child = node.firstChild + 1;
		Scalarm left_core_squared = DistanceToRSSCoreSquared(query_p, _nodes[left_child].bound);
		Scalarm right_core_squared = DistanceToRSSCoreSquared(query_p, _nodes[right_child].bound);
		if (left_core_squared > right_core_squared)
		{
			std::swap(left_child, right_child);
			std::swap(left_core_squared, right_core_squared);
		}

		if (result.faceId < 0)
		{
			QueryNode(left_child, query_p, result);
		}
		else
		{
			Scalarm limit = result.unsigned_distance + _nodes[left_child].bound.sweep_radius;
			if (left_core_squared <= limit * limit)
			{
				QueryNode(left_child, query_p, result);
			}
		}

		// 左子树查询后，result 可能已经更新，所以这里必须重新计算右子树的 limit。
		if (result.faceId < 0)
		{
			QueryNode(right_child, query_p, result);
		}
		else
		{
			Scalarm limit =result.unsigned_distance + _nodes[right_child].bound.sweep_radius;
			if (right_core_squared <= limit * limit)
			{
				QueryNode(right_child, query_p, result);
			}
		}
	}

	Scalarm SSVH::DistanceToRSSCoreSquared(const Point3m& query_p,const RSS& bound) const
	{
		Point3m delta = query_p - bound.origin;

		Scalarm u = delta * bound.axisU;
		Scalarm v = delta * bound.axisV;
		Scalarm n = delta * bound.axisN;

		Scalarm du = 0;
		Scalarm dv = 0;

		if (u < 0)
			du = -u;
		else if (u > bound.lengthU)
			du = u - bound.lengthU;

		if (v < 0)
			dv = -v;
		else if (v > bound.lengthV)
			dv = v - bound.lengthV;

		return du * du + dv * dv + n * n;
	}


	void SSVH::UpdateClosestPointWithTriangle(int face_index, const Point3m& query_p, ClosestPointResult& result) const
	{
		const CachedTriangle& triangle = _triangles[face_index];
		const Point3m& a = triangle.a;
		const Point3m& ab = triangle.ab;
		const Point3m& ac = triangle.ac;
		const Point3m b = a + ab;
		const Point3m c = a + ac;

		auto setResult = [&](Scalarm s, Scalarm t, ClosestFeature feature)
			{
				Point3m closest_point = a + ab * s + ac * t;
				Point3m diff = query_p - closest_point;
				Scalarm distance_squared = diff * diff;

				if (distance_squared < result.unsigned_distance_squared)
				{
					result.feature = feature;
					result.faceId = face_index;
					result.closest_point = closest_point;
					result.unsigned_distance_squared = distance_squared;
					result.unsigned_distance = std::sqrt(distance_squared);
				}
			};

		const Point3m ap = query_p - a;

		const Scalarm d1 = ab * ap;
		const Scalarm d2 = ac * ap;
		if (d1 <= 0 && d2 <= 0)
		{
			setResult(0, 0, ClosestFeature::VertexA);
			return;
		}

		const Point3m bp = query_p - b;
		const Scalarm d3 = ab * bp;
		const Scalarm d4 = ac * bp;
		if (d3 >= 0 && d4 <= d3)
		{
			setResult(1, 0, ClosestFeature::VertexB);
			return;
		}

		const Scalarm vc = d1 * d4 - d3 * d2;
		if (vc <= 0 && d1 >= 0 && d3 <= 0)
		{
			const Scalarm denom = d1 - d3;
			Scalarm v = denom > kScalarTolerance ? d1 / denom : 0;
			v = std::clamp(v, Scalarm(0), Scalarm(1));
			ClosestFeature feature = ClosestFeature::EdgeAB;
			if (v == Scalarm(0))
				feature = ClosestFeature::VertexA;
			else if (v == Scalarm(1))
				feature = ClosestFeature::VertexB;
			setResult(v, 0, feature);
			return ;
		}

		const Point3m cp = query_p - c;
		const Scalarm d5 = ab * cp;
		const Scalarm d6 = ac * cp;
		if (d6 >= 0 && d5 <= d6)
		{
			setResult(0, 1, ClosestFeature::VertexC);
			return;
		}

		const Scalarm vb = d5 * d2 - d1 * d6;
		if (vb <= 0 && d2 >= 0 && d6 <= 0)
		{
			const Scalarm denom = d2 - d6;
			Scalarm w = denom > kScalarTolerance ? d2 / denom : 0;
			w = std::clamp(w, Scalarm(0), Scalarm(1));
			ClosestFeature feature = ClosestFeature::EdgeCA;
			if (w == Scalarm(0))
				feature = ClosestFeature::VertexA;
			else if (w == Scalarm(1))
				feature = ClosestFeature::VertexC;
			setResult(0, w, feature);
			return;
		}

		const Scalarm va = d3 * d6 - d5 * d4;
		if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
		{
			const Scalarm denom = (d4 - d3) + (d5 - d6);
			Scalarm w = denom > kScalarTolerance ? (d4 - d3) / denom : 0;
			w = std::clamp(w, Scalarm(0), Scalarm(1));
			ClosestFeature feature = ClosestFeature::EdgeBC;
			if (w == Scalarm(0))
				feature = ClosestFeature::VertexB;
			else if (w == Scalarm(1))
				feature = ClosestFeature::VertexC;
			setResult(1 - w, w, feature);
			return;
		}

		const Scalarm denom = va + vb + vc;
		if (std::abs(denom) > kScalarTolerance)
		{
			const Scalarm v = vb / denom;
			const Scalarm w = vc / denom;
			setResult(v, w, ClosestFeature::Face);
		}
		else
		{
			auto setSegmentResult = [&](const Point3m& p0,
				const Point3m& p1,
				ClosestFeature start_feature,
				ClosestFeature end_feature,
				ClosestFeature edge_feature)
				{
					const Point3m edge = p1 - p0;
					const Scalarm len2 = edge * edge;
					Scalarm u = len2 > kScalarTolerance ? ((query_p - p0) * edge) / len2 : 0;
					u = std::clamp(u, Scalarm(0), Scalarm(1));

					ClosestPointResult candidate{};
					candidate.faceId = face_index;
					if (u == Scalarm(0))
					{
						candidate.feature = start_feature;
					}
					else if (u == Scalarm(1))
					{
						candidate.feature = end_feature;
					}
					else
					{
						candidate.feature = edge_feature;
					}
					candidate.closest_point = p0 + edge * u;
					const Point3m diff = query_p - candidate.closest_point;
					const Scalarm distance_squared = diff * diff;
					if (distance_squared < result.unsigned_distance_squared)
					{
						candidate.unsigned_distance_squared = distance_squared;
						candidate.unsigned_distance = std::sqrt(distance_squared);
						result = candidate;
					}
				};

			setSegmentResult(a, b,
				ClosestFeature::VertexA,
				ClosestFeature::VertexB,
				ClosestFeature::EdgeAB);
			setSegmentResult(a, c,
				ClosestFeature::VertexA,
				ClosestFeature::VertexC,
				ClosestFeature::EdgeCA);
			setSegmentResult(b, c,
				ClosestFeature::VertexB,
				ClosestFeature::VertexC,
				ClosestFeature::EdgeBC);
		}
	}

	bool SSVH::Empty() const
	{
		return _root < 0 || _nodes.empty();
	}

	void SSVH::BuildNode(int nodeIndex, std::size_t begin, std::size_t end)
	{
		Point3m mean;
		RSS bound = FitRSS(begin, end, mean);
		_nodes[nodeIndex].bound = bound;
		std::size_t faceCount = end - begin;
		if (faceCount == 1)
		{
			_nodes[nodeIndex].firstChild = static_cast<int>(begin);
			_nodes[nodeIndex].faceCount = 1;
			return;
		}

		Point3m splitAxis = bound.axisU;
		auto first = _faceIndices.begin() + begin;
		auto last = _faceIndices.begin() + end;
		auto middle = std::partition(first, last,
			[&](int id)
			{
				const CachedTriangle& triangle = _triangles[id];
				const Point3m center = triangle.a +
					(triangle.ab + triangle.ac) / Scalarm(3);
				return (center - mean) * splitAxis <= 0;
			});
		std::size_t mid = middle - _faceIndices.begin();

		// 所有三角形都被分到同一侧时，PQP 直接从数量中间强制分开
		if (mid == begin || mid == end)
		{
			mid = begin + faceCount / 2;
		}
		int firstChild = static_cast<int>(_nodes.size());
		_nodes.emplace_back();
		_nodes.emplace_back();
		_nodes[nodeIndex].firstChild = firstChild;
		_nodes[nodeIndex].faceCount = 0;
		BuildNode(firstChild, begin, mid);
		BuildNode(firstChild + 1, mid, end);
	}

	RSS SSVH::FitRSS(std::size_t begin, std::size_t end, Point3m& mean_p)
	{
		Eigen::Vector3<Scalarm> mean = Eigen::Vector3<Scalarm>::Zero();
		Eigen::Matrix3<Scalarm> covariance = Eigen::Matrix3<Scalarm>::Zero();
		int vn = 0;
		for (int i = begin; i < end; ++i)
		{
			int fid = _faceIndices[i];
			const CFaceO& fi = _mesh->face[fid];
			for (int j = 0; j < fi.VN(); ++j)
			{
				const Point3m& p = fi.P(j);
				Eigen::Vector3<Scalarm> pi(p.X(), p.Y(), p.Z());
				mean += pi;
				covariance += (pi * pi.transpose());
				++vn;
			}
		}
		mean /= vn;
		covariance = covariance / vn - mean * mean.transpose();

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3<Scalarm>> solver(covariance);
		Eigen::Vector3<Scalarm> axis[3];
		axis[0] = solver.eigenvectors().col(2).normalized();
		axis[1] = solver.eigenvectors().col(1).normalized();
		axis[2] = axis[0].cross(axis[1]).normalized();

		RSS bound;
		Scalarm height = std::numeric_limits<Scalarm>::max();
		Scalarm lowest = std::numeric_limits<Scalarm>::lowest();
		Scalarm max_coord[3] = { lowest, lowest, lowest };
		Scalarm min_coord[3] = { height, height, height };
		for (int i = begin; i < end; ++i)
		{
			int fid = _faceIndices[i];
			const CFaceO& fi = _mesh->face[fid];
			for (int j = 0; j < fi.VN(); ++j)
			{
				const Point3m& p = fi.P(j);
				Eigen::Vector3<Scalarm> pi(p.X(), p.Y(), p.Z());

				for (int k = 0; k < 3; ++k)
				{
					min_coord[k] = std::min(min_coord[k], (pi - mean).dot(axis[k]));
					max_coord[k] = std::max(max_coord[k], (pi - mean).dot(axis[k]));
				}
			}
		}
		bound.axisU = Point3m{ axis[0].x(),axis[0].y(),axis[0].z() };
		bound.axisV = Point3m{ axis[1].x(),axis[1].y(),axis[1].z() };
		bound.axisN = Point3m{ axis[2].x(),axis[2].y(),axis[2].z() };
		bound.lengthU = max_coord[0] - min_coord[0];
		bound.lengthV = max_coord[1] - min_coord[1];
		bound.sweep_radius = (max_coord[2] - min_coord[2]) / 2;
		Scalarm centerZ = (max_coord[2] + min_coord[2]) / 2;
		mean_p = Point3m(mean.x(), mean.y(), mean.z());
		bound.origin = mean_p
			+ bound.axisU * min_coord[0] 
			+ bound.axisV * min_coord[1]
			+ bound.axisN * centerZ;
		return bound;
	}
}
