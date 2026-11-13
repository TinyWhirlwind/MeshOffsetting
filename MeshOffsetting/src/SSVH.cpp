#include "SSVH.h"
#include "Tolerance.h"
#include <algorithm>
#include <Eigen/Dense>
namespace meshoffsetting
{
	SSVH::SSVH(const CMeshO& mesh)
		:_mesh(mesh)
	{
		Build();
	}

	void SSVH::Build()
	{
		_nodes.clear();
		_faceIndices.clear();
		_faceIndices.reserve(_mesh.face.size());
		for (int i = 0; i < _mesh.fn; ++i)
		{
			if (_mesh.face[i].IsD())
				continue;
			_faceIndices.push_back(i);

		}
		if (_faceIndices.empty())
		{
			_root = -1;
			return;
		}
		int fn = _faceIndices.size();
		_nodes.reserve(fn * 2 - 1);
		_nodes.emplace_back();
		_root = 0;
		BuildNode(_root, 0, fn);
	}

	ClosestPointResult SSVH::QueryClosestPoint(const Point3m& q)
	{
		ClosestPointResult result;
		if (!Empty())
		{
			QueryNode(_root, q, result);
		}
		return result;
	}

	void SSVH::QueryNode(int node_index, const Point3m& query_p, ClosestPointResult& result)
	{
		const SSVHNode& node = _nodes[node_index];
		if (DistanceToRSS(query_p, node.bound) >= result.unsigned_distance)
		{
			return;
		}
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
		Scalarm left_dis = DistanceToRSS(query_p, _nodes[left_child].bound);
		Scalarm right_dis = DistanceToRSS(query_p, _nodes[right_child].bound);
		if (left_dis > right_dis)
		{
			std::swap(left_child, right_child);
			std::swap(left_dis, right_dis);
		}

		if (left_dis < result.unsigned_distance)
			QueryNode(left_child, query_p, result);
		if (right_dis < result.unsigned_distance)
			QueryNode(right_child, query_p, result);
	}

	Scalarm SSVH::DistanceToRSS(const Point3m& query_p, const RSS& bound) const
	{
		Point3m delta = query_p - bound.origin;
		Scalarm u = delta * bound.axisU;
		Scalarm v = delta * bound.axisV;
		Scalarm n = delta * bound.axisN;

		Scalarm du = 0.0;
		Scalarm dv = 0.0;
		if (u < 0)
			du = -u;
		else if (u > bound.lengthU)
			du = u - bound.lengthU;

		if (v < 0)
			dv = -v;
		else if (v > bound.lengthV)
			dv = v - bound.lengthV;
		Scalarm rectangleDistance = std::sqrt(du * du + dv * dv + n * n);
		return std::max(Scalarm(0), rectangleDistance - bound.sweep_radius);
	}


	void SSVH::UpdateClosestPointWithTriangle(int face_index, const Point3m& query_p, ClosestPointResult& result)
	{
		const auto& face = _mesh.face[face_index];
		const Point3m a = face.P(0);
		const Point3m b = face.P(1);
		const Point3m c = face.P(2);

		auto setResult = [&](Scalarm s, Scalarm t)
			{
				Point3m closest_point =
					a + (b - a) * s + (c - a) * t;

				Point3m diff = query_p - closest_point;
				Scalarm distance = std::sqrt(diff * diff);

				if (distance < result.unsigned_distance)
				{
					result.found = true;
					result.u = s;
					result.v = t;
					result.faceId = face_index;
					result.closest_point = closest_point;
					result.unsigned_distance = distance;
				}
			};

		const Point3m ab = b - a;
		const Point3m ac = c - a;
		const Point3m ap = query_p - a;

		const Scalarm d1 = ab * ap;
		const Scalarm d2 = ac * ap;
		if (d1 <= 0 && d2 <= 0)
		{
			setResult(0, 0);
			//result.closestType = ClosestType::Vertex;
			//result.closestFeature = 0; // A
			return;
		}

		const Point3m bp = query_p - b;
		const Scalarm d3 = ab * bp;
		const Scalarm d4 = ac * bp;
		if (d3 >= 0 && d4 <= d3)
		{
			setResult(1, 0);
			//result.closestType = ClosestType::Vertex;
			//result.closestFeature = 1; // B
			return;
		}

		const Scalarm vc = d1 * d4 - d3 * d2;
		if (vc <= 0 && d1 >= 0 && d3 <= 0)
		{
			const Scalarm denom = d1 - d3;
			const Scalarm v = denom > kScalarTolerance ? d1 / denom : 0;
			setResult(v, 0);
			//result.closestType = ClosestType::Edge;
			//result.closestFeature = 0; // AB
			return ;
		}

		const Point3m cp = query_p - c;
		const Scalarm d5 = ab * cp;
		const Scalarm d6 = ac * cp;
		if (d6 >= 0 && d5 <= d6)
		{
			setResult(0, 1);
			//result.closestType = ClosestType::Vertex;
			//result.closestFeature = 2; // C
			return;
		}

		const Scalarm vb = d5 * d2 - d1 * d6;
		if (vb <= 0 && d2 >= 0 && d6 <= 0)
		{
			const Scalarm denom = d2 - d6;
			const Scalarm w = denom > kScalarTolerance ? d2 / denom : 0;
			setResult(0, w);
			//result.closestType = ClosestType::Edge;
			//result.closestFeature = 2; // AC
			return;
		}

		const Scalarm va = d3 * d6 - d5 * d4;
		if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
		{
			const Scalarm denom = (d4 - d3) + (d5 - d6);
			const Scalarm w = denom > kScalarTolerance ? (d4 - d3) / denom : 0;
			setResult(1 - w, w);
			//result.closestType = ClosestType::Edge;
			//result.closestFeature = 1; // BC
			return;
		}

		const Scalarm denom = va + vb + vc;
		if (std::abs(denom) > kScalarTolerance)
		{
			const Scalarm v = vb / denom;
			const Scalarm w = vc / denom;
			setResult(v, w);
			//result.closestType = ClosestType::Face;
			//result.closestFeature = 0;
		}
		else
		{
			auto setSegmentResult = [&](const Point3m& p0, const Point3m& p1, Scalarm s0, Scalarm t0, Scalarm s1, Scalarm t1)
				{
					const Point3m edge = p1 - p0;
					const Scalarm len2 = edge * edge;
					Scalarm u = len2 > kScalarTolerance ? ((query_p - p0) * edge) / len2 : 0;
					u = std::clamp(u, Scalarm(0), Scalarm(1));

					ClosestPointResult candidate{};
					candidate.faceId = face_index;
					//candidate._sign = 0;
					//candidate.intersected = false;
					candidate.u = s0 + (s1 - s0) * u;
					candidate.v = t0 + (t1 - t0) * u;
					candidate.closest_point = p0 + edge * u;
					const Point3m diff = query_p - candidate.closest_point;
					candidate.unsigned_distance = std::sqrt(diff * diff);
					//candidate.intersected = candidate.unsigned_distance <= kScalarTolerance;
					if (candidate.unsigned_distance < result.unsigned_distance)
					{
						result = candidate;
					}
				};

			setSegmentResult(a, b, 0, 0, 1, 0);
			setSegmentResult(a, c, 0, 0, 0, 1);
			setSegmentResult(b, c, 1, 0, 0, 1);
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
				Point3m center = FaceCentroid(id);
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
			const CFaceO& fi = _mesh.face[fid];
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
			const CFaceO& fi = _mesh.face[fid];
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

	Point3m SSVH::FaceCentroid(int fid)
	{
		auto& f = _mesh.face[fid];
		Point3m center{ 0,0,0 };
		for (int i = 0; i < f.VN(); ++i)
		{
			center += f.P(i);
		}
		center /= (Scalarm)(f.VN());
		return center;
	}
}