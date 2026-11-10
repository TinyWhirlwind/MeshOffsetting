#include "SSVH.h"
#include <set>
#include <vector>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>
namespace meshoffsetting
{
	SSVH::SSVH(const CMeshO& mesh)
		:_mesh(mesh)
	{
		_nodes.clear();
		_root = -1;
		Build(mesh);
	}

	void SSVH::Build(const CMeshO& mesh)
	{
		_mesh = mesh;
		_faceIndices.clear();
		_faceIndices.reserve(mesh.face.size());
		for (int i = 0; i < mesh.fn; ++i)
		{
			if (mesh.face[i].IsD())
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
		//暂时未实现
		return ClosestPointResult();
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
		auto middle = std::partition(
			_faceIndices.begin(),
			_faceIndices.end(),
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

	RSS SSVH::FitRSS(std::size_t begin, std::size_t end, Point3m mean_p)
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

				Eigen::Vector3<Scalarm> delta = pi - mean;
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
		Scalarm height = std::numeric_limits<Scalarm>::lowest();
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
		bound.sweepRadius = (max_coord[2] - min_coord[2]) / 2;
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