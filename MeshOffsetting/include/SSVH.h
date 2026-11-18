/**
 * @file SSVH.h
 * @author zhangpingping
 * @date 2026-11-07
 * @brief
 * 最近三角形查询,包围体：矩形扫掠球
 */
#ifndef MESHOFFSETTING_SSVH_H
#define MESHOFFSETTING_SSVH_H
#include "CMeshO.h"
#include <vector>
#include <limits>
namespace meshoffsetting
{
	enum class ClosestFeature
	{
		Unknown,
		Face,
		EdgeAB,
		EdgeBC,
		EdgeCA,
		VertexA,
		VertexB,
		VertexC
	};

	struct ClosestPointResult
	{
		int faceId = -1;
		Point3m closest_point;
		Scalarm unsigned_distance = std::numeric_limits<Scalarm>::max();
		Scalarm unsigned_distance_squared = std::numeric_limits<Scalarm>::max();
		ClosestFeature feature = ClosestFeature::Unknown;
	};
	
	struct RSS
	{
		Scalarm sweep_radius = 0;
		Point3m origin;

		Point3m axisU;
		Point3m axisV;
		Point3m axisN;

		Scalarm lengthU = 0;
		Scalarm lengthV = 0;
	};

	struct CachedTriangle
	{
		Point3m a;
		Point3m ab;
		Point3m ac;
	};

	struct SSVHNode
	{
		RSS bound;
		int firstChild = -1;
		int faceCount = 0;
		bool isLeaf() const
		{
			return faceCount > 0;
		}
	};
	class SSVH
	{
	public:
		explicit SSVH(const CMeshO& mesh);
		SSVH() = default;
		void Initialize(const CMeshO& mesh);

		ClosestPointResult QueryClosestPoint(const Point3m& q) const;
		bool Empty() const;

	private:
		void Build();
		void BuildNode(int nodeIndex, std::size_t begin, std::size_t end);
		RSS FitRSS(std::size_t begin, std::size_t end, Point3m& mean);
		void QueryNode(int node_index, const Point3m& query_p, ClosestPointResult& result) const;
		void UpdateClosestPointWithTriangle(int face_index, const Point3m& query_p, ClosestPointResult& result) const;
		Scalarm DistanceToRSSCoreSquared(const Point3m& query_p, const RSS& bound) const;
	private:
		int _root = -1;
		const CMeshO* _mesh = nullptr;
		std::vector<SSVHNode> _nodes;
		std::vector<int> _faceIndices;
		std::vector<CachedTriangle> _triangles;
	};
}

#endif
