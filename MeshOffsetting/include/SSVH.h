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
	enum ClosestType
	{
		Unknown = 0,
		Face,//0
		Edge,//01,12,02
		Vertex//0,1,2
	};

	struct ClosestPointResult
	{
		bool found = false;
		int faceId = -1;
		Point3m closest_point;
		Scalarm unsigned_distance = std::numeric_limits<Scalarm>::max();
		Scalarm u = 0.0;
		Scalarm	v = 0.0;

		ClosestType closestType = ClosestType::Unknown;
		int closestFeature = -1;
		bool intersected = false;
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
		~SSVH() = default;

		ClosestPointResult QueryClosestPoint(const Point3m& q);
		bool Empty() const;

	private:
		void Build();
		void BuildNode(int nodeIndex, std::size_t begin, std::size_t end);
		RSS FitRSS(std::size_t begin, std::size_t end, Point3m& mean);
		Point3m FaceCentroid(int fid);
		void QueryNode(int node_index, const Point3m& query_p, ClosestPointResult& result);
		void UpdateClosestPointWithTriangle(int face_index, const Point3m& query_p, ClosestPointResult& result);
		ClosestPointResult CalcDistancePointToPrimitive(const Point3m& p, int face_index) const;
		Scalarm DistanceToRSS(const Point3m& query_p, const RSS& bound) const;
	private:
		int _root = -1;
		CMeshO _mesh;
		std::vector<SSVHNode> _nodes;
		std::vector<int> _faceIndices;
	};
}

#endif