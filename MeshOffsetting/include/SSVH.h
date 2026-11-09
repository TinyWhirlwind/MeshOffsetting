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
#include "SignedDistance.h"
#include <vector>
namespace meshoffsetting
{
	struct ClosestPointResult
	{
		bool found = false;
		int faceId = -1;
		Point3m closestPoint;
		Scalarm distance = std::numeric_limits<Scalarm>::max();
		Scalarm u, v = 0.0;
	};
	
	struct RSS
	{
		Scalarm sweepRadius = 0;
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

		void Build(const CMeshO& mesh);
		ClosestPointResult QueryClosestPoint(const Point3m& q);
		bool Empty() const;

	private:
		void BuildNode(int nodeIndex, std::size_t begin, std::size_t end);
		RSS FitRSS(std::size_t begin, std::size_t end, Point3m mean);
		Point3m FaceCentroid(int fid);
	private:
		int _root = 0;
		CMeshO _mesh;
		std::vector<SSVHNode> _nodes;
		std::vector<int> _faceIndices;
	};
}

#endif