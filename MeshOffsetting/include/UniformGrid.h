/**
 * @file UniformGrid.h
 * @author zhangpingping
 * @date 2026-11-06
 * @brief
 *
 */
#ifndef MESHOFFSETTING_UNIFORM_GRID_H
#define MESHOFFSETTING_UNIFORM_GRID_H
#include "CMeshO.h"
#include <vcg/space/index/grid_util.h>
namespace meshoffsetting
{
	enum class GridNodeState
	{
		Unknown,
		NodeDistanceReady,
		NodeInvalid
	};

	struct GridNode
	{
		GridNodeState state = GridNodeState::Unknown;
		Scalarm signedDistance = 0;
		/*Scalarm unSignedDistance = 0;
		Point3m position = Point3m(0, 0, 0);
		int id = -1;*/
	};

	class UniformGrid
	{
	public:
		UniformGrid(const Box3m& box, const Scalarm& offset, const int count = 512);
		~UniformGrid();

		void GetSizeAndVoxel(Point3m& size, Point3m& voxel) const;
		GridNode& Node(const Point3i& node);
		std::size_t NodeIndex(const Point3i& node);
		Point3i NodeCoord(std::size_t index);
		Point3m NodePosition(const Point3i& node);
	private:
		Box3m _bbox;
		Point3m _dim;
		Point3i _size;
		Point3m _voxel;
		int _cell_count = 0;
		std::vector<GridNode> _nodes;

		void ComputeDimAndVoxel(const Scalarm& cell_width);
	};
}
#endif