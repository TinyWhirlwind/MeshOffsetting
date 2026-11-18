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
#include <array>
#include <cstddef>
#include <vector>
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
		Scalarm signed_distance = 0;
		/*Scalarm unSignedDistance = 0;
		Point3m position = Point3m(0, 0, 0);
		int id = -1;*/
	};

	class UniformGrid
	{
	public:
		static constexpr int kBrickSize = 8;
		static constexpr std::size_t kBrickNodeCount =
			kBrickSize * kBrickSize * kBrickSize;

		UniformGrid(const Box3m& box, const Scalarm& offset, int block_size = 8, const int count = 512);
		~UniformGrid() = default;
		UniformGrid() = default;

		Point3i GetSize() const;
		Scalarm GetCellWidth() const;
		Box3m GetBox() const;

		void ResetNodes();
		void ActivateNodeRange(const Point3i& min_node, const Point3i& max_node);
		void DeactivateNodeRange(const Point3i& min_node, const Point3i& max_node);
		bool TryGetNodeDistance(const Point3i& node, Scalarm& distance) const;
		bool IsInside(const Point3i& node) const;

		std::size_t StoredBrickCount() const;
		GridNode& StoredNode(std::size_t brick_index, std::size_t local_index);
		const GridNode& StoredNode(std::size_t brick_index, std::size_t local_index) const;
		Point3i StoredNodeCoord(std::size_t brick_index, std::size_t local_index) const;

		std::size_t NodeIndex(const Point3i& node) const;
		Point3i NodeCoord(std::size_t index) const;
		Point3m NodePosition(const Point3i& node) const;
		Point3m NodePosition(std::size_t index) const;
	private:
		struct NodeBrick
		{
			Point3i brick_coord = Point3i(0, 0, 0);
			std::array<GridNode, kBrickNodeCount> nodes;
		};

		Box3m _bbox;
		Point3m _dim;
		Point3i _size;
		Point3i _brick_grid_size;
		Point3m _voxel;
		int _cell_count = 0;
		int _block_size = 0;	
		Scalarm _cell_width = 0;
		std::vector<int> _brick_lookup;
		std::vector<NodeBrick> _bricks;

		void ComputeDimAndVoxel();
		std::size_t BrickLookupIndex(const Point3i& brick_coord) const;
		std::size_t LocalNodeIndex(const Point3i& node) const;
		int FindBrick(const Point3i& brick_coord) const;
		int GetOrCreateBrick(const Point3i& brick_coord);
	};
}
#endif
