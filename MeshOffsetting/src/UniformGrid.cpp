#include "UniformGrid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace meshoffsetting
{
	UniformGrid::UniformGrid(
		const Box3m& box,
		const Scalarm& offset,
		int block_size,
		const int count)
		: _bbox(box), _cell_count(count), _block_size(block_size)
	{
		if (count <= 3)
		{
			throw std::invalid_argument("UniformGrid node count must be greater than 3");
		}
		if (block_size <= 0)
		{
			throw std::invalid_argument("UniformGrid block size must be positive");
		}
		if (!std::isfinite(offset))
		{
			throw std::invalid_argument("UniformGrid offset must be finite");
		}

		_bbox.Offset(std::abs(offset));
		Scalarm dim_x = _bbox.DimX();
		Scalarm dim_y = _bbox.DimY();
		Scalarm dim_z = _bbox.DimZ();
		Scalarm max_side = std::max({ dim_x, dim_y, dim_z });
		if (!std::isfinite(dim_x) || !std::isfinite(dim_y) ||
			!std::isfinite(dim_z) || !std::isfinite(max_side) ||
			max_side <= Scalarm(0))
		{
			throw std::invalid_argument("UniformGrid box must have a finite positive size");
		}

		_cell_width = max_side / (count - 3);
		_bbox.Offset(_cell_width);
		ComputeDimAndVoxel();
	}

	Point3i UniformGrid::GetSize() const
	{
		return _size;
	}

	Scalarm UniformGrid::GetCellWidth() const
	{
		return _cell_width;
	}

	Box3m UniformGrid::GetBox() const
	{
		return _bbox;
	}

	void UniformGrid::ResetNodes()
	{
		_bricks.clear();
		std::fill(_brick_lookup.begin(), _brick_lookup.end(), -1);
	}

	void UniformGrid::ActivateNodeRange(
		const Point3i& min_node,
		const Point3i& max_node)
	{
		int min_x = std::max(0, min_node.X());
		int min_y = std::max(0, min_node.Y());
		int min_z = std::max(0, min_node.Z());
		int max_x = std::min(_size.X() - 1, max_node.X());
		int max_y = std::min(_size.Y() - 1, max_node.Y());
		int max_z = std::min(_size.Z() - 1, max_node.Z());

		if (min_x > max_x || min_y > max_y || min_z > max_z)
		{
			return;
		}

		int min_bx = min_x / kBrickSize;
		int min_by = min_y / kBrickSize;
		int min_bz = min_z / kBrickSize;
		int max_bx = max_x / kBrickSize;
		int max_by = max_y / kBrickSize;
		int max_bz = max_z / kBrickSize;

		for (int bz = min_bz; bz <= max_bz; ++bz)
		{
			for (int by = min_by; by <= max_by; ++by)
			{
				for (int bx = min_bx; bx <= max_bx; ++bx)
				{
					Point3i brick_coord(bx, by, bz);
					int brick_index = GetOrCreateBrick(brick_coord);
					NodeBrick& brick = _bricks[brick_index];

					int begin_x = std::max(min_x, bx * kBrickSize);
					int begin_y = std::max(min_y, by * kBrickSize);
					int begin_z = std::max(min_z, bz * kBrickSize);
					int end_x = std::min(max_x, bx * kBrickSize + kBrickSize - 1);
					int end_y = std::min(max_y, by * kBrickSize + kBrickSize - 1);
					int end_z = std::min(max_z, bz * kBrickSize + kBrickSize - 1);

					for (int z = begin_z; z <= end_z; ++z)
					{
						for (int y = begin_y; y <= end_y; ++y)
						{
							for (int x = begin_x; x <= end_x; ++x)
							{
								GridNode& node = brick.nodes[
									LocalNodeIndex(Point3i(x, y, z))];
								if (node.state == GridNodeState::NodeInvalid)
								{
									node.state = GridNodeState::Unknown;
									node.signed_distance = Scalarm(0);
								}
							}
						}
					}
				}
			}
		}
	}

	void UniformGrid::DeactivateNodeRange(
		const Point3i& min_node,
		const Point3i& max_node)
	{
		int min_x = std::max(0, min_node.X());
		int min_y = std::max(0, min_node.Y());
		int min_z = std::max(0, min_node.Z());
		int max_x = std::min(_size.X() - 1, max_node.X());
		int max_y = std::min(_size.Y() - 1, max_node.Y());
		int max_z = std::min(_size.Z() - 1, max_node.Z());

		if (min_x > max_x || min_y > max_y || min_z > max_z)
		{
			return;
		}

		int min_bx = min_x / kBrickSize;
		int min_by = min_y / kBrickSize;
		int min_bz = min_z / kBrickSize;
		int max_bx = max_x / kBrickSize;
		int max_by = max_y / kBrickSize;
		int max_bz = max_z / kBrickSize;

		for (int bz = min_bz; bz <= max_bz; ++bz)
		{
			for (int by = min_by; by <= max_by; ++by)
			{
				for (int bx = min_bx; bx <= max_bx; ++bx)
				{
					int brick_index = FindBrick(Point3i(bx, by, bz));
					if (brick_index < 0)
					{
						continue;
					}

					NodeBrick& brick = _bricks[brick_index];
					int begin_x = std::max(min_x, bx * kBrickSize);
					int begin_y = std::max(min_y, by * kBrickSize);
					int begin_z = std::max(min_z, bz * kBrickSize);
					int end_x = std::min(max_x, bx * kBrickSize + kBrickSize - 1);
					int end_y = std::min(max_y, by * kBrickSize + kBrickSize - 1);
					int end_z = std::min(max_z, bz * kBrickSize + kBrickSize - 1);

					for (int z = begin_z; z <= end_z; ++z)
					{
						for (int y = begin_y; y <= end_y; ++y)
						{
							for (int x = begin_x; x <= end_x; ++x)
							{
								GridNode& node = brick.nodes[
									LocalNodeIndex(Point3i(x, y, z))];
								node.state = GridNodeState::NodeInvalid;
								node.signed_distance = Scalarm(0);
							}
						}
					}
				}
			}
		}
	}

	bool UniformGrid::TryGetNodeDistance(
		const Point3i& node,
		Scalarm& distance) const
	{
		if (!IsInside(node))
		{
			return false;
		}

		Point3i brick_coord(
			node.X() / kBrickSize,
			node.Y() / kBrickSize,
			node.Z() / kBrickSize);
		int brick_index = FindBrick(brick_coord);
		if (brick_index < 0)
		{
			return false;
		}

		const GridNode& stored_node =
			_bricks[brick_index].nodes[LocalNodeIndex(node)];
		if (stored_node.state != GridNodeState::NodeDistanceReady ||
			!std::isfinite(stored_node.signed_distance))
		{
			return false;
		}

		distance = stored_node.signed_distance;
		return true;
	}

	bool UniformGrid::IsInside(const Point3i& node) const
	{
		return node.X() >= 0 && node.X() < _size.X() &&
			node.Y() >= 0 && node.Y() < _size.Y() &&
			node.Z() >= 0 && node.Z() < _size.Z();
	}

	std::size_t UniformGrid::StoredBrickCount() const
	{
		return _bricks.size();
	}

	GridNode& UniformGrid::StoredNode(
		std::size_t brick_index,
		std::size_t local_index)
	{
		return _bricks[brick_index].nodes[local_index];
	}

	const GridNode& UniformGrid::StoredNode(
		std::size_t brick_index,
		std::size_t local_index) const
	{
		return _bricks[brick_index].nodes[local_index];
	}

	Point3i UniformGrid::StoredNodeCoord(
		std::size_t brick_index,
		std::size_t local_index) const
	{
		const Point3i& brick_coord = _bricks[brick_index].brick_coord;
		int local_x = static_cast<int>(local_index % kBrickSize);
		int local_y = static_cast<int>((local_index / kBrickSize) % kBrickSize);
		int local_z = static_cast<int>(
			local_index / (kBrickSize * kBrickSize));

		return Point3i(
			brick_coord.X() * kBrickSize + local_x,
			brick_coord.Y() * kBrickSize + local_y,
			brick_coord.Z() * kBrickSize + local_z);
	}

	std::size_t UniformGrid::NodeIndex(const Point3i& node) const
	{
		return static_cast<std::size_t>(node.X()) +
			static_cast<std::size_t>(node.Y()) *
			static_cast<std::size_t>(_size.X()) +
			static_cast<std::size_t>(node.Z()) *
			static_cast<std::size_t>(_size.X()) *
			static_cast<std::size_t>(_size.Y());
	}

	Point3i UniformGrid::NodeCoord(std::size_t index) const
	{
		int x = static_cast<int>(index % static_cast<std::size_t>(_size.X()));
		int y = static_cast<int>(
			(index / static_cast<std::size_t>(_size.X())) %
			static_cast<std::size_t>(_size.Y()));
		int z = static_cast<int>(
			index /
			(static_cast<std::size_t>(_size.X()) *
			 static_cast<std::size_t>(_size.Y())));
		return Point3i(x, y, z);
	}

	Point3m UniformGrid::NodePosition(const Point3i& node) const
	{
		Point3m pos;
		pos.X() = _bbox.min.X() + node.X() * _voxel.X();
		pos.Y() = _bbox.min.Y() + node.Y() * _voxel.Y();
		pos.Z() = _bbox.min.Z() + node.Z() * _voxel.Z();
		return pos;
	}

	Point3m UniformGrid::NodePosition(std::size_t index) const
	{
		return NodePosition(NodeCoord(index));
	}

	void UniformGrid::ComputeDimAndVoxel()
	{
		auto paddedCount = [this](Scalarm dim, Scalarm cell_width) -> int
		{
			int count = static_cast<int>(std::ceil(dim / cell_width));
			count = ((count + _block_size - 1) / _block_size) * _block_size;
			return count;
		};

		int count_x = paddedCount(_bbox.DimX(), _cell_width);
		int count_y = paddedCount(_bbox.DimY(), _cell_width);
		int count_z = paddedCount(_bbox.DimZ(), _cell_width);
		_size = Point3i(count_x + 1, count_y + 1, count_z + 1);
		_voxel = Point3m(_cell_width, _cell_width, _cell_width);
		_dim = Point3m{
			_voxel.X() * count_x,
			_voxel.Y() * count_y,
			_voxel.Z() * count_z };
		_bbox.max.X() = _bbox.min.X() + _dim.X();
		_bbox.max.Y() = _bbox.min.Y() + _dim.Y();
		_bbox.max.Z() = _bbox.min.Z() + _dim.Z();

		_brick_grid_size = Point3i(
			(_size.X() + kBrickSize - 1) / kBrickSize,
			(_size.Y() + kBrickSize - 1) / kBrickSize,
			(_size.Z() + kBrickSize - 1) / kBrickSize);

		std::size_t lookup_size =
			static_cast<std::size_t>(_brick_grid_size.X()) *
			static_cast<std::size_t>(_brick_grid_size.Y()) *
			static_cast<std::size_t>(_brick_grid_size.Z());
		_brick_lookup.assign(lookup_size, -1);
		_bricks.clear();
	}

	std::size_t UniformGrid::BrickLookupIndex(
		const Point3i& brick_coord) const
	{
		return static_cast<std::size_t>(brick_coord.X()) +
			static_cast<std::size_t>(brick_coord.Y()) *
			static_cast<std::size_t>(_brick_grid_size.X()) +
			static_cast<std::size_t>(brick_coord.Z()) *
			static_cast<std::size_t>(_brick_grid_size.X()) *
			static_cast<std::size_t>(_brick_grid_size.Y());
	}

	std::size_t UniformGrid::LocalNodeIndex(const Point3i& node) const
	{
		int local_x = node.X() % kBrickSize;
		int local_y = node.Y() % kBrickSize;
		int local_z = node.Z() % kBrickSize;
		return static_cast<std::size_t>(
			local_x +
			local_y * kBrickSize +
			local_z * kBrickSize * kBrickSize);
	}

	int UniformGrid::FindBrick(const Point3i& brick_coord) const
	{
		if (brick_coord.X() < 0 || brick_coord.X() >= _brick_grid_size.X() ||
			brick_coord.Y() < 0 || brick_coord.Y() >= _brick_grid_size.Y() ||
			brick_coord.Z() < 0 || brick_coord.Z() >= _brick_grid_size.Z())
		{
			return -1;
		}

		return _brick_lookup[BrickLookupIndex(brick_coord)];
	}

	int UniformGrid::GetOrCreateBrick(const Point3i& brick_coord)
	{
		std::size_t lookup_index = BrickLookupIndex(brick_coord);
		int brick_index = _brick_lookup[lookup_index];
		if (brick_index >= 0)
		{
			return brick_index;
		}

		NodeBrick brick;
		brick.brick_coord = brick_coord;
		for (GridNode& node : brick.nodes)
		{
			node.state = GridNodeState::NodeInvalid;
			node.signed_distance = Scalarm(0);
		}

		brick_index = static_cast<int>(_bricks.size());
		_bricks.push_back(std::move(brick));
		_brick_lookup[lookup_index] = brick_index;
		return brick_index;
	}
}
