#include "UniformGrid.h"
namespace meshoffsetting
{
	UniformGrid::UniformGrid(const Box3m& box, const Scalarm& offset, int block_size, const int count)
		:_bbox(box), _cell_count(count), _block_size(block_size)
	{
		_bbox.Offset(std::abs(offset));
		Scalarm max_side = std::max({ _bbox.DimX(), _bbox.DimY() ,_bbox.DimZ() });
		_cell_width = max_side / (count - 3); //count个节点count-1个cell 左右两边再预留一个cell 模型占count-3
		_bbox.Offset(_cell_width);
		ComputeDimAndVoxel();
	}

	UniformGrid::~UniformGrid()
	{
	}

	void UniformGrid::GetSizeAndVoxel(Point3m& size, Point3m& voxel) const
	{
		size = _size;
		voxel = _voxel;
	}

	Scalarm UniformGrid::GetCellWidth() const
	{
		return _cell_width;
	}

	GridNode& UniformGrid::Node(const Point3i& node)
	{
		size_t index = node.X() + node.Y() * _size.X() + node.Z() * _size.X() * _size.Y();
		return _nodes[index];
	}

	std::size_t UniformGrid::NodeIndex(const Point3i& node)
	{
		return node.X() + node.Y() * _size.X() + node.Z() * _size.X() * _size.Y();
	}

	Point3i UniformGrid::NodeCoord(std::size_t index)
	{
		int x = index % _size.X();
		int y = (index / _size.X()) % _size.Y();
		int z = index / (_size.X() * _size.Y());
		return Point3i(x, y, z);
	}

	Point3m UniformGrid::NodePosition(const Point3i& node)
	{
		Point3m pos;
		pos.X() = _bbox.min.X() + node.X() * _voxel.X();
		pos.Y() = _bbox.min.Y() + node.Y() * _voxel.Y();
		pos.Z() = _bbox.min.Z() + node.Z() * _voxel.Z();
		return pos;
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
		_dim = Point3m{ _voxel.X() * count_x, _voxel.Y() * count_y, _voxel.Z() * count_z };
		_nodes.resize(_size.X() * _size.Y() * _size.Z());
	}
}
