#include "UniformGrid.h"
namespace meshoffsetting
{
	UniformGrid::UniformGrid(const Box3m& box, const Scalarm& offset, const int count)
		:_bbox(box), _cell_count(count)
	{
		Scalarm max_side = std::max({ _bbox.DimX(), _bbox.DimY() ,_bbox.DimZ() });
		Scalarm cell_width = (max_side + offset * 2)/ count;
		ComputeDimAndVoxel(cell_width);
	}

	UniformGrid::~UniformGrid()
	{
	}

	void UniformGrid::GetSizeAndVoxel(Point3m& size, Point3m& voxel) const
	{
		size = _size;
		voxel = _voxel;
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

	void UniformGrid::ComputeDimAndVoxel(const Scalarm& cell_width)
	{
		_dim = _bbox.max - _bbox.min;
		_size.X() = std::ceil(_dim.X() / cell_width);
		_size.Y() = std::ceil(_dim.Y() / cell_width);
		_size.Z() = std::ceil(_dim.Z() / cell_width);
		_voxel.X() = cell_width;
		_voxel.Y() = cell_width;
		_voxel.Z() = cell_width;
		_dim = Point3m(_size.X() * _voxel.X(),
			_size.Y() * _voxel.Y(),
			_size.Z() * _voxel.Z());
		_nodes.resize(_size.X() * _size.Y() * _size.Z());
	}
}
