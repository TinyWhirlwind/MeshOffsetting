#include "OffsetGrid.h"
#include <algorithm>
template <class ScalarType>
OffsetGrid<ScalarType>::OffsetGrid(ScalarType offset_value, ScalarType cell_width) :
    _offsetValue(offset_value), _cellWidth(cell_width)
{
}

template <class ScalarType>
OffsetGrid<ScalarType>::~OffsetGrid() {};

template <class ScalarType>
void OffsetGrid<ScalarType>::Init(const Box3m& bbox, ScalarType cell_width, ScalarType offset_value)
{
    assert(cell_width > ScalarType(0));

    Point3i cell_count(1,1,1);
    Box3m real_bbox = bbox;
    const Point3<ScalarType> size = bbox.max - bbox.min;
    for (int i = 0; i < 3; ++i)
    {
        cell_count[i] = std::max(1, int(std::ceil(size[i] / cell_width)));
        real_bbox.max[i] = real_bbox.min[i] + ScalarType(cell_count[i]) * cell_width;
    }
    Init(real_bbox, cell_count, offset_value);
}

template <class ScalarType>
void OffsetGrid<ScalarType>::Init(const Box3m& bbox, const Point3i cell_count, ScalarType offset_value)
{
    assert(cell_count[0] > 0);
    assert(cell_count[1] > 0);
    assert(cell_count[2] > 0);

    this->bbox = bbox;
    this->siz = cell_count;
    this->ComputeDimAndVoxel();
    _offsetValue = offset_value;
    _cellWidth = this->voxel;
    _nodes.assign(NodeCount(), NodeData());
}

template <class ScalarType>
ScalarType OffsetGrid<ScalarType>::OffsetValue() const
{
    return _offsetValue;
}

template <class ScalarType>
ScalarType OffsetGrid<ScalarType>::CellWidth() const
{
    return _cellWidth;
}

template <class ScalarType>
Point3i OffsetGrid<ScalarType>::NodeDims() const
{
    return Point3i(this->siz[0] + 1, this->siz[1] + 1, this->siz[2] + 1);
}

template <class ScalarType>
size_t OffsetGrid<ScalarType>::NodeCount() const
{
    Point3i dim = NodeDims();
    return size_t(dim[0] * dim[1] * dim[2]);
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsValidNodeCoord(const Point3i& p) const
{
    const Point3i dim = NodeDims();
    return p[0] >= 0 && p[0] < dim[0] &&
           p[1] >= 0 && p[1] < dim[1] &&
        p[2] >= 0 && p[2] < dim[2];
}

template <class ScalarType>
size_t OffsetGrid<ScalarType>::NodeIndex(const Point3i& p) const
{
    assert(IsValidNodeCoord(p));
    const Point3i dim = NodeDims();//x变化最快，y次之，z变化最慢，
    return size_t(dim[0] + this->siz[0] * dim[1] + this->siz[0] * this->siz[1] * dim[2]);
}


template <class ScalarType>
typename OffsetGrid<ScalarType>::NodeData& OffsetGrid<ScalarType>::Node(const Point3i& p)
{
    return _nodes[NodeIndex(p)];
}

template <class ScalarType>
const typename OffsetGrid<ScalarType>::NodeData& OffsetGrid<ScalarType>::Node(const Point3i& p) const
{
    return _nodes[NodeIndex(p)];
}

template <class ScalarType>
void OffsetGrid<ScalarType>::ResetNodes()
{
    std::fill(_nodes.begin(), _nodes.end(), NodeData());
}

template <class ScalarType>
Point3m OffsetGrid<ScalarType>::NodePosition(const Point3i& p) const
{
    return Point3m(
        this->bbox.min[0] + p[0] * this->voxel,
        this->bbox.min[1] + p[1] * this->voxel,
        this->bbox.min[2] + p[2] * this->voxel
    );
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetUnsignedDistance(const Point3i& p, ScalarType d)
{
    NodeData node = Node(p);
    node.unsignedDistanced = d;
    //node.state = NodeState::NodeDistanceReady;
    node.hasUnsignedDistance = true;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetSignedDistance(const Point3i& p, ScalarType d)
{
    NodeData node = Node(p);
    node.signedDistanced = d;
    node.hasSignedDistance = true;
    node.state = NodeState::NodeDistanceReady;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::HasUnsignedDistance(const Point3i& p) const
{
    return Node(p).hasUnsignedDistance;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::HasSignedDistance(const Point3i& p) const
{
    return Node(p).hasSignedDistance;
}

template <class ScalarType>
ScalarType OffsetGrid<ScalarType>::UnsignedDistance(const Point3i& p) const
{
    return Node(p).unsignedDistance;
}

template <class ScalarType>
ScalarType OffsetGrid<ScalarType>::SignedDistance(const Point3i& p) const
{
    return Node(p).signedDistance;
}

template <class ScalarType>
//判断 node 是否在 offset 窄带附近
bool OffsetGrid<ScalarType>::IsNodeInNarrowBand(const Point3i& p) const
{
    const NodeData& n = Node(p);
    assert(n.hasSignedDistance);
    return std::abs(n.signedDistance - _offsetValue) <= _cellWidth;
}

template <class ScalarType>
template<class Func>
void OffsetGrid<ScalarType>::ForEachNode(Func fn) const
{
    for (int i = 0; i < this->siz[0]; ++i)
    {
        for (int j = 0; j < this->siz[1]; ++j)
        {
            for (int k = 0; k < this->siz[2]; ++k)
            {
                fn(Point3i(i, j, k));
            }
        }
    }
}
