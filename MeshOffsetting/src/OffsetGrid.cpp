#include "OffsetGrid.h"
#include <algorithm>
#include <math.h>
template <class ScalarType>
OffsetGrid<ScalarType>::OffsetGrid(const Box3m& bbox, ScalarType cell_width, ScalarType offset_value)
{
    assert(cell_width > ScalarType(0));

    Point3i cell_count(1, 1, 1);
    Box3m real_bbox = bbox;
    const Point3<ScalarType> size = bbox.max - bbox.min;
    for (int i = 0; i < 3; ++i)
    {
        cell_count[i] = std::max(1, int(std::ceil(size[i] / cell_width)));
        real_bbox.max[i] = real_bbox.min[i] + ScalarType(cell_count[i]) * cell_width;
    }
    assert(cell_count[0] > 0);
    assert(cell_count[1] > 0);
    assert(cell_count[2] > 0);
    this->bbox.max = real_bbox.max;
    this->bbox.min = real_bbox.min;
    this->siz = cell_count;
    this->ComputeDimAndVoxel();
    _offsetValue = offset_value;
    _cellWidth = this->voxel[0];
    _nodes.assign(NodeCount(), NodeData());
}

template <class ScalarType>
OffsetGrid<ScalarType>::~OffsetGrid() {};

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
std::vector<typename OffsetGrid<ScalarType>::NodeData>& OffsetGrid<ScalarType>::Nodes()
{
    return _nodes;
}

template <class ScalarType>
const std::vector<typename OffsetGrid<ScalarType>::NodeData>& OffsetGrid<ScalarType>::Nodes() const
{
    return _nodes;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsVaildNodeCoord(const Point3i& p) const
{
    const Point3i dim = NodeDims();
    return p[0] >= 0 && p[0] < dim[0] &&
           p[1] >= 0 && p[1] < dim[1] &&
        p[2] >= 0 && p[2] < dim[2];
}

template <class ScalarType>
size_t OffsetGrid<ScalarType>::NodeIndex(const Point3i& p) const
{
    assert(IsVaildNodeCoord(p));
    const Point3i dim = NodeDims();//x变化最快，y次之，z变化最慢，
    //return size_t(dim[0] + p->siz[0] * dim[1] + this->siz[0] * this->siz[1] * dim[2]);
    return size_t(p[0] + p[1] * dim[0] + p[2] * dim[0] * dim[1]);
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
typename OffsetGrid<ScalarType>::NodeData& OffsetGrid<ScalarType>::Node(const int& id)
{
    return _nodes[id];
}

template <class ScalarType>
const typename OffsetGrid<ScalarType>::NodeData& OffsetGrid<ScalarType>::Node(const int& id) const
{
    return _nodes[id];
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
        this->bbox.min[0] + p[0] * this->voxel[0],
        this->bbox.min[1] + p[1] * this->voxel[1],
        this->bbox.min[2] + p[2] * this->voxel[2]
    );
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetUnsignedDistance(const Point3i& p, QueryResult qr)
{
    NodeData& node = Node(p);
    node.unsignedDistance = qr._dist;
    //node.state = NodeState::NodeDistanceReady;
    node.hasUnsignedDistance = true;
    //node.triId = tid;
    node.query = qr;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetSignedDistance(const Point3i& p, QueryResult qr)
{
    NodeData& node = Node(p);
    node.signedDistance = static_cast<float>(qr._sign) * static_cast<float>(qr._dist);
    node.hasSignedDistance = true;
    node.state = NodeState::NodeDistanceReady;
    //node.triId = qr.id;
    node.query = qr;
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
void OffsetGrid<ScalarType>::SetNodeState(const Point3i& p, NodeState s)
{
    Node(p).state = s;
}

template <class ScalarType>
typename OffsetGrid<ScalarType>::NodeState OffsetGrid<ScalarType>::GetNodeState(const Point3i& p) const
{
    return Node(p).state;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetNodeInvalid(const Point3i& p)
{
    NodeData& node = Node(p);
    node.state = NodeState::NodeInvalid;
    node.hasUnsignedDistance = false;
    node.hasSignedDistance = false;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsNodeInvalid(const Point3i& p) const
{
    return GetNodeState(p) == NodeState::NodeInvalid;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsNodeInNarrowBand(const Point3i& p) const
{
    const NodeData& n = Node(p);
    assert(n.hasSignedDistance);
    return std::abs(n.signedDistance - _offsetValue) <= _cellWidth;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::BuildBlocks()
{
    _blockSide = EstimateBlockCount(_cellWidth, _offsetValue);
    _blockDim[0] = std::max(1, (this->siz[0] + _blockSide - 1) / _blockSide);
    _blockDim[1] = std::max(1, (this->siz[1] + _blockSide - 1) / _blockSide);
    _blockDim[2] = std::max(1, (this->siz[2] + _blockSide - 1) / _blockSide);
    _blockCount = static_cast<size_t>(_blockDim[0]) * _blockDim[1] * _blockDim[2];
    _blocks.clear();
    _blocks.resize(_blockCount);
    ScalarType radius = std::sqrt(3.0) * _blockSide * _cellWidth * 0.5;
    _blocks.assign(_blockCount, BlockData(radius));

    for (int i = 0; i < _blockDim[0]; ++i)
    {
        for (int j = 0; j < _blockDim[1]; ++j)
        {
            for (int k = 0; k < _blockDim[2]; ++k)
            {
                BlockData& b = Block(Point3i(i, j, k));
                b.minCell = Point3i(i * _blockSide, j * _blockSide, k * _blockSide);
                b.cellSpan = Point3i(_blockSide, _blockSide, _blockSide);
                Point3i maxCell;
                maxCell[0] = b.minCell[0] + _blockSide;
                maxCell[1] = b.minCell[1] + _blockSide;
                maxCell[2] = b.minCell[2] + _blockSide;
                bool isCalcRadius = false;
                if (maxCell[0] > this->siz[0])
                {
                    isCalcRadius = true;
                    b.cellSpan[0] = this->siz[0] - b.minCell[0];
                }
                if (maxCell[1] > this->siz[1])
                {
                    isCalcRadius = true;
                    b.cellSpan[1] = this->siz[1] - b.minCell[1];
                }
                if (maxCell[2] > this->siz[2])
                {
                    isCalcRadius = true;
                    b.cellSpan[2] = this->siz[2] - b.minCell[2];
                }

                b.center = Point3m(
                    this->bbox.min[0] + (ScalarType(b.minCell[0]) + ScalarType(b.cellSpan[0]) * ScalarType(0.5)) * _cellWidth,
                    this->bbox.min[1] + (ScalarType(b.minCell[1]) + ScalarType(b.cellSpan[1]) * ScalarType(0.5)) * _cellWidth,
                    this->bbox.min[2] + (ScalarType(b.minCell[2]) + ScalarType(b.cellSpan[2]) * ScalarType(0.5)) * _cellWidth
                );
                if (isCalcRadius)
                {
                    const ScalarType hx = ScalarType(b.cellSpan[0]) * _cellWidth * ScalarType(0.5);
                    const ScalarType hy = ScalarType(b.cellSpan[1]) * _cellWidth * ScalarType(0.5);
                    const ScalarType hz = ScalarType(b.cellSpan[2]) * _cellWidth * ScalarType(0.5);

                    b.radius = std::sqrt(hx * hx + hy * hy + hz * hz);
                }
            }
        }
    }
}

template <class ScalarType>
void OffsetGrid<ScalarType>::ResetBlocks()
{
    std::fill(_blocks.begin(), _blocks.end(), BlockData());
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsVaildBlockCoord(const Point3i& p) const
{
    const Point3i dim = BlockDims();
    return p[0] >= 0 && p[0] < dim[0] &&
        p[1] >= 0 && p[1] < dim[1] &&
        p[2] >= 0 && p[2] < dim[2];
}

template <class ScalarType>
size_t OffsetGrid<ScalarType>::BlockIndex(const Point3i& p) const
{
    assert(IsVaildBlockCoord(p));
    const Point3i dim = BlockDims();//x变化最快，y次之，z变化最慢，
    return size_t(p[0] + p[1] * dim[0] + p[2] * dim[0] * dim[1]);
}

template <class ScalarType>
typename OffsetGrid<ScalarType>::BlockData& OffsetGrid<ScalarType>::Block(const Point3i& p)
{
    return _blocks[BlockIndex(p)];
}

template <class ScalarType>
const typename OffsetGrid<ScalarType>::BlockData& OffsetGrid<ScalarType>::Block(const Point3i& p) const
{
    return _blocks[BlockIndex(p)];
}

template <class ScalarType>
int OffsetGrid<ScalarType>::BlockSide() const
{
    return _blockSide;
}

template <class ScalarType>
Point3i OffsetGrid<ScalarType>::BlockDims() const
{
    return _blockDim;
}

template <class ScalarType>
size_t OffsetGrid<ScalarType>::BlockCount() const
{
    return _blockCount;
}

template <class ScalarType>
Point3m OffsetGrid<ScalarType>::BlockCenter(const BlockData& block)
{
    return block.center;
}

template <class ScalarType>
ScalarType OffsetGrid<ScalarType>::BlockCircumsphereRadius(const BlockData& block) const
{
    return block.radius;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::SetBlockState(const Point3i& p, BlockState s)
{
    BlockData& b = Block(p);
    b.state = s;
}

template <class ScalarType>
typename OffsetGrid<ScalarType>::BlockState OffsetGrid<ScalarType>::GetBlockState(const Point3i& p) const
{
    return Block(p).state;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsBlockInvalid(const Point3i& p) const
{
    return GetBlockState(p) == BlockState::BlockInvalid;
}

template <class ScalarType>
bool OffsetGrid<ScalarType>::IsBlockRetained(const Point3i& p) const
{
    return GetBlockState(p) == BlockState::BlockRetained;
}

template <class ScalarType>
std::set<typename OffsetGrid<ScalarType>::GridEdge>& OffsetGrid<ScalarType>::GetIntersectEdges()
{
    BuildIntersectEdges();
    return &_gridEdges;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::BuildHermiteSamples()
{

}

template <class ScalarType>
void OffsetGrid<ScalarType>::FindIntersectionPointOnGridEdge()
{

}

template <class ScalarType>
int OffsetGrid<ScalarType>::EstimateBlockCount(const ScalarType cell_width, const ScalarType offset_value) const
{
    const ScalarType r = std::abs(offset_value);
    const ScalarType l = std::abs(cell_width);

    if (l <= 0.0) return 2;
    if (r <= 0.0) return 2;

    int m = static_cast<int>(std::ceil(r / l));

    if (m < 2)   m = 2;
    if (m > 128) m = 128;

    return m;
}

template <class ScalarType>
void OffsetGrid<ScalarType>::BuildIntersectEdges()
{
    _gridEdges.clear();
    Point3i nodeDims = NodeDims();
    for (int i = 0; i <= this->siz[0]; ++i)
    {
        for (int j = 0; j <= this->siz[1]; ++j)
        {
            for (int k = 0; k <= this->siz[2]; ++k)
            {
                ForEachAdjacentEdge(Point3i{ i,j,k });
            }
        }
    }
}

template <class ScalarType>
void OffsetGrid<ScalarType>::ForEachAdjacentEdge(const Point3i& start)
{
    int sx = start.X();
    int sy = start.Y();
    int sz = start.Z();
    const NodeData& sn = Node(start);
    Point3i end;
    for (int i = 1; i <= 2; ++i)
    {
        int step = std::pow(-1, i);
        for (int j = 0; j < 3; ++j)
        {
            end[j] = start[j] + step;
            if (end[j] < 0 || end[j] > this->siz[j])continue;
            const NodeData& en = Node(end);
            if (en.state != NodeState::NodeDistanceReady)continue;
            if (!sn.hasSignedDistance || !en.hasSignedDistance)continue;
            if (sn.signedDistance * en.signedDistance >= 0)continue;

            size_t startId = NodeIndex(start);
            size_t endId = NodeIndex(end);
            GridEdge edge(startId, endId);
            _gridEdges.insert(edge);
        }
    }
}

template class OffsetGrid<float>;
template class OffsetGrid<double>;
