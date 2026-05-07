#include "mymesh.h"
#include <base_type.h>
#include <vcg/space/index/grid_util.h>
using namespace vcg;
template <class ScalarType>
class OffsetGrid : public BasicGrid<ScalarType>
{
public:
    enum NodeState
    {
        NodeUnknown = 0,
        NodeDistanceReady

        //NodeInvalid,
        //BlockRetained
    };
    enum BlockState
    {
        BlockUnknown = 0,
        BlockRetained,
        BlockInvalid
    };

    struct NodeData
    {
        NodeState state;
        ScalarType unsignedDistance;
        ScalarType signedDistance;
        bool hasUnsignedDistance;
        bool hasSignedDistance;

        NodeData() :state(NodeState::NodeUnknown),
            unsignedDistance(0), 
            signedDistance(0), 
            hasUnsignedDistance(false), 
            hasSignedDistance(false) 
        {}
    };
    struct BlockData
    {
        BlockState state;
        Point3i minCell;
        Point3i cellSpan;
        ScalarType radiusSq;
        Point3m center;
        BlockData() :
            state(BlockState::BlockUnknown),
            cellSpan(0,0,0),
            minCell(0,0,0),
            center(0, 0, 0),
            radiusSq(0)
        {
        }

       BlockData(ScalarType rsq) :
            state(BlockState::BlockUnknown),
            cellSpan(0,0,0),
            minCell(0,0,0),
           radiusSq(rsq)
        {
        }
    };

    OffsetGrid(ScalarType offset_value, ScalarType cell_width);
    ~OffsetGrid();
    //------------------------------Grid-------------------------------------------
    void Init(const Box3m& bbox, ScalarType cell_width, ScalarType offset_value);
    void Init(const Box3m& bbox, const Point3i cell_count, ScalarType offset_value);
    ScalarType OffsetValue() const;
    ScalarType CellWidth() const;
    Point3i NodeDims() const;
    size_t NodeCount() const;
    bool IsVaildNodeCoord(const Point3i& p) const;
    size_t NodeIndex(const Point3i& p) const;
    NodeData& Node(const Point3i& p);
    const NodeData& Node(const Point3i& p) const;
    void ResetNodes();
    // node 整数坐标 -> 世界坐标。这个位置就是 SSVH 最近距离查询里的 q。
    Point3m NodePosition(const Point3i& p) const;
    // 缓存一次 unsigned distance 查询结果
    void SetUnsignedDistance(const Point3i& p, ScalarType d);
    // 缓存一次 signed distance 查询结果
    void SetSignedDistance(const Point3i& p, ScalarType d);
    bool HasUnsignedDistance(const Point3i& p) const;
    bool HasSignedDistance(const Point3i& p) const;
    ScalarType UnsignedDistance(const Point3i& p) const;
    ScalarType SignedDistance(const Point3i& p) const;
    //判断 node 是否在 offset 窄带附近
    bool IsNodeInNarrowBand(const Point3i& p) const;
    template<class Func>
    void ForEachNode(Func fn) const;

    //------------------------------Block-------------------------------------------
    void BuildBlocks();
    void ResetBlocks();

    bool IsVaildBlockCoord(const Point3i& p) const; 
    size_t BlockIndex(const Point3i& p) const;

    BlockData& Block(const Point3i& p);
    const BlockData& Block(const Point3i& p) const;

    int BlockSide() const;
    Point3i BlockDims() const;
    size_t BlockCount() const;

    Point3m BlockCenter(const BlockData& block);
    ScalarType BlockCircumsphereRadiusSq(const BlockData& block) const;

    void SetBlockState(const Point3i& p, BlockState s);
    BlockState GetBlockState(const Point3i& p) const;
    bool IsBlockInvalid(const Point3i& p) const;
    bool IsBlockRetained(const Point3i& p) const;

    template<class Func>
    void ForEachNodeInBlock(const Point3i& bp, Func fn) const;

private:
    int EstimateBlockCount(const ScalarType cell_width, const ScalarType offset_value) const;

private:
    ScalarType _offsetValue;
    ScalarType _cellWidth;
    std::vector<NodeData> _nodes;

    int _blockSide;
    size_t _blockCount;
    Point3i _blockDim;
    std::vector<BlockData> _blocks;
    
};