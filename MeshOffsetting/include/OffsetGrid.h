#pragma once
#include "mymesh.h"
#include <base_type.h>
#include <set>
#include <vcg/space/index/grid_util.h>
using namespace vcg;
template <class ScalarType>
class OffsetGrid : public BasicGrid<ScalarType>
{
public:
    enum NodeState
    {
        NodeUnknown = 0,
        NodeDistanceReady,
        NodeInvalid

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

    struct GridEdge
    {
        int _p0, _p1;
        bool isV;
        GridEdge(int p0, int p1) :
            _p0(std::min(p0, p1)), _p1(std::max(p0, p1)) 
        {
        };
        
        /*void SetV()
        {
            isV = true;
        }
        bool IsV()
        {
            return isV;
        }
        void ClearV()
        {
            isV = false;
        }*/
    };

    struct BlockData
    {
        BlockState state;
        Point3i minCell;
        Point3i cellSpan;
        ScalarType radius;
        Point3m center;
        BlockData() :
            state(BlockState::BlockUnknown),
            cellSpan(0,0,0),
            minCell(0,0,0),
            center(0, 0, 0),
            radius(0)
        {
        }

       BlockData(ScalarType rsq) :
            state(BlockState::BlockUnknown),
            cellSpan(0,0,0),
            minCell(0,0,0),
           radius(rsq)
        {
        }
    };

    OffsetGrid(const Box3m& bbox, ScalarType cell_width, ScalarType offset_value);
    ~OffsetGrid();
    //------------------------------Grid-------------------------------------------
    ScalarType OffsetValue() const;
    ScalarType CellWidth() const;
    Point3i NodeDims() const;
    size_t NodeCount() const;
    std::vector<NodeData>& Nodes();
    const std::vector<NodeData>& Nodes() const;
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

    void SetNodeState(const Point3i& p, NodeState s);
    NodeState GetNodeState(const Point3i& p) const;
    void SetNodeInvalid(const Point3i& p);
    bool IsNodeInvalid(const Point3i& p) const;

    template<class Func>
    void ForEachNode(Func fn) const
    {
        const Point3i dim = NodeDims();
        for (int i = 0; i < dim[0]; ++i)
        {
            for (int j = 0; j < dim[1]; ++j)
            {
                for (int k = 0; k < dim[2]; ++k)
                {
                    fn(Point3i(i, j, k));
                }
            }
        }
    }
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
    ScalarType BlockCircumsphereRadius(const BlockData& block) const;

    void SetBlockState(const Point3i& p, BlockState s);
    BlockState GetBlockState(const Point3i& p) const;
    bool IsBlockInvalid(const Point3i& p) const;
    bool IsBlockRetained(const Point3i& p) const;

    template<class Func>
    void ForEachBlock(Func fn) const
    {
        const Point3i dim = BlockDims();
        for (int i = 0; i < dim[0]; ++i)
        {
            for (int j = 0; j < dim[1]; ++j)
            {
                for (int k = 0; k < dim[2]; ++k)
                {
                    fn(Point3i(i, j, k));
                }
            }
        }
    }

    template<class Func>
    void ForEachNodeInBlock(const Point3i& bp, Func fn) const
    {
        const BlockData& b = Block(bp);

        for (int i = b.minCell[0]; i <= b.minCell[0] + b.cellSpan[0]; ++i)
        {
            for (int j = b.minCell[1]; j <= b.minCell[1] + b.cellSpan[1]; ++j)
            {
                for (int k = b.minCell[2]; k <= b.minCell[2] + b.cellSpan[2]; ++k)
                {
                    fn(Point3i(i, j, k));
                }
            }
        }
    }
    //----------------------------Grid Edge-----------------------------------
    void ClearIntersectEdge();
    const std::set<GridEdge>&  GetIntersectEdge() const;

private:
    int EstimateBlockCount(const ScalarType cell_width, const ScalarType offset_value) const;
    void ForEachAdjacentEdge(const Point3i& start);
private:
    ScalarType _offsetValue;
    ScalarType _cellWidth;
    std::vector<NodeData> _nodes;

    int _blockSide;
    size_t _blockCount;
    Point3i _blockDim;
    std::vector<BlockData> _blocks;
    std::set<GridEdge> _gridEdges;
    
};
