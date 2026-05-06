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
            //NodeValid
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
            unsignedDistance(false) 
        {}
    };

    OffsetGrid(ScalarType offset_value, ScalarType cell_width);
    ~OffsetGrid();

    void Init(const Box3m& bbox, ScalarType cell_width, ScalarType offset_value);
    void Init(const Box3m& bbox, const Point3i cell_count, ScalarType offset_value);

    ScalarType OffsetValue() const;
    ScalarType CellWidth() const;
    Point3i NodeDims() const;
    size_t NodeCount() const;
    bool IsValidNodeCoord(const Point3i& p) const;
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

private:
    ScalarType _offsetValue;
    ScalarType _cellWidth;
    std::vector<NodeData> _nodes;
};