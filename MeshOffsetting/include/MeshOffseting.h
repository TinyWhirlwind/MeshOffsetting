#ifndef _MESH_OFFSETING_H_
#define _MESH_OFFSETING_H_
#include <optional>
#include "mymesh.h"
#include "OffsetGrid.h"
#include "BVH.h"

enum class TriangleRegion
{
    Face = 0,

    Edge01,
    Edge12,
    Edge20,

    Vertex0,
    Vertex1,
    Vertex2,

    Unknown
};

struct CanonicalFeature
{
    int _faceId = -1;
    ClosestType _type = ClosestType::Unknown;
    int _feature = -1; // 这个 feature 是相对 sameFace 的局部编号
};

class MeshOffseting
{
public:
    struct Params
    {
        float voxelSize = 1.0;
        float offsetDistance = 0.0;
        //int bvhLeafSize = 10;
    };

    explicit MeshOffseting(const CMeshO& mesh, const Params& params);
    ~MeshOffseting();

public:
    void SetMesh(const CMeshO& mesh);
    void SetVoxelSize(float cellWidth);
    void SetOffsetDistance(float offset);
    void SetParams(const Params& params);
    
    void RebuildBVH();
    void RebuildGrid();
    void Run();

private:
    bool QuerySignedDistance(const Point3m& p, QueryResult& qr) const;
    bool IsInvalidByUnsignedDistance(float unsignedDistance, float radius) const;
    bool IsInvalidBySignedDistance(float signedDistance, float radius) const;
    void MarkNodesInWorldBoxInvalid(const Box3m& box);

    void ApplySignedDistanceFilter();
    void ApplyOctreeFilter();
    void FindIntersectionPointOnGridEdge();
    //找到相同三角形，并将相交边的两个最近点映射到同一个三角形上
    bool FindSameTriangle(const OffsetGrid<float>::NodeData& n0, const OffsetGrid<float>::NodeData& n1, CanonicalFeature& cf0, CanonicalFeature& cf1);
    TriangleRegion ToTriangleRegion(const QueryResult& qr);
    void SolveByBisection(OffsetGrid<float>::GridEdge& edge);
    void SolveByAnalyticalSolution(OffsetGrid<float>::GridEdge& edge, CanonicalFeature cf0, CanonicalFeature cf1);

    bool judgeEdgeInTriangleOrder(int v0, int v1, CFaceO* f,int& order);
    bool judgeVertexInTriangleOrder(int v, CFaceO* f,int& order);
private:
    CMeshO _mesh;
    Params _params;
    std::vector<CFaceO*> _faces;
    std::optional<BVH<CFaceO*>> _bvh;
    std::optional<OffsetGrid<float>> _grid;
};
#endif 
