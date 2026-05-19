#ifndef _MESH_OFFSETING_H_
#define _MESH_OFFSETING_H_
#include <optional>
#include "mymesh.h"
#include "OffsetGrid.h"
#include "BVH.h"

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
    CFaceO* FindSameTriangle(const OffsetGrid<float>::NodeData& n0, const OffsetGrid<float>::NodeData& n1);

private:
    CMeshO _mesh;
    Params _params;
    std::vector<CFaceO*> _faces;
    std::optional<BVH<CFaceO*>> _bvh;
    std::optional<OffsetGrid<float>> _grid;
};
#endif 
