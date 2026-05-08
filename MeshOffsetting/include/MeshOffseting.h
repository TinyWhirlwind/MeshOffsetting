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
    float QuerySignedDistance(const Point3m& p) const;
    bool IsInvalidByUnsignedDistance(float unsignedDistance, float radius) const;
    bool IsInvalidBySignedDistance(float signedDistance, float radius) const;
    void MarkNodesInWorldBoxInvalid(const Box3m& box);
    void ApplyOctreeFilter();

private:
    CMeshO _mesh;
    Params _params;
    std::vector<CFaceO*> _faces;
    std::optional<BVH<CFaceO*>> _bvh;
    std::optional<OffsetGrid<float>> _grid;
};
#endif 
