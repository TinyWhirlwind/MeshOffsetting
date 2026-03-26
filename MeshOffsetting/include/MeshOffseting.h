#ifndef _MESH_OFFSETING_H_
#
#include "mymesh.h"
#include "Grid.h"
#define Scalarm double
class MeshOffseting
{
public:
    MeshOffseting(const CMeshO& mesh);
    ~MeshOffseting();

public:
    void setVoxelSize();
    void SetVolumeDim();
    void setOffsetingDistance(Scalarm distance);
private:
    CMeshO offsetMesh;
};