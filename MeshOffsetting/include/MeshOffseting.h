#ifndef _MESH_OFFSETING_H_
#define _MESH_OFFSETING_H_
#include "mymesh.h"
#include "OffsetGrid.h"
#define Scalarm double
class MeshOffseting
{
public:
    MeshOffseting(CMeshO mesh);
    ~MeshOffseting();

public:
    void setVoxelSize();
    void SetVolumeDim();
    void setOffsetingDistance(Scalarm distance);
private:
    CMeshO offsetMesh;
};
#endif 