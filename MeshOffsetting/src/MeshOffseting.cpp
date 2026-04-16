#include "MeshOffseting.h"
MeshOffseting::MeshOffseting(CMeshO mesh)
{
    vcg::tri::UpdateBounding<CMeshO>::Box(mesh);

}

MeshOffseting::~MeshOffseting()
{
}