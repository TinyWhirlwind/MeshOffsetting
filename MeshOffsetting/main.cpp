#include "CMeshO.h"
#include "NarrowBandSDF.h"
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/import_stl.h>
using namespace meshoffsetting;
int main()
{
    const char* input_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/bunny_10k.stl";
    const char* output_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/bunny_10k.ply";
    CMeshO mesh;
    int load_mask = 0;
    vcg::tri::io::ImporterSTL<CMeshO>::Open(mesh, input_path, load_mask);
    NarrowBandSDF sdf(mesh, 0.01,129,8);
    auto grid = sdf.GetGird();
    auto nodes = grid->GetNodes();
    CMeshO point_cloud;
    for (int i = 0; i < nodes.size(); ++i)
    {
        Point3m pos = grid->NodePosition(i);
        if (nodes[i].state == GridNodeState::NodeDistanceReady)
        {
            auto vt = vcg::tri::Allocator<CMeshO>::AddVertex(point_cloud, pos);
            vt->C() = vcg::Color4b(255, 0, 0, 255);
        }
    }
    vcg::tri::io::ExporterPLY<CMeshO>::Save(point_cloud, output_path, vcg::tri::io::Mask::IOM_VERTCOLOR);
    return 0;
}
