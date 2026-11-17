#include "CMeshO.h"
#include "SlicerSDF.h"
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/import_stl.h>
#include <wrap/io_trimesh/import_obj.h>
using namespace meshoffsetting;
int main()
{
    const char* input_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/8916_9cm_x1.obj";
    const char* output_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/8916_9cm_x1.ply";
    CMeshO mesh;
    mesh.face.EnableWedgeNormal();
    int load_mask = 0;
    //vcg::tri::io::ImporterSTL<CMeshO>::Open(mesh, input_path, load_mask);
    vcg::tri::io::ImporterOBJ<CMeshO>::Open(mesh, input_path, load_mask);
    CMeshO point_cloud;

    NarrowBandSDF sdf(mesh, 0,129,8);
    auto grid = sdf.GetGird();
    auto nodes = grid->GetNodes();
    for (int i = 0; i < nodes.size(); ++i)
    {
        Point3m pos = grid->NodePosition(i);
        if (nodes[i].state == GridNodeState::NodeDistanceReady)
        {
            auto vt = vcg::tri::Allocator<CMeshO>::AddVertex(point_cloud, pos);
            vt->C() = vcg::Color4b(255, 0, 0, 255);
        }
    } 

#if 0
    mesh.updateBoxAndNormals();
    Scalarm thickness = 0.02;
    int n_layer = (mesh.bbox.max.Z() - mesh.bbox.min.Z())/ thickness;
    SlicerSDF sdf(mesh, -0.01,65);
    sdf.StartSlicerAll(n_layer, thickness);
    std::vector<SliceResult> result = sdf.GetLayers();
    for (int i =0;i< result.size();++i)
    {
        if (i % 2 == 0)
        {
            for (auto& ct : result[i].contours)
            {
                for (auto& pt : ct)
                {
                    auto vt = vcg::tri::Allocator<CMeshO>::AddVertex(point_cloud, pt);
                    vt->C() = vcg::Color4b(255, 0, 0, 255);
                }
            }
        }
        else
        {
            for (auto& ct : result[i].contours)
            {
                for (auto& pt : ct)
                {
                    auto vt = vcg::tri::Allocator<CMeshO>::AddVertex(point_cloud, pt);
                    vt->C() = vcg::Color4b(0, 255, 0, 255);
                }
            }
        }
    }
#endif
    vcg::tri::io::ExporterPLY<CMeshO>::Save(point_cloud, output_path, vcg::tri::io::Mask::IOM_VERTCOLOR);
    return 0;
}
