#include "CMeshO.h"
#include "SlicerSDF.h"
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/import_stl.h>
#include <wrap/io_trimesh/import_obj.h>
#include <vcg/complex/algorithms/clean.h>
#include <chrono>
#include <iostream>
#include <omp.h>
using namespace meshoffsetting;
int main()
{
    //omp_set_num_threads(1);
    const char* input_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/8916_9cm_x1.obj";
    const char* output_path ="D:/AppData/MyProject/MeshOffsetting/MeshOffsetting/model/8916_9cm_x1.ply";
    CMeshO mesh;
    mesh.face.EnableWedgeNormal();
    int load_mask = 0;
    //vcg::tri::io::ImporterSTL<CMeshO>::Open(mesh, input_path, load_mask);
    vcg::tri::io::ImporterOBJ<CMeshO>::Open(mesh, input_path, load_mask);
    CMeshO point_cloud;
#if 0
    NarrowBandSDF sdf(mesh, 0.,129,8);
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
        else if(nodes[i].state != GridNodeState::Unknown)
        {
            //auto vt = vcg::tri::Allocator<CMeshO>::AddVertex(point_cloud, pos);
            //vt->C() = vcg::Color4b(0, 255, 0, 255);
        }
    }

#endif
    /*vcg::tri::Clean<CMeshO>::RemoveDuplicateVertex(mesh);
    vcg::tri::Clean<CMeshO>::RemoveUnreferencedVertex(mesh);
    vcg::tri::Allocator<CMeshO>::CompactEveryVector(mesh);*/
    mesh.updateBoxAndNormals();
    /*mesh.face.EnableMark();
    vcg::tri::InitFaceIMark(mesh);
    std::vector<CFaceO*> intersectFaces;
    bool hasIntersection = vcg::tri::Clean<CMeshO>::SelfIntersections(mesh, intersectFaces);
    if (!hasIntersection)
    {
        std::cout << "没有自相交" << std::endl;
    }*/
    Scalarm thickness = 0.02;
    int n_layer = (mesh.bbox.max.Z() - mesh.bbox.min.Z())/ thickness;
    auto start_time = std::chrono::steady_clock::now();
    SlicerSDF sdf(mesh, 0, 513);
    sdf.StartSlicerAll(n_layer, thickness);
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "SDF and slicing time: "<< elapsed_time.count() << " seconds\n";
    const auto& result = sdf.GetLayers();
    for (int i =0;i< result.size();++i)
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
    vcg::tri::io::ExporterPLY<CMeshO>::Save(point_cloud, output_path, vcg::tri::io::Mask::IOM_VERTCOLOR);
    return 0;
}
