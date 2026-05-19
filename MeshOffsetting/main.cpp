#include "mymesh.h"
#include "BVH.h"
using namespace vcg;
int main()
{
    CMeshO mesh;
    tri::Allocator<CMeshO>::AddVertices(mesh, 3);
    tri::Allocator<CMeshO>::AddFaces(mesh, 1);
    mesh.vert[0].P() = Point3f(0, 0, 0);
    mesh.vert[1].P() = Point3f(1, 0, 0);
    mesh.vert[2].P() = Point3f(0, 1, 0);
    
    mesh.face[0].V(0) = &mesh.vert[0];
    mesh.face[0].V(1) = &mesh.vert[1];
    mesh.face[0].V(2) = &mesh.vert[2];
    tri::UpdateBounding<CMeshO>::Box(mesh);
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
    tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(mesh);

	Point3f p0 = Point3f(0.5, 0.5, 1);
    std::vector<CFaceO*> faces;
    for (int i = 0; i < mesh.fn; ++i) 
    {
        faces.push_back(&mesh.face[i]);
    }
    BVH<CFaceO*> bvh(faces);
    QueryResult result;
	result._dist = std::numeric_limits<float>::max();
	bvh.QueryClosestPoint(p0, result);
	std::cout << "Closest point: " << result._closestPoint[0] << ", " << result._closestPoint[1] << ", " << result._closestPoint[2] << std::endl;
	std::cout << "Closest dist: " << result._dist << std::endl;
    std::cout << "Hello World!" << std::endl;
}
