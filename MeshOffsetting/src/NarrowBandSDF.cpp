#include "CMeshO.h"
#include "NarrowBandSDF.h"
namespace meshoffsetting
{
    NarrowBandSDF::NarrowBandSDF(const CMeshO& mesh, const Scalarm offsetDistance, int nodeCount):
        _ssvh(mesh),
		_grid(mesh.bbox, offsetDistance, nodeCount)
    {
    }
    
    //Bounding Box Filter
    void NarrowBandSDF::ApplyBoundingBoxFilter()
    {
        
    }

    //Signed Distance Filter
    void NarrowBandSDF::ApplySignedDistanceFilter()
    {

    }

    //Octree Filter
    void NarrowBandSDF::ApplyOctreeFilter()
    {

    }
}