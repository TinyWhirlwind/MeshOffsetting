#include "mymesh.h"
#include <vcg/space/point3.h>
template <class ScalarType>
class Grid
{
public:
    typdef vcg::Point3<ScalarType> Point3T;
    typdef vcg::Point3i Point3i;
    Box3m _box;
    Point3T _dim[3];
    Point3T _voxel[3];
    Point3i _size[3];

    void GetDimAndVoxel()
    {
        _dim = _box.max - _box.min;
        for (int i = 0; i < 3; ++i)
        {
            _voxel[0] = _dim[0] / _size[0];
        }
    }

    Point3i GetGridId(const Point3T& p)
    {
    }

};