#include "mymesh.h"
#include <base_type.h>
template <class ScalarType>
class Grid
{
public:
    Grid(const Box3m& box, ScalarType voxel_size):_box(box)
    {}
    ~Grid();
public:

    void GetDimAndVoxel()
    {
    }

    Point3m GetGridId(const Point3m& p)
    {

    }
public:
    Box3m _box;
    Point3m _dim[3];
    Point3m _voxel[3];


};