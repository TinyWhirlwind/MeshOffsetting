#include "mymesh.h"
#include <base_type.h>
using namespace vcg;
template <class ScalarType>
struct GridNode
{
    float dist = FLT_MAX;
    int faceId = -1;
    int sign = 0;            // -1 inside, +1 outside, 0 unknown
    bool isV = false;
    float s, t;
};
struct GridCell
{
    int node[8];

};
class Grid
{
public:
    Grid(const Box3m& box, Point3i voxel_size):_box(box), _size(voxel_size)
    {
        GetDimAndVoxel();
    }
    ~Grid();
public:
    void GetDimAndVoxel()
    {
        _dim = _box.max - _box.min;
        for (int i = 0; i < 3; ++i)
        {
            _voxel[i] = _dim[i] / _size[i];
        }
    }

    Point3i GetGridId(const Point3m& p)
    {
        Point3m coord = p - _box.min;
        Point3i result;
        for (int i = 0; i < 3; ++i)
        {
            result[i] = (int)coord[i] / _voxel[i];
        }
        return result;
    }

    /*Matrix44m GetCoord()
    {
        Matrix44 m,t;
        m.SetScale(_voxel);
        t.SetTranslate(_box.min);
    }*/

    void GetGridPBox(const Point3i& id, Box3m box)
    {
        Point3m coord;
        for (int i = 0; i < 3; ++i)
        {
            coord[i] = id[i] * _voxel[i];
        }
        coord += box.min;
        box.min = coord;
        box.max = coord + _voxel;
    }




public:
    Box3m _box;
    Point3m _dim;
    Point3m _voxel;
    Point3i _size;


};