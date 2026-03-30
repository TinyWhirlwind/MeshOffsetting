#include "mymesh.h"
class AABBBox :public Box3m
{
public:
    unsigned int GetMaxAxis()
    {
        _dim = _diag.X() > _diag.Y() ? 0 : (_diag.Y() > _diag.Z() ? 1 : 2);
        return _dim;
    }

    Scalarm GetSAH()
    {
        return (_diag.X() * _diag.Y() + _diag.X() * _diag.Z() + _diag.Z() * _diag.Y()) * 2;
    }

    Scalarm CalcCost()
    {

    }

    //相交测试
    bool checkOverlapAABB(const AABBBox& other)
    {
        if (other.IsEmpty()|| other.IsNull())
            return false;
        if (other.max.X() < this->min.X() || other.min.X() > this->max.X())
        {
            return false;
        }
        if (other.max.Y() < this->min.Y() || other.min.Y() > this->max.Y())
        {
            return false;
        }
        if (other.max.Z() < this->min.Z() || other.min.Z() > this->max.Z())
        {
            return false;
        }
        return true;
    }

    bool checkOverlapShpere(const AABBBox& other)
    {
        return true;
    }

private:
    Point3m _diag;
    Point3m _center;
    unsigned int _dim;
};
