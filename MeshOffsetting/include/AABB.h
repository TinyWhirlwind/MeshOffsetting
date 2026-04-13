#include "mymesh.h"
#include <algorithm>
class AABBBox :public Box3m
{
public:
    unsigned int GetMaxAxis()
    {
        _dim = _diag.X() > _diag.Y() ? 0 : (_diag.Y() > _diag.Z() ? 1 : 2);
        return _dim;
    }

   /* void Merge(const AABBBox& other)
    {
        this->min.X() = std::min(this->min.X(), other.min.X());
        this->min.Y() = std::min(this->min.Y(), other.min.Y());
        this->min.Z() = std::min(this->min.Z(), other.min.Z());

        this->max.X() = std::max(this->max.X(), other.max.X());
        this->max.Y() = std::max(this->max.Y(), other.max.Y());
        this->max.Z() = std::max(this->max.Z(), other.max.Z());
    }*/

    AABBBox Merge(const AABBBox& box, const AABBBox& other)
    {
        AABBBox result;
        result.min.X() = std::min(box.min.X(), other.min.X());
        result.min.Y() = std::min(box.min.Y(), other.min.Y());
        result.min.Z() = std::min(box.min.Z(), other.min.Z());

        result.max.X() = std::max(box.max.X(), other.max.X());
        result.max.Y() = std::max(box.max.Y(), other.max.Y());
        result.max.Z() = std::max(box.max.Z(), other.max.Z());
        return result;
    }

    Point3m NormalizePointToBounds(Point3m p)
    {
        Point3m vec = p - this->min;
        if (this->max.X() > this->min.X())
            vec.X() /= this->max.X() - this->min.X();
        if (this->max.Y() > this->min.Y())
            vec.Y() /= this->max.Y() - this->min.Y();
        if (this->max.Z() > this->min.Z())
            vec.Z() /= this->max.Z() - this->min.Z();
        return vec;
    }

    Scalarm SurfaceArea()
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
