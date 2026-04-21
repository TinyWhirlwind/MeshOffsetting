#ifndef _SSBVH_H
#define _SSBVH_H
#include"mymesh.h"
#include "BVH.h"
//最近距离查询
enum BV
{
    PSS = 0,
    LSS,
    RSS,
};
template <class Primitive>
class SSBVH : public BVH<Primitive>
{
public:
    SSBVH(const std::vector<Primitive>& primitives, int maxLeafSize = 4) : BVH<Primitive>(primitives, maxLeafSize) {}
    virtual ~SSBVH() {}
};
#endif