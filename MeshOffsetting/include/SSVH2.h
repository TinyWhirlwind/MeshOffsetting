#ifndef MESHOFFSETTING_SSVH_H
#define MESHOFFSETTING_SSVH_H
#include"CMeshO.h"
#include "BVH.h"
//最近距离查询
enum BV
{
    PSS = 0,
    LSS,
    RSS,
};
template <class Primitive>
class SSVH : public BVH<Primitive>
{
public:
    SSVH(const std::vector<Primitive>& primitives, int maxLeafSize = 4) : BVH<Primitive>(primitives, maxLeafSize) {}
    virtual ~SSBVH() {}
};
#endif    