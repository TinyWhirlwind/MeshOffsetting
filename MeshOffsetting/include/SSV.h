/**************************************************************************
 *   Copyright (c) 2025
 *
 *   创建日期:2026-04-24
 *   作者:zhangpingping
 *   摘要:swept sphere volume
 *   参考：Flexible Collision Library
 ***************************************************************************/
#include "mymesh.h"
struct RSS
{
    Point3m origin;//矩形的某个点（原点）
    Point3m axis[3];//正交基
    float length[2];//矩形边长
    float radius;
};
struct SSVNode
{
    RSS bv;
    union
    {
        int leafOffset;   // leaf
        int interiorOffset; // interior
    };
    uint16_t nPrimitives;  // 0 => interior, >0 => leaf
    uint8_t axis;          
};

struct BuildSSVNode
{
    RSS bv;
    BuildSSVNode* child[2] = { nullptr,nullptr };

    int leafOffset = -1;
    int nPrimitives = 0;
    int splitAxis = -1;

};