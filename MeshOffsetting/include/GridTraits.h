#pragma once

#include "mymesh.h"

struct FaceBBoxGetter
{
    Box3m operator()(const CFaceO& f) const
    {
        Box3m box;
        box.Set(f.cP(0));
        box.Add(f.cP(1));
        box.Add(f.cP(2));
        return box;
    }
};

struct PointBBoxGetter
{
    Box3m operator()(const Point3m& p) const
    {
        Box3m box;
        box.Set(p);
        return box;
    }
};

struct VertexBBoxGetter
{
    Box3m operator()(const CVertexO& v) const
    {
        Box3m box;
        box.Set(v.cP());
        return box;
    }
};

struct EdgeBBoxGetter
{
    Box3m operator()(const CEdgeO& e) const
    {
        Box3m box;
        box.Set(e.V(0)->cP());
        box.Add(e.V(1)->cP());
        return box;
    }
};

