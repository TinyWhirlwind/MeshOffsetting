/**
 * @file FCPWClosestPointQuery.h
 * @author zhangpingping
 * @date 2026-11-18
 * @brief FCPW accelerated closest-point query for triangle meshes.
 */
#ifndef MESHOFFSETTING_FCPW_CLOSEST_POINT_QUERY_H
#define MESHOFFSETTING_FCPW_CLOSEST_POINT_QUERY_H

#include "SSVH.h"

#include <memory>

namespace meshoffsetting
{
    class FCPWClosestPointQuery
    {
    public:
        explicit FCPWClosestPointQuery(const CMeshO& mesh);
        FCPWClosestPointQuery();
        ~FCPWClosestPointQuery();

        void Initialize(const CMeshO& mesh);

        ClosestPointResult QueryClosestPoint(const Point3m& q) const;
        bool Empty() const;

    private:
        struct Impl;
        std::unique_ptr<Impl> _impl;
    };
}

#endif
