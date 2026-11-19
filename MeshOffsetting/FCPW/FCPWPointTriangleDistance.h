/**
 * @file FCPWPointTriangleDistance.h
 * @author zhangpingping
 * @date 2026-11-18
 * @brief Point-to-triangle closest point and distance calculation.
 *
 * The Voronoi-region tests are adapted from FCPW v1.2.0
 * (include/fcpw/geometry/triangles.inl), which cites
 * Real-Time Collision Detection by Christer Ericson.
 * FCPW is distributed under the MIT License. See LICENSE-FCPW.txt.
 */
#ifndef MESHOFFSETTING_FCPW_POINT_TRIANGLE_DISTANCE_H
#define MESHOFFSETTING_FCPW_POINT_TRIANGLE_DISTANCE_H

#include "base_type.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace meshoffsetting
{
    enum class FCPWTriangleFeature
    {
        Face,
        EdgeAB,
        EdgeBC,
        EdgeCA,
        VertexA,
        VertexB,
        VertexC
    };

    struct FCPWPointTriangleResult
    {
        Point3m closest_point = Point3m(0, 0, 0);

        // The three values are the weights of vertices A, B and C.
        std::array<Scalarm, 3> barycentric = {
            Scalarm(1), Scalarm(0), Scalarm(0)
        };

        Scalarm squared_distance =
            std::numeric_limits<Scalarm>::max();
        FCPWTriangleFeature feature =
            FCPWTriangleFeature::VertexA;

        Scalarm Distance() const
        {
            return std::sqrt(std::max(Scalarm(0), squared_distance));
        }
    };

    namespace fcpw_detail
    {
        inline FCPWPointTriangleResult MakeResult(
            const Point3m& query,
            const Point3m& closest_point,
            const std::array<Scalarm, 3>& barycentric,
            FCPWTriangleFeature feature)
        {
            FCPWPointTriangleResult result;
            result.closest_point = closest_point;
            result.barycentric = barycentric;
            result.feature = feature;

            const Point3m difference = query - closest_point;
            result.squared_distance = difference * difference;
            return result;
        }

        inline void TestSegment(
            const Point3m& query,
            const Point3m& start,
            const Point3m& end,
            const std::array<Scalarm, 3>& start_barycentric,
            const std::array<Scalarm, 3>& end_barycentric,
            FCPWTriangleFeature start_feature,
            FCPWTriangleFeature end_feature,
            FCPWTriangleFeature edge_feature,
            FCPWPointTriangleResult& best)
        {
            const Point3m edge = end - start;
            const Scalarm length_squared = edge * edge;

            Scalarm t = Scalarm(0);
            if (length_squared > Scalarm(0))
            {
                t = ((query - start) * edge) / length_squared;
                t = std::clamp(t, Scalarm(0), Scalarm(1));
            }

            std::array<Scalarm, 3> barycentric;
            for (int i = 0; i < 3; ++i)
            {
                barycentric[i] =
                    start_barycentric[i] * (Scalarm(1) - t) +
                    end_barycentric[i] * t;
            }

            FCPWTriangleFeature feature = edge_feature;
            if (t == Scalarm(0))
            {
                feature = start_feature;
            }
            else if (t == Scalarm(1))
            {
                feature = end_feature;
            }

            const Point3m closest_point = start + edge * t;
            FCPWPointTriangleResult candidate = MakeResult(
                query, closest_point, barycentric, feature);

            if (candidate.squared_distance < best.squared_distance)
            {
                best = candidate;
            }
        }

        inline FCPWPointTriangleResult ClosestPointOnDegenerateTriangle(
            const Point3m& query,
            const Point3m& a,
            const Point3m& b,
            const Point3m& c)
        {
            const std::array<Scalarm, 3> barycentric_a = {
                Scalarm(1), Scalarm(0), Scalarm(0)
            };
            const std::array<Scalarm, 3> barycentric_b = {
                Scalarm(0), Scalarm(1), Scalarm(0)
            };
            const std::array<Scalarm, 3> barycentric_c = {
                Scalarm(0), Scalarm(0), Scalarm(1)
            };

            FCPWPointTriangleResult best;
            TestSegment(
                query, a, b,
                barycentric_a, barycentric_b,
                FCPWTriangleFeature::VertexA,
                FCPWTriangleFeature::VertexB,
                FCPWTriangleFeature::EdgeAB,
                best);
            TestSegment(
                query, a, c,
                barycentric_a, barycentric_c,
                FCPWTriangleFeature::VertexA,
                FCPWTriangleFeature::VertexC,
                FCPWTriangleFeature::EdgeCA,
                best);
            TestSegment(
                query, b, c,
                barycentric_b, barycentric_c,
                FCPWTriangleFeature::VertexB,
                FCPWTriangleFeature::VertexC,
                FCPWTriangleFeature::EdgeBC,
                best);
            return best;
        }
    }

    inline FCPWPointTriangleResult ClosestPointOnTriangle(
        const Point3m& query,
        const Point3m& a,
        const Point3m& b,
        const Point3m& c)
    {
        const Point3m ab = b - a;
        const Point3m ac = c - a;
        const Point3m bc = c - b;

        const Scalarm ab_squared = ab * ab;
        const Scalarm ac_squared = ac * ac;
        const Scalarm bc_squared = bc * bc;
        const Scalarm ab_dot_ac = ab * ac;

        // |AB x AC|^2, computed using dot products.
        const Scalarm area_squared = std::max(
            Scalarm(0),
            ab_squared * ac_squared - ab_dot_ac * ab_dot_ac);
        const Scalarm longest_edge_squared = std::max({
            ab_squared, ac_squared, bc_squared
        });
        const Scalarm degenerate_limit =
            std::numeric_limits<Scalarm>::epsilon() * Scalarm(64) *
            longest_edge_squared * longest_edge_squared;

        if (area_squared <= degenerate_limit)
        {
            return fcpw_detail::ClosestPointOnDegenerateTriangle(
                query, a, b, c);
        }

        const Point3m aq = query - a;
        const Scalarm d1 = ab * aq;
        const Scalarm d2 = ac * aq;
        if (d1 <= Scalarm(0) && d2 <= Scalarm(0))
        {
            return fcpw_detail::MakeResult(
                query, a,
                { Scalarm(1), Scalarm(0), Scalarm(0) },
                FCPWTriangleFeature::VertexA);
        }

        const Point3m bq = query - b;
        const Scalarm d3 = ab * bq;
        const Scalarm d4 = ac * bq;
        if (d3 >= Scalarm(0) && d4 <= d3)
        {
            return fcpw_detail::MakeResult(
                query, b,
                { Scalarm(0), Scalarm(1), Scalarm(0) },
                FCPWTriangleFeature::VertexB);
        }

        const Point3m cq = query - c;
        const Scalarm d5 = ab * cq;
        const Scalarm d6 = ac * cq;
        if (d6 >= Scalarm(0) && d5 <= d6)
        {
            return fcpw_detail::MakeResult(
                query, c,
                { Scalarm(0), Scalarm(0), Scalarm(1) },
                FCPWTriangleFeature::VertexC);
        }

        const Scalarm vc = d1 * d4 - d3 * d2;
        if (vc <= Scalarm(0) &&
            d1 >= Scalarm(0) && d3 <= Scalarm(0))
        {
            const Scalarm v = d1 / (d1 - d3);
            return fcpw_detail::MakeResult(
                query, a + ab * v,
                { Scalarm(1) - v, v, Scalarm(0) },
                FCPWTriangleFeature::EdgeAB);
        }

        const Scalarm vb = d5 * d2 - d1 * d6;
        if (vb <= Scalarm(0) &&
            d2 >= Scalarm(0) && d6 <= Scalarm(0))
        {
            const Scalarm w = d2 / (d2 - d6);
            return fcpw_detail::MakeResult(
                query, a + ac * w,
                { Scalarm(1) - w, Scalarm(0), w },
                FCPWTriangleFeature::EdgeCA);
        }

        const Scalarm va = d3 * d6 - d5 * d4;
        if (va <= Scalarm(0) &&
            (d4 - d3) >= Scalarm(0) &&
            (d5 - d6) >= Scalarm(0))
        {
            const Scalarm w =
                (d4 - d3) /
                ((d4 - d3) + (d5 - d6));
            return fcpw_detail::MakeResult(
                query, b + bc * w,
                { Scalarm(0), Scalarm(1) - w, w },
                FCPWTriangleFeature::EdgeBC);
        }

        const Scalarm inverse_denom =
            Scalarm(1) / (va + vb + vc);
        const Scalarm v = vb * inverse_denom;
        const Scalarm w = vc * inverse_denom;
        return fcpw_detail::MakeResult(
            query,
            a + ab * v + ac * w,
            { Scalarm(1) - v - w, v, w },
            FCPWTriangleFeature::Face);
    }
}

#endif
