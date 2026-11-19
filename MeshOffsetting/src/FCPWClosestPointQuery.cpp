#include "FCPWClosestPointQuery.h"

#include <fcpw/fcpw.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace meshoffsetting
{
    struct FCPWClosestPointQuery::Impl
    {
        fcpw::Scene<3> scene;
        std::vector<int> primitive_to_face;
    };

    namespace
    {
        bool IsFinitePoint(const Point3m& point)
        {
            return std::isfinite(point.X()) &&
                std::isfinite(point.Y()) &&
                std::isfinite(point.Z());
        }

        ClosestFeature ClassifyFeature(const fcpw::Vector2& uv)
        {
            // FCPW stores the A and B barycentric weights in uv.
            const float weight_a = uv[0];
            const float weight_b = uv[1];
            const float weight_c = 1.0f - weight_a - weight_b;
            const float tolerance =
                std::numeric_limits<float>::epsilon() * 32.0f;

            const bool a_is_zero = std::abs(weight_a) <= tolerance;
            const bool b_is_zero = std::abs(weight_b) <= tolerance;
            const bool c_is_zero = std::abs(weight_c) <= tolerance;

            if (b_is_zero && c_is_zero)
                return ClosestFeature::VertexA;
            if (a_is_zero && c_is_zero)
                return ClosestFeature::VertexB;
            if (a_is_zero && b_is_zero)
                return ClosestFeature::VertexC;
            if (c_is_zero)
                return ClosestFeature::EdgeAB;
            if (a_is_zero)
                return ClosestFeature::EdgeBC;
            if (b_is_zero)
                return ClosestFeature::EdgeCA;

            return ClosestFeature::Face;
        }
    }

    FCPWClosestPointQuery::FCPWClosestPointQuery() = default;

    FCPWClosestPointQuery::FCPWClosestPointQuery(const CMeshO& mesh)
    {
        Initialize(mesh);
    }

    FCPWClosestPointQuery::~FCPWClosestPointQuery() = default;

    void FCPWClosestPointQuery::Initialize(const CMeshO& mesh)
    {
        _impl.reset();

        std::vector<fcpw::Vector3> positions;
        std::vector<int> vertex_map(mesh.vert.size(), -1);
        positions.reserve(static_cast<std::size_t>(mesh.vn));

        for (std::size_t vertex_index = 0;
             vertex_index < mesh.vert.size();
             ++vertex_index)
        {
            const CVertexO& vertex = mesh.vert[vertex_index];
            if (vertex.IsD() || !IsFinitePoint(vertex.P()))
                continue;

            vertex_map[vertex_index] =
                static_cast<int>(positions.size());
            positions.emplace_back(
                static_cast<float>(vertex.P().X()),
                static_cast<float>(vertex.P().Y()),
                static_cast<float>(vertex.P().Z()));
        }

        std::vector<fcpw::Vector3i> triangles;
        std::vector<int> primitive_to_face;
        triangles.reserve(static_cast<std::size_t>(mesh.fn));
        primitive_to_face.reserve(static_cast<std::size_t>(mesh.fn));

        for (std::size_t face_index = 0;
             face_index < mesh.face.size();
             ++face_index)
        {
            const CFaceO& face = mesh.face[face_index];
            if (face.IsD())
                continue;

            int indices[3] = { -1, -1, -1 };
            bool valid_face = true;
            for (int corner = 0; corner < 3; ++corner)
            {
                const CVertexO* face_vertex = face.cV(corner);
                if (face_vertex == nullptr)
                {
                    valid_face = false;
                    break;
                }

                const std::size_t original_vertex_index =
                    vcg::tri::Index(mesh, face_vertex);
                if (original_vertex_index >= vertex_map.size() ||
                    vertex_map[original_vertex_index] < 0)
                {
                    valid_face = false;
                    break;
                }

                indices[corner] = vertex_map[original_vertex_index];
            }

            if (!valid_face)
                continue;

            const Point3m edge_ab = face.P(1) - face.P(0);
            const Point3m edge_ac = face.P(2) - face.P(0);
            const Point3m edge_bc = face.P(2) - face.P(1);
            const Scalarm longest_edge_squared = std::max({
                edge_ab.SquaredNorm(),
                edge_ac.SquaredNorm(),
                edge_bc.SquaredNorm()
            });
            const Scalarm area_squared =
                (edge_ab ^ edge_ac).SquaredNorm();
            const Scalarm degenerate_limit =
                std::numeric_limits<Scalarm>::epsilon() *
                longest_edge_squared * longest_edge_squared;
            if (!std::isfinite(area_squared) ||
                area_squared <= degenerate_limit)
            {
                continue;
            }

            triangles.emplace_back(
                indices[0], indices[1], indices[2]);
            primitive_to_face.push_back(
                static_cast<int>(face_index));
        }

        if (triangles.empty())
            return;

        std::unique_ptr<Impl> new_impl = std::make_unique<Impl>();
        new_impl->primitive_to_face = std::move(primitive_to_face);
        new_impl->scene.setObjectCount(1);
        new_impl->scene.setObjectVertices(positions, 0);
        new_impl->scene.setObjectTriangles(triangles, 0);
        new_impl->scene.build(
            fcpw::AggregateType::Bvh_SurfaceArea,
            false,
            false,
            false);

        _impl = std::move(new_impl);
    }

    ClosestPointResult FCPWClosestPointQuery::QueryClosestPoint(
        const Point3m& q) const
    {
        ClosestPointResult result;
        if (Empty())
            return result;

        const fcpw::Vector3 query_point(
            static_cast<float>(q.X()),
            static_cast<float>(q.Y()),
            static_cast<float>(q.Z()));
        fcpw::Interaction<3> interaction;

        const bool found =
            _impl->scene.findClosestPoint(query_point, interaction);
        if (!found ||
            interaction.primitiveIndex < 0 ||
            interaction.primitiveIndex >=
                static_cast<int>(_impl->primitive_to_face.size()) ||
            !std::isfinite(interaction.d) ||
            !interaction.p.allFinite() ||
            !interaction.uv.allFinite())
        {
            return result;
        }

        result.faceId =
            _impl->primitive_to_face[interaction.primitiveIndex];
        result.closest_point = Point3m(
            static_cast<Scalarm>(interaction.p[0]),
            static_cast<Scalarm>(interaction.p[1]),
            static_cast<Scalarm>(interaction.p[2]));
        result.unsigned_distance =
            static_cast<Scalarm>(interaction.d);
        result.unsigned_distance_squared =
            result.unsigned_distance * result.unsigned_distance;
        result.feature = ClassifyFeature(interaction.uv);
        return result;
    }

    bool FCPWClosestPointQuery::Empty() const
    {
        return _impl == nullptr ||
            _impl->primitive_to_face.empty();
    }
}
