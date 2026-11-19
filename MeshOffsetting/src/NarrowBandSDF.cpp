#include "CMeshO.h"
#include "Tolerance.h"
#include "NarrowBandSDF.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace meshoffsetting
{
    NarrowBandSDF::NarrowBandSDF(
        const CMeshO& mesh,
        const Scalarm offsetDistance,
        int nodeCount,
        int blockWidth) :
        _offset_distance(offsetDistance),
        _mesh(mesh),
        _block_width(blockWidth)
    {
        if (nodeCount <= 3)
        {
            throw std::invalid_argument("NarrowBandSDF node count must be greater than 3");
        }
        if (blockWidth <= 0 ||
            (blockWidth & (blockWidth - 1)) != 0)
        {
            throw std::invalid_argument("NarrowBandSDF block width must be a positive power of two");
        }
        if (!std::isfinite(_offset_distance))
        {
            throw std::invalid_argument("NarrowBandSDF offset distance must be finite");
        }

        vcg::tri::UpdateBounding<CMeshO>::Box(_mesh);
        if (_mesh.fn > 0)
        {
            _mesh.face.EnableFFAdjacency();
            vcg::tri::UpdateTopology<CMeshO>::FaceFace(_mesh);
            vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
            vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(_mesh);
        }

        _closest_point_query.Initialize(_mesh);
        _grid = UniformGrid(_mesh.bbox, _offset_distance, _block_width, nodeCount);
        _cell_width = _grid.GetCellWidth();
        Build();
    }

    Point3i NarrowBandSDF::GetGridSize() const
    {
        return _grid.GetSize();
    }

    Box3m NarrowBandSDF::GetGridBox() const
    {
        return _grid.GetBox();
    }

    Scalarm NarrowBandSDF::GetCellWidth() const
    {
        return _grid.GetCellWidth();
    }

    bool NarrowBandSDF::TryGetNodeDistance(
        const Point3i& node_id,
        Scalarm& signed_distance) const
    {
        return _grid.TryGetNodeDistance(node_id, signed_distance);
    }

    bool NarrowBandSDF::EvaluateNodeDistance(
        const Point3i& node_id,
        Scalarm& signed_distance) const
    {
        if (TryGetNodeDistance(node_id, signed_distance))
        {
            return true;
        }

        if (!_grid.IsInside(node_id))
        {
            return false;
        }

        return QuerySignedDistance(
            _grid.NodePosition(node_id), signed_distance);
    }

    bool NarrowBandSDF::QuerySignedDistance(
        const Point3m& query_point,
        Scalarm& signed_distance) const
    {
        ClosestPointResult closest_point = _closest_point_query.QueryClosestPoint(query_point);
        bool is_valid = true;
        int sign = CheckSign(closest_point, query_point, is_valid);
        if (!is_valid || !std::isfinite(closest_point.unsigned_distance))
        {
            return false;
        }
        signed_distance = sign * closest_point.unsigned_distance;
        return true;
    }

    void NarrowBandSDF::AngleWeightedPseudoNormal(CVertexO* v, Point3m& normal) const
    {
        vcg::face::VFIterator<CMeshO::FaceType> vfi(v);
        for (; !vfi.End(); ++vfi)
        {
            auto* nf = vfi.F();
            int vid = vfi.I();
            auto v0 = nf->P(vid);
            auto v1 = nf->P((vid + 1) % 3);
            auto v2 = nf->P((vid + 2) % 3);
            auto e1 = v1 - v0;
            auto e2 = v2 - v0;
            if (e1.SquaredNorm() <= kScalarToleranceSquared ||
                e2.SquaredNorm() <= kScalarToleranceSquared)
            {
                continue; // 退化面
            }
            e1.Normalize();
            e2.Normalize();

            float cosTheta = std::clamp(e1.dot(e2), -1.0f, 1.0f);
            float theta = std::acos(cosTheta);
            normal += nf->N() * theta;
        }
    }

    void NarrowBandSDF::AngleWeightedPseudoNormal(const CFaceO* f, int edge_id, Point3m& normal) const
    {
        auto nf = f->FFp(edge_id);
        if (nf != f)
        {
            normal = f->N() + nf->N();
        }
        else
        {
            //边界边
            normal = f->N();
        }
    }

    int NarrowBandSDF::CheckSign(const ClosestPointResult& query_result, const Point3m& center, bool& is_valid) const
    {
        is_valid = true;
        if (query_result.faceId < 0)
        {
            is_valid = false;
            return 0;
        }
        //Angle-weighted pseudonormal
        //根据查询点到网格的最近点，根据最近点的位置（点、边、面），采用三种法向，判断查询点在模型内部还是外部
        int sign = 0;
        const CFaceO* f = &_mesh.face[query_result.faceId];
        Point3m nf(0., 0., 0.);
        ClosestFeature type = query_result.feature;
        switch (type)
        {
        case ClosestFeature::VertexA:
            nf = f->V(0)->N();
            break;
        case ClosestFeature::VertexB:
            nf = f->V(1)->N();
            break;
        case ClosestFeature::VertexC:
            nf = f->V(2)->N();
            break;
        case ClosestFeature::EdgeAB:
            AngleWeightedPseudoNormal(f, 0, nf);
            break;
        case ClosestFeature::EdgeBC:
            AngleWeightedPseudoNormal(f, 1, nf);
            break;
        case ClosestFeature::EdgeCA:
            AngleWeightedPseudoNormal(f, 2, nf);
            break;
        case ClosestFeature::Face:
            nf = f->N();
            break;
        case ClosestFeature::Unknown:
        default:
            is_valid = false;
            return 0;
        }
        if (nf.SquaredNorm() <= kScalarToleranceSquared)
        {
            nf = f->N();
        }
        sign = (center - query_result.closest_point) * nf > 0 ? 1 : -1;
        return sign;
    }

    void NarrowBandSDF::SubdivideBlock(const Point3i& min_node, const Point3i& max_node, int depth)
    {
        if (depth == 0)
            return;
        Point3i mid_node = (min_node + max_node) / 2;
        Scalarm box_circumsphere_radius = (_grid.NodePosition(min_node) - _grid.NodePosition(mid_node)).Norm() / 2;
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                for (int k = 0; k < 2; ++k)
                {
                    Point3i child_min(
                        i == 0 ? min_node.X() : mid_node.X(),
                        j == 0 ? min_node.Y() : mid_node.Y(),
                        k == 0 ? min_node.Z() : mid_node.Z());
                    Point3i child_max(
                        i == 0 ? mid_node.X() : max_node.X(),
                        j == 0 ? mid_node.Y() : max_node.Y(),
                        k == 0 ? mid_node.Z() : max_node.Z());
                    Point3m box_center = (_grid.NodePosition(child_min) + _grid.NodePosition(child_max)) / 2;

                    //判断sub_box的符号距离
                    ClosestPointResult closest_point = _closest_point_query.QueryClosestPoint(box_center);
                    bool query_ok = closest_point.faceId >= 0 && std::isfinite(closest_point.unsigned_distance);
                    if (!query_ok)
                    {
                        continue;
                    }

                    bool sign_valid = true;
                    int sign = CheckSign(closest_point, box_center, sign_valid);
                    if (!sign_valid)
                    {
                        //符号不可信不能删除
                        //MarkNodeRangeInvalid(child_min, child_max);
                        continue;
                    }

                    Scalarm signed_distance = sign * closest_point.unsigned_distance;
                    if (std::abs(signed_distance - _offset_distance) > _cell_width + box_circumsphere_radius)
                    {
                        //设置 Invalid grid nodes
                        MarkNodeRangeInvalid(child_min, child_max);
                        continue;
                    }
                    SubdivideBlock(child_min, child_max, depth - 1);
                }
            }
        }
    }

    bool NarrowBandSDF::IsBoxValidBySignedDistance(Scalarm signed_distance, Scalarm circumsphere_radius) const
    {
        Scalarm d = std::abs(_offset_distance);
        if (_offset_distance > 0)
        {
            Scalarm upper = std::min(
                -d + _cell_width + circumsphere_radius,
                d - _cell_width - circumsphere_radius);
            Scalarm lower = -d - _cell_width - circumsphere_radius;
            if (signed_distance >= lower && signed_distance <= upper)
            {
                return false;
            }
        }
        else
        {
            Scalarm upper = d + _cell_width + circumsphere_radius;
            Scalarm lower = std::max(
                d - _cell_width - circumsphere_radius,
                -d + _cell_width + circumsphere_radius);
            if (signed_distance >= lower && signed_distance <= upper)
            {
                return false;
            }
        }
        return true;
    }

    void NarrowBandSDF::MarkNodeRangeInvalid(const Point3i& min_node, const Point3i& max_node)
    {
        _grid.DeactivateNodeRange(min_node, max_node);
    }

    void NarrowBandSDF::InitializeActiveNodes()
    {
        _grid.ResetNodes();

        for (const GridBlock& block : _blocks)
        {
            if (block.valid)
            {
                _grid.ActivateNodeRange(block.min_node, block.max_node);
            }
        }

        // Invalid wins for nodes shared by valid and invalid blocks.
        for (const GridBlock& block : _blocks)
        {
            if (!block.valid)
            {
                _grid.DeactivateNodeRange(block.min_node, block.max_node);
            }
        }
    }

    void NarrowBandSDF::ComputeRemainingNodeDistances()
    {
        const std::ptrdiff_t brick_count = static_cast<std::ptrdiff_t>(_grid.StoredBrickCount());

#pragma omp parallel for schedule(guided, 1)
        for (std::ptrdiff_t brick_index = 0; brick_index < brick_count; ++brick_index)
        {
            for (std::size_t local_index = 0; local_index < UniformGrid::kBrickNodeCount; ++local_index)
            {
                Point3i node_id = _grid.StoredNodeCoord(static_cast<std::size_t>(brick_index), local_index);
                if (!_grid.IsInside(node_id))
                {
                    continue;
                }

                GridNode& cur_node = _grid.StoredNode(static_cast<std::size_t>(brick_index), local_index);
                if (cur_node.state != GridNodeState::Unknown)
                {
                    continue;
                }

                Point3m node_pos = _grid.NodePosition(node_id);
                ClosestPointResult closest_point = _closest_point_query.QueryClosestPoint(node_pos);
                bool is_valid = true;
                int sign = CheckSign(closest_point, node_pos, is_valid);
                if (!is_valid || !std::isfinite(closest_point.unsigned_distance))
                {
                    cur_node.state = GridNodeState::NodeInvalid;
                    continue;
                }

                cur_node.state = GridNodeState::NodeDistanceReady;
                cur_node.signed_distance = sign * closest_point.unsigned_distance;
            }
        }
    }


    void NarrowBandSDF::Build()
    {
        CreateBlocks();
        ApplyBoundingBoxFilter();
        ApplySignedDistanceFilter();
        InitializeActiveNodes();
        ApplyOctreeFilter();
        ComputeRemainingNodeDistances();

        _blocks.clear();
        _blocks.shrink_to_fit();
    }

    void NarrowBandSDF::CreateBlocks()
    {
        Point3i size = _grid.GetSize();

        int count_x = (size.X() - 1) / static_cast<int>(_block_width);
        int count_y = (size.Y() - 1) / static_cast<int>(_block_width);
        int count_z = (size.Z() - 1) / static_cast<int>(_block_width);
        _blocks.clear();
        _blocks.reserve(
            static_cast<std::size_t>(count_x) *
            static_cast<std::size_t>(count_y) *
            static_cast<std::size_t>(count_z));
        for (int i = 0; i < count_x; ++i)
        {
            for (int j = 0; j < count_y; ++j)
            {
                for (int k = 0; k < count_z; ++k)
                {
                    GridBlock block;
                    block.min_node = Point3i(i * _block_width,
                                            j * _block_width,
                                            k * _block_width);
                    block.max_node = Point3i((i + 1) * _block_width,
                                            (j + 1) * _block_width,
                                            (k + 1) * _block_width);
					Point3m min_point = _grid.NodePosition(block.min_node);
					Point3m max_point = _grid.NodePosition(block.max_node);
					//block.radius = (max_point - min_point).Norm() * 0.5;
					block.center = (min_point + max_point) * 0.5;
                    _blocks.push_back(block);
                }
            }
        }
    }

    //Bounding Box Filter
    void NarrowBandSDF::ApplyBoundingBoxFilter()
    {
        _block_circumsphere_radius = std::sqrt(3.0) * _cell_width * _block_width * 0.5;
        const std::ptrdiff_t block_count = static_cast<std::ptrdiff_t>(_blocks.size());

#pragma omp parallel for schedule(guided, 8)
        for (std::ptrdiff_t block_index = 0; block_index < block_count; ++block_index)
        {
            GridBlock& cubic = _blocks[static_cast<std::size_t>(block_index)];
            ClosestPointResult closest_point = _closest_point_query.QueryClosestPoint(cubic.center);
            //计算block中心点到网格的无符号距离场剔除不在窄带内的block
            bool query_ok = closest_point.faceId >= 0 && std::isfinite(closest_point.unsigned_distance);
            if (!query_ok)
            {
                cubic.valid = true;
                continue;
            }
            if (std::abs(closest_point.unsigned_distance - std::abs(_offset_distance)) > _block_circumsphere_radius + _cell_width)
            {
                cubic.valid = false;
            }
            else
            {
				cubic.valid = true;
				cubic.closestPointResult = closest_point;
            }
        }

    }

    //Signed Distance Filter
    void NarrowBandSDF::ApplySignedDistanceFilter()
    {
        for (auto& cubic : _blocks)
        {
            if (cubic.closestPointResult.faceId < 0 ||
                !std::isfinite(cubic.closestPointResult.unsigned_distance))
            {
                continue;
            }

            if (!cubic.valid)continue;
            bool sign_valid = true;
            int sign = CheckSign(cubic.closestPointResult, cubic.center, sign_valid);
            if (!sign_valid)
            {
                // An uncertain sign cannot safely reject this block.
                continue;
            }
			Scalarm signed_distance = sign * cubic.closestPointResult.unsigned_distance;
            if (!IsBoxValidBySignedDistance(signed_distance, _block_circumsphere_radius))
            {
                cubic.valid = false;
            }
        }
    }

    /*
    Octree Filter : To further exclude the invalid grid nodes
    from the narrow - band signed distance - field.
    recursively subdivide the bounding box into eight sub-boxes
    */
    void NarrowBandSDF::ApplyOctreeFilter()
    {
        int depth = 0;
        for (int size = _block_width; size > 1; size /= 2)
            ++depth;
        // Skip the final subdivision; remaining nodes are evaluated afterward.
        if (depth > 0)
            --depth;
        if (depth == 0)
            return;

        // Adjacent blocks share boundary nodes. Process eight parity groups
        // in sequence so parallel blocks never write the same node.
        std::vector<std::size_t> groups[8];
        for (std::size_t index = 0; index < _blocks.size(); ++index)
        {
            const GridBlock& cubic = _blocks[index];
            if (!cubic.valid)
                continue;
            const int group =
                ((cubic.min_node.X() / _block_width) & 1) |
                (((cubic.min_node.Y() / _block_width) & 1) << 1) |
                (((cubic.min_node.Z() / _block_width) & 1) << 2);
            groups[group].push_back(index);
        }

#pragma omp parallel
        {
            for (int group = 0; group < 8; ++group)
            {
                const std::ptrdiff_t count =
                    static_cast<std::ptrdiff_t>(groups[group].size());
#pragma omp for schedule(guided, 1)
                for (std::ptrdiff_t index = 0; index < count; ++index)
                {
                    const GridBlock& cubic = _blocks[groups[group][index]];
                    SubdivideBlock(cubic.min_node, cubic.max_node, depth);
                }
                // The implicit barrier keeps neighboring groups apart.
            }
        }

        // Valid node distances are computed once after all filters finish.
    }

}
