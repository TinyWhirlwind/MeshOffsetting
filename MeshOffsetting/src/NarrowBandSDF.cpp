#include "CMeshO.h"
#include "Tolerance.h"
#include "NarrowBandSDF.h"
namespace meshoffsetting
{
    NarrowBandSDF::NarrowBandSDF(const CMeshO& mesh, const Scalarm offsetDistance,int nodeCount, int blockWidth) :
        _offset_distance(offsetDistance),
        _mesh(mesh),
        _block_width(blockWidth)
    {
        _mesh.updateBoxAndNormals();
        _ssvh = SSVH(_mesh);
        _grid = UniformGrid(_mesh.bbox,offsetDistance,blockWidth,nodeCount);
        _cell_width = _grid.GetCellWidth();

        _mesh.vert.EnableVFAdjacency();
        _mesh.face.EnableVFAdjacency();
        _mesh.face.EnableFFAdjacency();
        vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
        vcg::tri::UpdateTopology<CMeshO>::VertexFace(_mesh);
        vcg::tri::UpdateTopology<CMeshO>::FaceFace(_mesh);
        Build();
    }

    UniformGrid* NarrowBandSDF::GetGird()
    {
        return &_grid;
    }

    ClosestFeature NarrowBandSDF::ClassifyBarycentric(Scalarm u, Scalarm v)
    {
        // U:dir(B-A) V:dir(C-A) 
        // query_p = A + u(B-A) + v(C-A) = (1 - u - v)A + uB + vC = wA + uB + vC
        Scalarm w = 1 - u - v;
        bool zu = std::abs(u) < kScalarTolerance;
        bool zv = std::abs(v) < kScalarTolerance;
        bool zw = std::abs(w) < kScalarTolerance;
        if (zu && zv) return ClosestFeature::VertexA;
        if (zu && zw) return ClosestFeature::VertexC;
        if (zv && zw) return ClosestFeature::VertexB;
        
        if (zu) return ClosestFeature::EdgeCA;
        if (zv) return ClosestFeature::EdgeAB;
        if (zw) return ClosestFeature::EdgeBC;

        return ClosestFeature::Face;
    }

    void NarrowBandSDF::AngleWeightedPseudoNormal(CVertexO* v, Point3m& normal)
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

    void NarrowBandSDF::AngleWeightedPseudoNormal(CFaceO* f, int edge_id, Point3m& normal)
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

    int NarrowBandSDF::CheckSign(ClosestPointResult query_result, Point3m center, bool& is_valid)
    {
        is_valid = true;
        if (query_result.faceId == -1)
        {
            is_valid = false;
            return 0;
        }
        //Angle-weighted pseudonormal
        //根据查询点到网格的最近点，根据最近点的位置（点、边、面），采用三种法向，判断查询点在模型内部还是外部
        int sign = 0;
        auto* f = &(_mesh.face[query_result.faceId]);
        Point3m nf(0., 0., 0.);
        ClosestFeature type = ClassifyBarycentric(query_result.u, query_result.v);
        switch (type)
        {
        case ClosestFeature::VertexA:
            AngleWeightedPseudoNormal(f->V(0), nf);
            break;
        case ClosestFeature::VertexB:
            AngleWeightedPseudoNormal(f->V(1), nf);
            break;
        case ClosestFeature::VertexC:
            AngleWeightedPseudoNormal(f->V(2), nf);
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
        default:
            break;
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
                    ClosestPointResult closest_point = _ssvh.QueryClosestPoint(box_center);
                    bool is_valid = true;
                    int sign = CheckSign(closest_point, box_center, is_valid);
                    if (!is_valid)
                    {
                        SetNodeState(child_min, child_max, GridNodeState::NodeInvalid);
                        continue;
                    }

                    Scalarm signed_distance = sign * closest_point.unsigned_distance;
                    if (std::abs(signed_distance - _offset_distance) > _cell_width + box_circumsphere_radius)
                    {
                        //设置 Invalid grid nodes
                        SetNodeState(child_min, child_max, GridNodeState::NodeInvalid);
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

    void NarrowBandSDF::SetNodeState(const Point3i& min_node, const Point3i& max_node, GridNodeState state)
    {
        for (int x = min_node.X(); x <= max_node.X(); ++x)
        {
            for (int y = min_node.Y(); y <= max_node.Y(); ++y)
            {
                for (int z = min_node.Z(); z <= max_node.Z(); ++z)
                {
                    Point3i node_id(x, y, z);
                    GridNode& cur_node = _grid.Node(node_id);
                    //共享节点只要属于任意一个无效子块，就必须是 NodeInvalid
                    if (cur_node.state == GridNodeState::Unknown && state == GridNodeState::NodeDistanceReady)
                    {
                        Point3m node_pos = _grid.NodePosition(node_id);
                        ClosestPointResult closest_point = _ssvh.QueryClosestPoint(node_pos);
                        bool is_valid = true;
                        int sign = CheckSign(closest_point, node_pos, is_valid);
                        if (is_valid)
                        {
                            cur_node.state = GridNodeState::NodeDistanceReady;
                            cur_node.signed_distance = sign * closest_point.unsigned_distance;
                        }
                        else
                        {
                            cur_node.state = GridNodeState::NodeInvalid;
                        }
                    }
                    else if(cur_node.state == GridNodeState::NodeDistanceReady && state == GridNodeState::NodeDistanceReady)
                    {
                        continue;
                    }
                    else
                    {
                        cur_node.state = GridNodeState::NodeInvalid;
                    }
                }
            }
        }
    }

    void NarrowBandSDF::Build()
    {
        if (std::abs(_offset_distance) <= kScalarTolerance)
            return;
        CreateBlocks();
        ApplyBoundingBoxFilter();
        ApplySignedDistanceFilter();
        ApplyOctreeFilter(); 
    }

    void NarrowBandSDF::CreateBlocks()
    {
        Point3m size = _grid.GetSize();

        int count_x = std::ceil((size.X() - 1) / _block_width);
        int count_y = std::ceil((size.Y() - 1) / _block_width);
        int count_z = std::ceil((size.Z() - 1) / _block_width);
        _blocks.clear();
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
        for (auto& cubic : _blocks)
        {
            ClosestPointResult closest_point = _ssvh.QueryClosestPoint(cubic.center);
            //计算block中心点到网格的无符号距离场剔除不在窄带内的block
            if (std::abs(closest_point.unsigned_distance - std::abs(_offset_distance)) > _block_circumsphere_radius + _cell_width)
            {
                cubic.valid = false;
                SetNodeState(cubic.min_node, cubic.max_node, GridNodeState::NodeInvalid);
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
            if (!cubic.valid)continue;
            int sign = CheckSign(cubic.closestPointResult, cubic.center, cubic.valid);
            if (!cubic.valid)continue;//不是重复判断
			Scalarm signed_distance = sign * cubic.closestPointResult.unsigned_distance;
            if (!IsBoxValidBySignedDistance(signed_distance, _block_circumsphere_radius))
            {
                cubic.valid = false;
                SetNodeState(cubic.min_node, cubic.max_node, GridNodeState::NodeInvalid);
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
        for (const GridBlock& cubic : _blocks)
        {
            if (!cubic.valid)
                continue;
			//递归细分,设置 Invalid grid nodes
            SubdivideBlock(cubic.min_node, cubic.max_node, depth);
        }

        //设置 valid grid nodes
        for (const GridBlock& cubic : _blocks)
        {
            if (!cubic.valid)
                continue;
            SetNodeState(cubic.min_node, cubic.max_node, GridNodeState::NodeDistanceReady);
        }
    }

}