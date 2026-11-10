#include "CMeshO.h"
#include "NarrowBandSDF.h"
namespace meshoffsetting
{
    NarrowBandSDF::NarrowBandSDF(const CMeshO& mesh, const Scalarm offsetDistance,int nodeCount, int blockWidth):
        _mesh(mesh),
        _offsetDistance(offsetDistance),
        _block_width(blockWidth),
        _ssvh(mesh),
		_grid(mesh.bbox, offsetDistance, nodeCount)
    {
        _cell_width = _grid.GetCellWidth();
    }


    int NarrowBandSDF::CheckSign(ClosestPointResult query_result, Point3m center, bool& is_valid)
    {
        is_valid = true;
        if (query_result.faceId == -1)
        {
            is_valid = false;
            return 0;
        }
        auto& f = _mesh.face[query_result.faceId];
        int sign = (center - query_result.closestPoint) * f.N() > 0 ? 1 : -1;
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
                    if (!is_valid)continue;

                    Scalarm signed_distance = sign * closest_point.unsigned_distance;
                    if (!IsBoxValidBySignedDistance(signed_distance, box_circumsphere_radius))
                    {
                        //设置无效

                        continue;
                    }
                    SubdivideBlock(child_min, child_max, --depth);
                }
            }
        }
    }

    bool NarrowBandSDF::IsBoxValidBySignedDistance(Scalarm signed_distance, Scalarm circumsphere_radius) const
    {
        Scalarm d = std::abs(_offsetDistance);
        if (_offsetDistance > 0)
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
                    if (state == GridNodeState::NodeDistanceReady)
                    {
                        Point3m node_pos = _grid.NodePosition(node_id);
                        ClosestPointResult closest_point = _ssvh.QueryClosestPoint(node_pos);
                        bool is_valid = true;
                        int sign = CheckSign(closest_point, node_pos, is_valid);
                        if (is_valid)
                        {
                            cur_node.signed_distance = sign * closest_point.unsigned_distance;
                        }
                    }
                    cur_node.state = GridNodeState::NodeInvalid;
                }
            }
        }
    }

    void NarrowBandSDF::Build()
    {
        if (_offsetDistance <= kScalarTolerance)
            return;
        CreateBlocks();
        ApplyBoundingBoxFilter();
    }

    void NarrowBandSDF::CreateBlocks()
    {
        Point3m size;
        Point3m voxel;
        _grid.GetSizeAndVoxel(size, voxel);

        int count_x = std::ceil((size.X() - 1) / _block_width);
        int count_y = std::ceil((size.Y() - 1) / _block_width);
        int count_z = std::ceil((size.Z() - 1) / _block_width);
        _blocks.clear();
        _blocks.resize(count_x * count_y * count_z);
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
		_block_circumsphere_radius = std::sqrt(3.0) * voxel.X() * _block_width * 0.5;
    }
    
    //Bounding Box Filter
    void NarrowBandSDF::ApplyBoundingBoxFilter()
    {
        for (auto& cubic : _blocks)
        {
            ClosestPointResult closest_point = _ssvh.QueryClosestPoint(cubic.center);
            ////用block中心点到网格的无符号距离场剔除不在窄带内的block
            if (std::abs(closest_point.unsigned_distance - _offsetDistance) > _block_circumsphere_radius + _cell_width)
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
            if (!cubic.valid)continue;
            int sign = CheckSign(cubic.closestPointResult, cubic.center, cubic.valid);
            if (!cubic.valid)continue;//不是重复判断
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
        for (const GridBlock& block : _blocks)
        {
            if (!block.valid)
                continue;
			//递归细分
            SubdivideBlock(block.min_node, block.max_node, depth);
        }
    }

}