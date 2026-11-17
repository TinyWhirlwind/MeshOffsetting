#include "SlicerSDF.h"
#include "Tolerance.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace meshoffsetting
{
    SlicerSDF::SlicerSDF(const CMeshO& mesh, Scalarm offset, int nodeCount, int blockWidth)
        : _nsdf(mesh, offset, nodeCount, blockWidth),
        _offset(offset)
    {
        _grid = _nsdf.GetGird();
        if (_grid == nullptr)
        {
            throw std::runtime_error("SDF grid is null");
        }

        _bbox = _grid->GetBox();
        _cell_width = _grid->GetCellWidth();
        Point3i size = _grid->GetSize();
        _size_x = size.X();
        _size_y = size.Y();
        _size_z = size.Z();
        if (_size_x < 2 || _size_y < 2 || _size_z < 2)
        {
            throw std::invalid_argument(
                "Each grid axis needs at least two nodes");
        }

        if (!std::isfinite(_cell_width)
            || _cell_width <= 0
            || !std::isfinite(_offset))
        {
            throw std::invalid_argument(
                "Invalid cell width or offset");
        }

        _plane_grid.values.assign(
            _size_x,
            std::vector<Scalarm>(_size_y, Scalarm(0)));
        _plane_grid.valid.assign(
            _size_x,
            std::vector<bool>(_size_y, false));
    }

    void SlicerSDF::StartSlicerOne(int layer, Scalarm layer_thickness)
    {
        if (layer < 0 ||
            !std::isfinite(layer_thickness) ||
            layer_thickness <= 0) {
            throw std::invalid_argument(
                "Invalid layer or layer thickness");
        }

        Scalarm z = _bbox.min.Z() + (layer + 1.0) * layer_thickness;

        std::vector<SliceResult> results;
        results.push_back(SliceAtZ(static_cast<Scalarm>(z)));
        _layers.swap(results);
    }

    void SlicerSDF::StartSlicerAll(int layer_count, Scalarm layer_thickness)
    {
        if (layer_count < 0 ||!std::isfinite(layer_thickness) || layer_thickness <= 0) {
            throw std::invalid_argument(
                "Invalid layer count or layer thickness");
        }
        Scalarm max_z =_bbox.min.Z() + (_size_z - 1) * _cell_width;
        std::vector<SliceResult> results;
        for (int layer = 0; layer < layer_count; layer++)
        {
            Scalarm zi = _bbox.min.Z() + (layer + 1) * layer_thickness;
            if (zi > max_z)
            {
                break;
            }
            results.push_back(SliceAtZ(zi));
        }
        _layers.swap(results);
    }

    const std::vector<SliceResult>& SlicerSDF::GetLayers() const
    {
        return _layers;
    }

    const std::vector<Point3m>& SlicerSDF::GetErrorNode()
    {
        return _error_nodes;
    }

    SliceResult SlicerSDF::SliceAtZ(Scalarm zi)
    {
        SliceResult result;
        result.z = zi;
        _contour_nodes.clear();
        _crossing_ids.clear();
        BuildPlaneField(zi);
        for (int x = 0; x + 1 < _size_x; x++)
        {
            for (int y = 0; y + 1 < _size_y; y++)
            {
                ProcessCell(x, y, result);
            }
        }
        TraceContours(result);
        return result;
    }

    void SlicerSDF::BuildPlaneField(Scalarm zi)
    {
        for (int x = 0; x < _size_x; x++)
        {
            for (int y = 0; y < _size_y; y++)
            {
                _plane_grid.values[x][y] = 0;
                _plane_grid.valid[x][y] = false;
            }
        }

        Scalarm grid_z = (zi - _bbox.min.Z()) / _cell_width;
        Scalarm max_grid_z = static_cast<Scalarm>(_size_z - 1);
        if (!std::isfinite(grid_z) || grid_z < 0 || grid_z > max_grid_z)
        {
            return;
        }

        int lower_z = static_cast<int>(std::floor(grid_z));
        Scalarm weight = grid_z - lower_z;
        bool needs_upper = weight > kScalarTolerance;
        if (needs_upper && lower_z + 1 >= _size_z)
        {
            return;
        }

        for (int x = 0; x + 1 < _size_x; x++)
        {
            for (int y = 0; y + 1 < _size_y; y++)
            {
                // 论文窄带中只要这个体素已有一个有效节点，
                // 它就可能是切片需要检查的候选体素。
                if (!HasStoredNodeInCell(
                    x, y, lower_z, needs_upper))
                {
                    continue;
                }

                FillPlaneNode(x, y, lower_z, weight, needs_upper);
                FillPlaneNode(x + 1, y, lower_z, weight, needs_upper);
                FillPlaneNode(x + 1, y + 1, lower_z, weight, needs_upper);
                FillPlaneNode(x, y + 1, lower_z, weight, needs_upper);
            }
        }
    }

    bool SlicerSDF::HasStoredNodeInCell(
        int x,
        int y,
        int lower_z,
        bool needs_upper)
    {
        for (int dx = 0; dx <= 1; dx++)
        {
            for (int dy = 0; dy <= 1; dy++)
            {
                GridNode& lower =
                    _grid->Node(Point3i(x + dx, y + dy, lower_z));
                if (lower.state == GridNodeState::NodeDistanceReady)
                {
                    return true;
                }

                if (needs_upper)
                {
                    GridNode& upper = _grid->Node(
                        Point3i(x + dx, y + dy, lower_z + 1));
                    if (upper.state == GridNodeState::NodeDistanceReady)
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool SlicerSDF::FillPlaneNode(
        int x,
        int y,
        int lower_z,
        Scalarm weight,
        bool needs_upper)
    {
        if (_plane_grid.valid[x][y])
        {
            return true;
        }

        Point3i lower_id(x, y, lower_z);
        GridNode& lower = _grid->Node(lower_id);
        Scalarm lower_distance = 0;
        if (lower.state == GridNodeState::NodeDistanceReady)
        {
            lower_distance = lower.signed_distance;
        }
        else if (!_nsdf.QuerySignedDistance(_grid->NodePosition(lower_id), lower_distance))
        {
            return false;
        }

        Scalarm distance = lower_distance;
        if (needs_upper)
        {
            Point3i upper_id(x, y, lower_z + 1);
            GridNode& upper = _grid->Node(upper_id);
            Scalarm upper_distance = 0;
            if (upper.state == GridNodeState::NodeDistanceReady)
            {
                upper_distance = upper.signed_distance;
            }
            else if (!_nsdf.QuerySignedDistance(
                _grid->NodePosition(upper_id), upper_distance))
            {
                return false;
            }

            distance = lower_distance * (Scalarm(1) - weight)
                + upper_distance * weight;
        }

        _plane_grid.values[x][y] = distance - _offset;
        _plane_grid.valid[x][y] = true;
        return true;
    }

    void SlicerSDF::ProcessCell( int x,int y, SliceResult& result)
    {
        /*
            顶点编号：
            0(x,y) -------- 3(x,y+1)
            |               |
            |               |
            1(x+1,y) ------ 2(x+1,y+1)
            边 0：0 -> 1
            边 1：1 -> 2
            边 2：2 -> 3
            边 3：3 -> 0
        */
        int xs[4] = { x, x + 1, x + 1, x };
        int ys[4] = { y, y, y + 1, y + 1 };
        Scalarm values[4] = {};
        int valid_count = 0;
        for (int k = 0; k < 4; k++)
        {
            if (_plane_grid.valid[xs[k]][ys[k]])
            {
                values[k] = _plane_grid.values[xs[k]][ys[k]];
                valid_count++;
            }
        }

        // 四个角点都有标量值，才能确定16种情况中的一种。
        // 这不是Marching Squares的一种case，而是输入场是否完整的检查。
        if (valid_count != 4)
        {
            if (valid_count != 0)
            {
                result.skipped_cells++;
            }
            return;
        }

        bool inside[4] = {};
        int case_id = 0;
        for (int k = 0; k < 4; k++)
        {
            if (!std::isfinite(values[k]))
            {
                result.skipped_cells++;
                return;
            }

            // 统一把非常接近0的值放到外侧，避免同一节点在相邻cell中分类不同。
            if (std::abs(values[k]) <= kScalarTolerance)
            {
                values[k] = kScalarTolerance;
            }

            inside[k] = values[k] < 0;
            if (inside[k])
            {
                case_id |= (1 << k);
            }
        }

        // case 0和15没有交线。
        if (case_id == 0 || case_id == 15)
        {
            return;
        }

        std::size_t edge_nodes[4] = {};
        for (int edge = 0; edge < 4; edge++)
        {
            int a = edge;
            int b = (edge + 1) % 4;
            if (inside[a] == inside[b])
            {
                continue;
            }

            edge_nodes[edge] = GetOrCreateCrossing(
                xs[a], ys[a],
                xs[b], ys[b],
                result.z);
        }

        switch (case_id)
        {
        case 1:
        case 14:
            Connect(edge_nodes[3], edge_nodes[0]);
            break;
        case 2:
        case 13:
            Connect(edge_nodes[0], edge_nodes[1]);
            break;
        case 3:
        case 12:
            Connect(edge_nodes[1], edge_nodes[3]);
            break;
        case 4:
        case 11:
            Connect(edge_nodes[1], edge_nodes[2]);
            break;
        case 6:
        case 9:
            Connect(edge_nodes[0], edge_nodes[2]);
            break;
        case 7:
        case 8:
            Connect(edge_nodes[2], edge_nodes[3]);
            break;
        case 5:
        case 10:
        {
            // case 5和10有四个交点，使用双线性插值的渐近判别。
            Scalarm scale = 0;
            for (int k = 0; k < 4; k++)
            {
                scale = std::max(scale, std::abs(values[k]));
            }
            Scalarm f0 = values[0] / scale;
            Scalarm f1 = values[1] / scale;
            Scalarm f2 = values[2] / scale;
            Scalarm f3 = values[3] / scale;
            Scalarm q = f0 * f2 - f1 * f3;
            // q == 0 时固定选择第一种配对。
            if (q >= 0)
            {
                Connect(edge_nodes[0], edge_nodes[1]);
                Connect(edge_nodes[2], edge_nodes[3]);
            }
            else
            {
                Connect(edge_nodes[0], edge_nodes[3]);
                Connect(edge_nodes[1], edge_nodes[2]);
            }
            break;
        }
        default:
            break;
        }
    }

    std::size_t SlicerSDF::GetOrCreateCrossing(int ax,int ay,int bx,int by,Scalarm zi)
    {
        std::size_t id_a = static_cast<std::size_t>(ax) * _size_y + ay;
        std::size_t id_b = static_cast<std::size_t>(bx) * _size_y + by;
        std::pair<std::size_t, std::size_t> key;
        if (id_a < id_b)
        {
            key = std::make_pair(id_a, id_b);
        }
        else
        {
            key = std::make_pair(id_b, id_a);
        }
        std::map<std::pair<std::size_t, std::size_t>,std::size_t>::iterator found;
        found = _crossing_ids.find(key);
        if (found != _crossing_ids.end()) {
            return found->second;
        }
        Scalarm value_a = _plane_grid.values[ax][ay];
        Scalarm value_b = _plane_grid.values[bx][by];

        // 与ProcessCell使用同一个零值规则。
        if (std::abs(value_a) <= kScalarTolerance)
        {
            value_a = kScalarTolerance;
        }
        if (std::abs(value_b) <= kScalarTolerance)
        {
            value_b = kScalarTolerance;
        }

        // 同比例缩放不改变插值比例，避免差值溢出。
        Scalarm scale = std::max(std::abs(value_a), std::abs(value_b));
        value_a = value_a / scale;
        value_b = value_b / scale;
        Scalarm denominator = value_a - value_b;
        Scalarm t = Scalarm(0.5);
        if (std::abs(denominator) > kScalarTolerance)
        {
            t = value_a / denominator;
        }
        t = std::clamp(t, Scalarm(0), Scalarm(1));
        Scalarm world_x = _bbox.min.X()
            + (ax + (bx - ax) * t) * _cell_width;

        Scalarm world_y =_bbox.min.Y()
            + (ay + (by - ay) * t) * _cell_width;

        ContourNode node;
        node.point = Point3m(
            static_cast<Scalarm>(world_x),
            static_cast<Scalarm>(world_y),
            zi);
        std::size_t node_id = _contour_nodes.size();
        _contour_nodes.push_back(node);
        _crossing_ids[key] = node_id;
        return node_id;
    }

    void SlicerSDF::Connect(std::size_t a, std::size_t b)
    {
        _contour_nodes[a].neighbours.push_back(b);
        _contour_nodes[b].neighbours.push_back(a);
    }

    bool SlicerSDF::SamePoint(const Point3m& a, const Point3m& b)
    {
        return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z();
    }

    bool SlicerSDF::HasArea(const std::vector<Point3m>& contour)
    {
        if (contour.size() < 3)
        {
            return false;
        }

        Scalarm origin_x = contour[0].X();
        Scalarm origin_y = contour[0].Y();
        Scalarm twice_area = 0;
        for (std::size_t i = 0; i < contour.size(); i++)
        {
            std::size_t next = i + 1;
            if (next == contour.size())
            {
                next = 0;
            }
            Scalarm ax = contour[i].X() - origin_x;
            Scalarm ay = contour[i].Y() - origin_y;
            Scalarm bx = contour[next].X() - origin_x;
            Scalarm by = contour[next].Y() - origin_y;
            twice_area += ax * by - ay * bx;
        }

        return twice_area != 0;
    }

    void SlicerSDF::TraceContours(SliceResult& result)
    {
        bool is_error = false;
        // 每个交点必须有两个邻居。
        for (std::size_t i = 0; i < _contour_nodes.size(); i++)
        {
            if (_contour_nodes[i].neighbours.size() != 2)
            {
               /* throw std::runtime_error(
                    "Contour is open or branches; "
                    "check grid bounds and SDF validity");*/
                is_error = true;
                break;
            }
        }
        if (is_error)
        {
            for (std::size_t i = 0; i < _contour_nodes.size(); i++)
            {
                _error_nodes.push_back(_contour_nodes[i].point);
            }
            return;
        }

        std::vector<bool> visited(_contour_nodes.size(), false);
        for (std::size_t start = 0; start < _contour_nodes.size();start++)
        {
            if (visited[start])
            {
                continue;
            }
            std::vector<Point3m> contour;
            std::size_t current = start;
            std::size_t previous = _contour_nodes.size();
            do {
                if (visited[current])
                {
                    throw std::runtime_error("Contour traversal failed");
                }
                visited[current] = true;
                const Point3m& point = _contour_nodes[current].point;
                // 零值网格顶点可能让多个交点坐标重合。
                // 先保留图连接，输出时去掉连续重复坐标。
                if (contour.empty() ||!SamePoint(contour.back(), point))
                {
                    contour.push_back(point);
                }

                const std::vector<std::size_t>& neighbours =
                    _contour_nodes[current].neighbours;
                std::size_t next = neighbours[0];
                if (next == previous)
                {
                    next = neighbours[1];
                }
                previous = current;
                current = next;

            } while (current != start);

            if (contour.size() > 1 &&
                SamePoint(contour.front(), contour.back()))
            {
                contour.pop_back();
            }

            // 去掉退化的点环、线环
            if (!HasArea(contour))
            {
                continue;
            }
            result.contours.push_back(contour);
        }
    }
}
