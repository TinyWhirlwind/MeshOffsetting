#include "SlicerSDF.h"
#include "Tolerance.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
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

        const std::size_t node_count =
            static_cast<std::size_t>(_size_x) *
            static_cast<std::size_t>(_size_y) *
            static_cast<std::size_t>(_size_z);
        _distance_cache.assign(node_count, Scalarm(0));
        _distance_cached.assign(node_count, false);
    }

    void SlicerSDF::StartSlicerOne(int layer, Scalarm layer_thickness)
    {
        if (layer < 0 ||
            !std::isfinite(layer_thickness) ||
            layer_thickness <= 0) {
            throw std::invalid_argument(
                "Invalid layer or layer thickness");
        }

        //_error_nodes.clear();
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
        //_error_nodes.clear();
        Scalarm max_z =_bbox.min.Z() + (_size_z - 1) * _cell_width;
        std::vector<Scalarm> slice_z_values;
        slice_z_values.reserve(static_cast<std::size_t>(layer_count));
        for (int layer = 0; layer < layer_count; layer++)
        {
            Scalarm zi = _bbox.min.Z() + (layer + 1) * layer_thickness;
            if (zi > max_z)
            {
                break;
            }
            slice_z_values.push_back(zi);
        }

        std::vector<SliceResult> results(slice_z_values.size());
        if (slice_z_values.empty())
        {
            _layers.swap(results);
            return;
        }

        if (!PrepareDistanceCache(slice_z_values))
        {
            for (std::size_t i = 0; i < results.size(); ++i)
            {
                results[i].z = slice_z_values[i];
                results[i].skipped_cells =
                    static_cast<std::size_t>(_size_x - 1) *
                    static_cast<std::size_t>(_size_y - 1);
            }
            _layers.swap(results);
            return;
        }

        std::vector<std::vector<Point3m>> layer_errors(
            slice_z_values.size());
        std::atomic<bool> failed(false);
        std::exception_ptr first_exception;

#pragma omp parallel
        {
            LayerWorkspace workspace;
            bool workspace_ready = true;
            try
            {
                InitializeWorkspace(workspace);
            }
            catch (...)
            {
#pragma omp critical(SlicerSDFException)
                {
                    if (!first_exception)
                    {
                        first_exception = std::current_exception();
                    }
                }
                failed.store(true);
                workspace_ready = false;
            }

#pragma omp for schedule(dynamic)
            for (int layer = 0;
                 layer < static_cast<int>(slice_z_values.size());
                 ++layer)
            {
                if (!workspace_ready || failed.load())
                {
                    continue;
                }

                try
                {
                    results[layer] = SliceAtZImpl(
                        slice_z_values[layer], workspace);
                    layer_errors[layer] = workspace.error_nodes;
                }
                catch (...)
                {
#pragma omp critical(SlicerSDFException)
                    {
                        if (!first_exception)
                        {
                            first_exception = std::current_exception();
                        }
                    }
                    failed.store(true);
                }
            }
        }

        if (first_exception)
        {
            std::rethrow_exception(first_exception);
        }

        for (const std::vector<Point3m>& errors : layer_errors)
        {
            //_error_nodes.insert(_error_nodes.end(), errors.begin(), errors.end());
        }
        _layers.swap(results);
    }

    const std::vector<SliceResult>& SlicerSDF::GetLayers() const
    {
        return _layers;
    }

    /*const std::vector<Point3m>& SlicerSDF::GetErrorNode()
    {
        return _error_nodes;
    }*/

    SliceResult SlicerSDF::SliceAtZ(Scalarm zi)
    {
        //_error_nodes.clear();
        SliceResult result;
        result.z = zi;

        std::vector<Scalarm> slice_z_values(1, zi);
        if (!PrepareDistanceCache(slice_z_values))
        {
            result.skipped_cells =
                static_cast<std::size_t>(_size_x - 1) *
                static_cast<std::size_t>(_size_y - 1);
            return result;
        }

        LayerWorkspace workspace;
        InitializeWorkspace(workspace);
        result = SliceAtZImpl(zi, workspace);
        //_error_nodes = workspace.error_nodes;
        return result;
    }

    void SlicerSDF::InitializeWorkspace(
        LayerWorkspace& workspace) const
    {
        workspace.plane_grid.values.assign(
            _size_x,
            std::vector<Scalarm>(_size_y, Scalarm(0)));
        workspace.plane_grid.valid.assign(
            _size_x,
            std::vector<bool>(_size_y, false));

        const std::size_t edge_count =
            static_cast<std::size_t>(_size_x - 1) * _size_y +
            static_cast<std::size_t>(_size_y - 1) * _size_x;
        const std::size_t initial_capacity =
            std::min(edge_count, std::size_t(16384));
        workspace.contour_nodes.reserve(initial_capacity);
    }

    bool SlicerSDF::ComputeSlicePosition(
        Scalarm zi,
        int& lower_z,
        Scalarm& weight,
        bool& needs_upper) const
    {
        Scalarm grid_z = (zi - _bbox.min.Z()) / _cell_width;
        Scalarm max_grid_z = static_cast<Scalarm>(_size_z - 1);
        if (!std::isfinite(grid_z) || grid_z < 0 || grid_z > max_grid_z)
        {
            return false;
        }

        lower_z = static_cast<int>(std::floor(grid_z));
        weight = grid_z - lower_z;
        needs_upper = weight > kScalarTolerance;
        return !needs_upper || lower_z + 1 < _size_z;
    }

    bool SlicerSDF::PrepareDistanceCache(
        const std::vector<Scalarm>& slice_z_values)
    {
        const std::size_t node_count =
            static_cast<std::size_t>(_size_x) *
            static_cast<std::size_t>(_size_y) *
            static_cast<std::size_t>(_size_z);
        std::vector<unsigned char> needed_nodes(node_count, 0);

        for (Scalarm zi : slice_z_values)
        {
            int lower_z = 0;
            Scalarm weight = 0;
            bool needs_upper = false;
            if (!ComputeSlicePosition(
                zi, lower_z, weight, needs_upper))
            {
                return false;
            }

            for (int x = 0; x + 1 < _size_x; x++)
            {
                for (int y = 0; y + 1 < _size_y; y++)
                {
                    bool candidate_cell = false;
                    for (int dx = 0; dx <= 1; dx++)
                    {
                        for (int dy = 0; dy <= 1; dy++)
                        {
                            GridNode& lower = _grid->Node(
                                Point3i(x + dx, y + dy, lower_z));
                            if (lower.state == GridNodeState::NodeDistanceReady)
                            {
                                candidate_cell = true;
                            }

                            if (needs_upper)
                            {
                                GridNode& upper = _grid->Node(
                                    Point3i(x + dx, y + dy, lower_z + 1));
                                if (upper.state == GridNodeState::NodeDistanceReady)
                                {
                                    candidate_cell = true;
                                }
                            }
                        }
                    }

                    if (!candidate_cell)
                    {
                        continue;
                    }

                    for (int dx = 0; dx <= 1; dx++)
                    {
                        for (int dy = 0; dy <= 1; dy++)
                        {
                            Point3i lower_id(x + dx, y + dy, lower_z);
                            needed_nodes[_grid->NodeIndex(lower_id)] = 1;

                            if (needs_upper)
                            {
                                Point3i upper_id(
                                    x + dx, y + dy, lower_z + 1);
                                needed_nodes[_grid->NodeIndex(upper_id)] = 1;
                            }
                        }
                    }
                }
            }
        }

        std::vector<std::size_t> node_indices;
        for (std::size_t index = 0; index < needed_nodes.size(); index++)
        {
            if (needed_nodes[index])
            {
                node_indices.push_back(index);
            }
        }

        const std::ptrdiff_t task_count =
            static_cast<std::ptrdiff_t>(node_indices.size());
#pragma omp parallel for schedule(guided, 64)
        for (std::ptrdiff_t task = 0; task < task_count; ++task)
        {
            const std::size_t index =
                node_indices[static_cast<std::size_t>(task)];
            Point3i node = _grid->NodeCoord(index);

            Scalarm distance = 0;
            GetGridNodeDistance(
                node.X(), node.Y(), node.Z(), distance);
        }

        return true;
    }

    SliceResult SlicerSDF::SliceAtZImpl(
        Scalarm zi,
        LayerWorkspace& workspace) const
    {
        SliceResult result;
        result.z = zi;
        workspace.contour_nodes.clear();
        workspace.crossing_ids.clear();
        workspace.error_nodes.clear();

        if (!BuildPlaneField(zi, workspace))
        {
            result.skipped_cells =
                static_cast<std::size_t>(_size_x - 1) *
                static_cast<std::size_t>(_size_y - 1);
            return result;
        }

        for (int x = 0; x + 1 < _size_x; x++)
        {
            for (int y = 0; y + 1 < _size_y; y++)
            {
                ProcessCell(x, y, result, workspace);
            }
        }
        TraceContours(result, workspace);
        return result;
    }

    bool SlicerSDF::BuildPlaneField(
        Scalarm zi,
        LayerWorkspace& workspace) const
    {
        for (int x = 0; x < _size_x; x++)
        {
            for (int y = 0; y < _size_y; y++)
            {
                workspace.plane_grid.values[x][y] = 0;
                workspace.plane_grid.valid[x][y] = false;
            }
        }

        int lower_z = 0;
        Scalarm weight = 0;
        bool needs_upper = false;
        if (!ComputeSlicePosition(
            zi, lower_z, weight, needs_upper))
        {
            return false;
        }

        for (int x = 0; x < _size_x; x++)
        {
            for (int y = 0; y < _size_y; y++)
            {
                FillPlaneNode(
                    x, y, lower_z, weight, needs_upper, workspace);
            }
        }

        return true;
    }

    bool SlicerSDF::GetGridNodeDistance(
        int x,
        int y,
        int z,
        Scalarm& distance)
    {
        Point3i node_id(x, y, z);
        const std::size_t index = _grid->NodeIndex(node_id);
        if (_distance_cached[index])
        {
            distance = _distance_cache[index];
            return true;
        }

        GridNode& node = _grid->Node(node_id);
        if (node.state == GridNodeState::NodeDistanceReady &&
            std::isfinite(node.signed_distance))
        {
            distance = node.signed_distance;
        }
        else if (!_nsdf.QuerySignedDistance(
                     _grid->NodePosition(node_id), distance) ||
                 !std::isfinite(distance))
        {
            return false;
        }

        _distance_cache[index] = distance;
        _distance_cached[index] = true;
        return true;
    }

    bool SlicerSDF::GetCachedGridNodeDistance(
        int x,
        int y,
        int z,
        Scalarm& distance) const
    {
        const std::size_t index =
            static_cast<std::size_t>(x) +
            static_cast<std::size_t>(y) * _size_x +
            static_cast<std::size_t>(z) * _size_x * _size_y;
        if (!_distance_cached[index])
        {
            return false;
        }

        distance = _distance_cache[index];
        return std::isfinite(distance);
    }

    bool SlicerSDF::FillPlaneNode(
        int x,
        int y,
        int lower_z,
        Scalarm weight,
        bool needs_upper,
        LayerWorkspace& workspace) const
    {
        if (workspace.plane_grid.valid[x][y])
        {
            return true;
        }

        Scalarm lower_distance = 0;
        if (!GetCachedGridNodeDistance(
            x, y, lower_z, lower_distance))
        {
            return false;
        }

        Scalarm distance = lower_distance;
        if (needs_upper)
        {
            Scalarm upper_distance = 0;
            if (!GetCachedGridNodeDistance(
                x, y, lower_z + 1, upper_distance))
            {
                return false;
            }

            distance = lower_distance * (Scalarm(1) - weight)
                + upper_distance * weight;
        }

        workspace.plane_grid.values[x][y] = distance - _offset;
        workspace.plane_grid.valid[x][y] = true;
        return true;
    }

    void SlicerSDF::ProcessCell(
        int x,
        int y,
        SliceResult& result,
        LayerWorkspace& workspace) const
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
            if (workspace.plane_grid.valid[xs[k]][ys[k]])
            {
                values[k] = workspace.plane_grid.values[xs[k]][ys[k]];
                valid_count++;
            }
        }

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

            inside[k] = values[k] < 0;
            if (inside[k])
            {
                case_id |= (1 << k);
            }
        }

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
                result.z,
                workspace);
        }

        switch (case_id)
        {
        case 1:
        case 14:
            Connect(edge_nodes[3], edge_nodes[0], workspace);
            break;
        case 2:
        case 13:
            Connect(edge_nodes[0], edge_nodes[1], workspace);
            break;
        case 3:
        case 12:
            Connect(edge_nodes[1], edge_nodes[3], workspace);
            break;
        case 4:
        case 11:
            Connect(edge_nodes[1], edge_nodes[2], workspace);
            break;
        case 6:
        case 9:
            Connect(edge_nodes[0], edge_nodes[2], workspace);
            break;
        case 7:
        case 8:
            Connect(edge_nodes[2], edge_nodes[3], workspace);
            break;
        case 5:
        case 10:
        {
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
            if (q >= 0)
            {
                Connect(edge_nodes[0], edge_nodes[1], workspace);
                Connect(edge_nodes[2], edge_nodes[3], workspace);
            }
            else
            {
                Connect(edge_nodes[0], edge_nodes[3], workspace);
                Connect(edge_nodes[1], edge_nodes[2], workspace);
            }
            break;
        }
        default:
            break;
        }
    }

    std::size_t SlicerSDF::GetOrCreateCrossing(
        int ax,
        int ay,
        int bx,
        int by,
        Scalarm zi,
        LayerWorkspace& workspace) const
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
        found = workspace.crossing_ids.find(key);
        if (found != workspace.crossing_ids.end()) {
            return found->second;
        }
        Scalarm value_a = workspace.plane_grid.values[ax][ay];
        Scalarm value_b = workspace.plane_grid.values[bx][by];

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
        std::size_t node_id = workspace.contour_nodes.size();
        workspace.contour_nodes.push_back(node);
        workspace.crossing_ids[key] = node_id;
        return node_id;
    }

    void SlicerSDF::Connect(
        std::size_t a,
        std::size_t b,
        LayerWorkspace& workspace) const
    {
        workspace.contour_nodes[a].neighbours.push_back(b);
        workspace.contour_nodes[b].neighbours.push_back(a);
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

    void SlicerSDF::TraceContours(
        SliceResult& result,
        LayerWorkspace& workspace) const
    {
        bool is_error = false;
        // 每个交点必须有两个邻居。
        for (std::size_t i = 0; i < workspace.contour_nodes.size(); i++)
        {
            if (workspace.contour_nodes[i].neighbours.size() != 2)
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
            for (std::size_t i = 0; i < workspace.contour_nodes.size(); i++)
            {
                if (workspace.contour_nodes[i].neighbours.size() != 2)
                {
                    workspace.error_nodes.push_back(
                        workspace.contour_nodes[i].point);
                }
            }
            return;
        }

        std::vector<bool> visited(workspace.contour_nodes.size(), false);
        for (std::size_t start = 0;
            start < workspace.contour_nodes.size();
            start++)
        {
            if (visited[start])
            {
                continue;
            }
            std::vector<Point3m> contour;
            std::size_t current = start;
            std::size_t previous = workspace.contour_nodes.size();
            do {
                if (visited[current])
                {
                    throw std::runtime_error("Contour traversal failed");
                }
                visited[current] = true;
                const Point3m& point = workspace.contour_nodes[current].point;
                // 零值网格顶点可能让多个交点坐标重合。
                // 先保留图连接，输出时去掉连续重复坐标。
                if (contour.empty() ||!SamePoint(contour.back(), point))
                {
                    contour.push_back(point);
                }

                const std::vector<std::size_t>& neighbours =
                    workspace.contour_nodes[current].neighbours;
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
