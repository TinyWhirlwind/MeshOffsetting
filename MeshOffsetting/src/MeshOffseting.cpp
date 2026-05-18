#include "MeshOffseting.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include <vcg/space/index/octree_template.h>

namespace
{
    template <class ScalarType>
    class BlockOctree : public vcg::OctreeTemplate<int, ScalarType>
    {
    public:
        using Base = vcg::OctreeTemplate<int, ScalarType>;
        using BoundingBoxType = typename Base::BoundingBoxType;
        using NodePointer = typename Base::NodePointer;

        void Init(const BoundingBoxType& bbox, int maxDepth)
        {
            this->boundingBox = bbox;
            this->Initialize(maxDepth);
        }

        NodePointer EnsureChild(NodePointer parent, int sonIndex)
        {
            NodePointer& child = this->Son(parent, sonIndex);
            if (child == nullptr)
            {
                child = this->NewNode(parent, sonIndex);
            }
            return child;
        }
    };

    inline int CeilLog2Positive(int value)
    {
        int depth = 0;
        int span = 1;
        while (span < value)
        {
            span <<= 1;//span = span * 2;
            ++depth;
        }
        return depth;
    }

    template <class ScalarType>
    bool InClosedInterval(ScalarType value, ScalarType lower, ScalarType upper)
    {
        return lower <= upper && value >= lower && value <= upper;
    }
}

MeshOffseting::MeshOffseting(const CMeshO& mesh, const Params& params)
    : _mesh(mesh), _params(params)
{
    vcg::tri::UpdateBounding<CMeshO>::Box(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(_mesh);

    RebuildBVH();
    RebuildGrid();
}

MeshOffseting::~MeshOffseting() = default;

void MeshOffseting::SetMesh(const CMeshO& mesh)
{
    _mesh = mesh;
    vcg::tri::UpdateBounding<CMeshO>::Box(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(_mesh);

    RebuildBVH();
    RebuildGrid();
}

void MeshOffseting::SetVoxelSize(float cellWidth)
{
    _params.voxelSize = cellWidth;
    RebuildGrid();
}

void MeshOffseting::SetOffsetDistance(float offset)
{
    _params.offsetDistance = offset;
    RebuildGrid();
}

void MeshOffseting::SetParams(const Params& params)
{
    _params = params;
    RebuildGrid();
}

void MeshOffseting::RebuildBVH()
{
    _faces.clear();
    _faces.reserve(_mesh.fn);
    for (auto& f : _mesh.face)
    {
        _faces.push_back(&f);
    }

    if (_faces.empty())
    {
        _bvh.reset();
        return;
    }

    _bvh.emplace(_faces);
}

void MeshOffseting::RebuildGrid()
{
    if (_mesh.bbox.IsNull() || _params.voxelSize <= 0.0f)
    {
        _grid.reset();
        return;
    }

    _grid.emplace(_mesh.bbox, _params.voxelSize, _params.offsetDistance);
}

bool MeshOffseting::QuerySignedDistance(const Point3m& p, float signDist) const
{
    assert(_bvh);

    QueryResult qr;
    _bvh->QueryClosestPoint(p, qr);
    if (qr.id == -1)
    {
        return false;
    }
    signDist = static_cast<float>(qr.sign) * static_cast<float>(qr.dist);
    return true;
}

bool MeshOffseting::IsInvalidByUnsignedDistance(float unsignedDistance, float radius) const
{
    const float d = std::abs(_params.offsetDistance);
    const float l = _params.voxelSize;
    return std::abs(unsignedDistance - d) > l + radius;
}

bool MeshOffseting::IsInvalidBySignedDistance(float signedDistance, float radius) const
{
    const float d = std::abs(_params.offsetDistance);
    const float l = _params.voxelSize;

    if (_params.offsetDistance > 0.0f)
    {
        const float lower = -d - l - radius;
        const float upper = std::min(-d + l + radius, d - l - radius);
        return InClosedInterval(signedDistance, lower, upper);
    }

    if (_params.offsetDistance < 0.0f)
    {
        const float lower = std::max(d - l - radius, -d + l + radius);
        const float upper = d + l + radius;
        return InClosedInterval(signedDistance, lower, upper);
    }

    return false;
}

void MeshOffseting::MarkNodesInWorldBoxInvalid(const Box3m& box)
{
    assert(_grid);

    const Point3i nodeDims = _grid->NodeDims();
    const Point3m origin = _grid->NodePosition(Point3i(0, 0, 0));
    const float cellWidth = _grid->CellWidth();
    const float eps = cellWidth * 1e-4f;

    Point3i minNode(0, 0, 0);
    Point3i maxNode(0, 0, 0);
    for (int axis = 0; axis < 3; ++axis)
    {
        minNode[axis] = std::max(
            0,
            static_cast<int>(std::ceil((box.min[axis] - origin[axis]) / cellWidth - eps)));
        maxNode[axis] = std::min(
            nodeDims[axis] - 1,
            static_cast<int>(std::floor((box.max[axis] - origin[axis]) / cellWidth + eps)));
    }

    if (minNode[0] > maxNode[0] || minNode[1] > maxNode[1] || minNode[2] > maxNode[2])
    {
        return;
    }

    for (int i = minNode[0]; i <= maxNode[0]; ++i)
    {
        for (int j = minNode[1]; j <= maxNode[1]; ++j)
        {
            for (int k = minNode[2]; k <= maxNode[2]; ++k)
            {
                _grid->SetNodeInvalid(Point3i(i, j, k));
            }
        }
    }
}

void MeshOffseting::ApplySignedDistanceFilter()
{
    assert(_bvh && _grid);

    _grid->ForEachBlock([&](const Point3i& block)
        {
            const OffsetGrid<float>::BlockData& bl = _grid->Block(block);
            const Point3m sc = _grid->BlockCenter(bl);
            const float r = _grid->BlockCircumsphereRadius(bl);
            const float d = std::abs(_params.offsetDistance);
            const float l = _params.voxelSize;

            QueryResult qr;
            _bvh->QueryClosestPoint(sc, qr);
            if (qr.id == -1)
            {
                return;
            }

            if (std::abs(static_cast<float>(qr.dist) - d) > l + r)
            {
                _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
                _grid->ForEachNodeInBlock(block, [&](const Point3i& nodeCoord)
                    {
                        _grid->SetNodeInvalid(nodeCoord);
                    });
                return;
            }

            _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockRetained);

            const float signedDist = static_cast<float>(qr.sign) * static_cast<float>(qr.dist);
            if (IsInvalidBySignedDistance(signedDist, r))
            {
                _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
                _grid->ForEachNodeInBlock(block, [&](const Point3i& nodeCoord)
                    {
                        _grid->SetNodeInvalid(nodeCoord);
                    });
            }
            else
            {
                _grid->ForEachNodeInBlock(block, [&](const Point3i& nodeCoord)
                    {
                        const Point3m nodePos = _grid->NodePosition(nodeCoord);
                        float nodeSignedDist = 0.0;
                        if (!QuerySignedDistance(nodePos, nodeSignedDist))
                        {
                            _grid->SetNodeInvalid(nodeCoord);
                        }
                        else
                        {
                            if (IsInvalidBySignedDistance(nodeSignedDist, 0.0f))
                            {
                                _grid->SetNodeInvalid(nodeCoord);
                            }
                            else
                            {
                                _grid->SetSignedDistance(nodeCoord, nodeSignedDist);
                            }
                        }
                    });
            }
        });
}

void MeshOffseting::ApplyOctreeFilter()
{
    assert(_bvh && _grid);

    const float cellWidth = _grid->CellWidth();
    _grid->ForEachBlock([&](const Point3i& bp)
        {
            if (!_grid->IsBlockRetained(bp))
            {
                return;
            }

            const auto& block = _grid->Block(bp);
            const int maxCellSpan = std::max({ block.cellSpan[0], block.cellSpan[1], block.cellSpan[2] });
            const int maxDepth = CeilLog2Positive(maxCellSpan);
            if (maxDepth <= 0)
            {
                return;
            }

            BlockOctree<Scalarm> octree;
            Box3m blockBox;
            blockBox.min = _grid->NodePosition(block.minCell);
            blockBox.max = _grid->NodePosition(block.minCell + block.cellSpan);
            octree.Init(blockBox, maxDepth);

            auto recurse = [&](auto&& self, typename BlockOctree<Scalarm>::NodePointer node) -> void
                {
                    Box3m nodeBox;
                    octree.BoundingBoxInWorldCoordinates(node, nodeBox);

                    const float signedDistance = QuerySignedDistance(nodeBox.Center());
                    const float radius = static_cast<float>(nodeBox.Dim().Norm() * 0.5);
                    if (IsInvalidBySignedDistance(signedDistance, radius))
                    {
                        MarkNodesInWorldBoxInvalid(nodeBox);
                        return;
                    }

                    const Point3m dim = nodeBox.Dim();
                    if (dim[0] <= cellWidth && dim[1] <= cellWidth && dim[2] <= cellWidth)
                    {
                        return;
                    }

                    for (int sonIndex = 0; sonIndex < 8; ++sonIndex)
                    {
                        self(self, octree.EnsureChild(node, sonIndex));
                    }
                };

            for (int sonIndex = 0; sonIndex < 8; ++sonIndex)
            {
                recurse(recurse, octree.EnsureChild(octree.Root(), sonIndex));
            }
        });

}

void MeshOffseting::BuildIntersectEdges()
{
    assert(_bvh && _grid);

    _grid->ClearIntersectEdge();
    _grid->BuildIntersectEdges();
}

void GetIntersectsPoint()
{
     
}

void MeshOffseting::Run()
{
    assert(_bvh && _grid);

    _grid->ResetNodes();
    _grid->ResetBlocks();
    _grid->BuildBlocks();

    ApplySignedDistanceFilter();
    ApplyOctreeFilter();
}
