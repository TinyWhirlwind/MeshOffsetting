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
    vcg::tri::UpdateTopology<CMeshO>::FaceFace(_mesh);
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

bool MeshOffseting::QuerySignedDistance(const Point3m& p, QueryResult& qr) const
{
    assert(_bvh);
    _bvh->QueryClosestPoint(p, qr);
    if (qr._id == -1)
    {
        return false;
    }
    //signDist = static_cast<float>(qr.sign) * static_cast<float>(qr.dist);
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
            if (qr._id == -1)
            {
                return;
            }

            if (std::abs(static_cast<float>(qr._dist) - d) > l + r)
            {
                _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
                _grid->ForEachNodeInBlock(block, [&](const Point3i& nodeCoord)
                    {
                        _grid->SetNodeInvalid(nodeCoord);
                    });
                return;
            }

            _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockRetained);
            float signedDist = static_cast<float>(qr._sign) * static_cast<float>(qr._dist);
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
                        QueryResult qr;
                        if (!QuerySignedDistance(nodePos, qr))
                        {
                            _grid->SetNodeInvalid(nodeCoord);
                        }
                        else
                        {
                            float nodeSignedDist = static_cast<float>(qr._sign) * static_cast<float>(qr._dist);
                            if (IsInvalidBySignedDistance(nodeSignedDist, 0.0f))
                            {
                                _grid->SetNodeInvalid(nodeCoord);
                            }
                            else
                            {
                                _grid->SetSignedDistance(nodeCoord, qr);
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

                    QueryResult qr;
                    if (!QuerySignedDistance(nodeBox.Center(), qr))
                        return;
                    float signedDistance = static_cast<float>(qr._sign) * static_cast<float>(qr._dist);
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

void MeshOffseting::FindIntersectionPointOnGridEdge()
{
    assert(_bvh && _grid);
    std::vector<OffsetGrid<float>::GridEdge>& interEdges = _grid->GetIntersectEdges();
    for (auto& ed : interEdges)
    {
        assert(ed._p0 == ed._p1);
        const OffsetGrid<float>::NodeData&  n0 = _grid->Node(ed._p0);
        const OffsetGrid<float>::NodeData&  n1 = _grid->Node(ed._p1);

        CFaceO* sameFace = FindSameTriangle(n0, n1);
        TriangleRegion r0 = ToTriangleRegion(n0.query);
        TriangleRegion r1 = ToTriangleRegion(n1.query);
        if (!sameFace)
        {
            //bisection search
            SolveByBisection(ed);
        }
        else
        {
            //analytical solution
            SolveByAnalyticalSolution(ed);
        }
    }
}

CFaceO* MeshOffseting::FindSameTriangle(const OffsetGrid<float>::NodeData& n0, const OffsetGrid<float>::NodeData& n1)
{
    QueryResult q0 = n0.query;
    QueryResult q1 = n1.query;
    assert(q0._id != -1 && q1._id != -1);
    if (q0._id == q1._id)
        return &_mesh.face[q0._id];
    CFaceO& f0 = _mesh.face[q0._id];
    CFaceO& f1 = _mesh.face[q1._id];
    for (int i = 0; i < 3; ++i)
    {
        if (f0.cFFp(i) != &f1)continue;
        if (q0._closestType == ClosestType::Vertex && (q0._closestFeature == i || q0._closestFeature == (i + 1) % 3))
        {
            return &f1;
        }
        if (q0._closestType == ClosestType::Edge && q0._closestFeature == i)
        {
            return &f1;
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        if (f1.cFFp(i) != &f0)continue;
        if (q1._closestType == ClosestType::Vertex && (q1._closestFeature == i || q1._closestFeature == (i + 1) % 3))
        {
            return &f0;
        }
        if (q1._closestType == ClosestType::Edge && q1._closestFeature == i)
        {
            return &f0;
        }
    }
    int share = 0;
    int sharePid = -1;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (f0.V(i)->Index() == f1.V(j)->Index())
            {
                share++;
                sharePid = f0.V(i)->Index();
            }
        }
    }
    if (share == 1)
    {
        if (q0._closestType == ClosestType::Vertex && f0.V(q0._closestFeature)->Index() == sharePid)
        {
            return &f1;
        }
        if (q1._closestType == ClosestType::Vertex && f1.V(q1._closestFeature)->Index() == sharePid)
        {
            return &f0;
        }
        if (q0._closestType == ClosestType::Edge && q1._closestType == ClosestType::Edge)
        {
            auto es0 = f0.V(q0._closestFeature)->Index();
            auto ee0 = f0.V((q0._closestFeature + 1) % 3)->Index();
            auto es1 = f1.V(q1._closestFeature)->Index();
            auto ee1 = f1.V((q1._closestFeature + 1) % 3)->Index();

            if (ee0 == es1)
            {
                return f0.cFFp(q0._closestFeature);
            }

            if (ee1 == es0)
            {
                return f1.cFFp(q1._closestFeature);
            }
        }
    }
    return nullptr;
} 

TriangleRegion MeshOffseting::ToTriangleRegion(const QueryResult& qr)
{
    if (qr._closestType == ClosestType::Face)
    {
        return TriangleRegion::Face;
    }

    if (qr._closestType == ClosestType::Edge)
    {
        if (qr._closestFeature == 0) return TriangleRegion::Edge01;
        if (qr._closestFeature == 1) return TriangleRegion::Edge02;
        if (qr._closestFeature == 2) return TriangleRegion::Edge12;
    }
    if (qr._closestType == ClosestType::Vertex)
    {
        if (qr._closestFeature == 0) return TriangleRegion::Vertex0;
        if (qr._closestFeature == 1) return TriangleRegion::Vertex1;
        if (qr._closestFeature == 2) return TriangleRegion::Vertex2;
    }

    return TriangleRegion::Unknown;
}

void MeshOffseting::SolveByBisection(OffsetGrid<float>::GridEdge& edge)
{
}

void MeshOffseting::SolveByAnalyticalSolution(OffsetGrid<float>::GridEdge& edge)
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
