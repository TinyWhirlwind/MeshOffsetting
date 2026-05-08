#include "MeshOffseting.h"
MeshOffseting::MeshOffseting(const CMeshO& mesh, const Params& params)
    :_mesh(mesh),_params(params)
{
    vcg::tri::UpdateBounding<CMeshO>::Box(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(_mesh);

    _faces.clear();
    _faces.reserve(_mesh.fn);
    for (auto& f : _mesh.face)
    {
        _faces.emplace_back(&f);
    }
    _bvh.emplace(_faces);
    _grid.emplace(_mesh.bbox, _params.voxelSize, _params.offsetDistance);
}

MeshOffseting::~MeshOffseting() = default;

void MeshOffseting::SetMesh(const CMeshO& mesh)
{
    _mesh = mesh;
    vcg::tri::UpdateBounding<CMeshO>::Box(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(_mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(_mesh);

    _faces.clear();
    _faces.reserve(_mesh.fn);
    for (auto& f : _mesh.face)
    {
        _faces.emplace_back(&f);
    }
    _bvh.emplace(_faces);
}

void MeshOffseting::SetVoxelSize(float cellWidth)
{
    _params.voxelSize = cellWidth;

    if (!_mesh.bbox.IsNull())
    {
        _grid.emplace(_mesh.bbox,
            static_cast<float>(_params.voxelSize),
            static_cast<float>(_params.offsetDistance));
    }
}
void MeshOffseting::SetOffsetDistance(float offset)
{
    _params.offsetDistance = offset;
    if (!_mesh.bbox.IsNull())
    {
        _grid.emplace(_mesh.bbox,
            static_cast<float>(_params.voxelSize),
            static_cast<float>(_params.offsetDistance));
    }
}
void MeshOffseting::SetParams(const Params& params)
{
    _params = params;

    if (!_mesh.bbox.IsNull()&& !_faces.empty())
    {
        _grid.emplace(_mesh.bbox,
            static_cast<float>(_params.voxelSize),
            static_cast<float>(_params.offsetDistance));
    }
}

void MeshOffseting::Run()
{
    assert(_bvh && _grid);

    _grid->BuildBlocks();
    _grid->ForEachBlock([&](const Point3i& block) {
        const  OffsetGrid<float>::BlockData& bl = _grid->Block(block);
        Point3f sc = _grid->BlockCenter(bl);
        float r = _grid->BlockCircumsphereRadius(bl);
        //1.unsigned distance filter
        QueryResult qr;
        _bvh->QueryClosestPoint(sc, qr);
        float d = std::abs(_params.offsetDistance);
        float l = _params.voxelSize;
        if (qr.id != -1)
        {
            if (std::abs(qr.dist - d) > l + r)
            {
                _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
            }
            else
            {
                _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockRetained);
            }
        }
        //2.signed distance filter
        float signedDist = qr.sign * qr.dist;
        if (_grid->IsBlockRetained(block))
        {
            if (_params.offsetDistance > 0)
            {
                float lower = -d - l - r;
                float upper = std::min(-d + l + r, d - l - r);
                if (signedDist >= lower && signedDist <= upper)
                {
                    _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
                }
            }
            if (_params.offsetDistance < 0)
            {
                float lower = std::max(d - l - r, -d + l + r);
                float upper = d + l + r;
                if (signedDist >= lower && signedDist <= upper)
                {
                    _grid->SetBlockState(block, OffsetGrid<float>::BlockState::BlockInvalid);
                }
            }
        }
        });
    //3.octree filter

}