/**
 * @file NarrowBandSDF.h
 * @author zhangpingping
 * @date 2026-11-07
 * @brief
 *
 */
#ifndef MESHOFFSETTING_NARROW_BAND_FILTER_H
#define MESHOFFSETTING_NARROW_BAND_FILTER_H
#include "CMeshO.h"
#include "FCPWClosestPointQuery.h"
#include "UniformGrid.h"
namespace meshoffsetting
{

    class NarrowBandSDF
    {
        struct GridBlock
        {
            Point3i min_node;
            Point3i max_node;
            Point3m center;
            bool valid = true;
            ClosestPointResult closestPointResult;
        };

    public:
        NarrowBandSDF(const CMeshO& mesh,const Scalarm offsetDistance, int nodeCount = 129, int blockWidth = 8);
        ~NarrowBandSDF() = default;
        NarrowBandSDF(const NarrowBandSDF&) = delete;
        NarrowBandSDF& operator=(const NarrowBandSDF&) = delete;
        NarrowBandSDF(NarrowBandSDF&&) = delete;
        NarrowBandSDF& operator=(NarrowBandSDF&&) = delete;

        Point3i GetGridSize() const;
        Box3m GetGridBox() const;
        Scalarm GetCellWidth() const;

        bool TryGetNodeDistance(
            const Point3i& node_id,
            Scalarm& signed_distance) const;
        bool EvaluateNodeDistance(
            const Point3i& node_id,
            Scalarm& signed_distance) const;
        bool QuerySignedDistance(
            const Point3m& query_point,
            Scalarm& signed_distance) const;
    private:
        void Build();
        void CreateBlocks();
        void InitializeActiveNodes();
        void AngleWeightedPseudoNormal(CVertexO* v,Point3m& normal) const;
        void AngleWeightedPseudoNormal(const CFaceO* f,int edge_id, Point3m& normal) const;
        void SubdivideBlock(const Point3i& min_node, const Point3i& max_node,int depth);
        int CheckSign(const ClosestPointResult& query_result,const Point3m& center,bool& is_valid) const;
        bool IsBoxValidBySignedDistance(Scalarm signed_distance, Scalarm circumsphere_radius) const;
        void MarkNodeRangeInvalid(const Point3i& min_node, const Point3i& max_node);

        void ApplyBoundingBoxFilter();
        void ApplySignedDistanceFilter();
        void ApplyOctreeFilter();
        void ComputeRemainingNodeDistances();
    private:
        Scalarm _offset_distance = 0;
        Scalarm _cell_width = 0;
        CMeshO _mesh;

        const std::size_t _block_width = 0;
        Scalarm _block_circumsphere_radius = 0;
        std::vector<GridBlock> _blocks;

        FCPWClosestPointQuery _closest_point_query;
        UniformGrid _grid;
    };
}
#endif
