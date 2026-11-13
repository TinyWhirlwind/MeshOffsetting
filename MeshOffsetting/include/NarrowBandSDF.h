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
#include "SSVH.h"
#include "UniformGrid.h"
namespace meshoffsetting
{
    enum class ClosestFeature
    {
        Face,
        EdgeAB,
        EdgeBC,
        EdgeCA,
        VertexA,
        VertexB,
        VertexC
    };

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

        //²âÊÔ½Ó¿Ú
        UniformGrid* GetGird();
    private:
        void Build();
		void CreateBlocks();
        ClosestFeature ClassifyBarycentric(Scalarm u, Scalarm v);
        void AngleWeightedPseudoNormal(CVertexO* v,Point3m& normal);
        void AngleWeightedPseudoNormal(CFaceO* f,int edge_id, Point3m& normal);
		void SubdivideBlock(const Point3i& min_node, const Point3i& max_node,int depth);
        int CheckSign(ClosestPointResult query_result, Point3m center, bool& is_valid);
        bool IsBoxValidBySignedDistance(Scalarm signed_distance, Scalarm circumsphere_radius) const;
        void SetNodeState(const Point3i& min_node, const Point3i& max_node, GridNodeState state);

        void ApplyBoundingBoxFilter();
        void ApplySignedDistanceFilter();
        void ApplyOctreeFilter();
    private:
        Scalarm _offset_distance = 0;
        Scalarm _cell_width = 0;
        CMeshO _mesh;

        const std::size_t _block_width = 0;
        Scalarm _block_circumsphere_radius = 0;
        std::vector<GridBlock> _blocks;
        
        SSVH _ssvh;
        UniformGrid _grid;
    };
}
#endif