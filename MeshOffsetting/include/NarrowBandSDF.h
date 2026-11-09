/**
 * @file NarrowBandFilter.h
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
    class NarrowBandSDF
    {
    public:
        NarrowBandSDF(const CMeshO& mesh,const Scalarm offsetDistance, int nodeCount = 512);
        ~NarrowBandSDF() = default;
    private:
		SSVH _ssvh;
        UniformGrid _grid;
        const std::size_t _blockSize = 8;
        void ApplyBoundingBoxFilter();
        void ApplySignedDistanceFilter();
        void ApplyOctreeFilter();
    };
}
#endif