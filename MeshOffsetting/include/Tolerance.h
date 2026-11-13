/**
 * @file Tolerance.h
 * @author zhangpingping
 * @date 2026-11-06
 * @brief
 * »›≤Ó…Ë÷√
 */
#ifndef MESHOFFSETTING_TOLERANCE_H
#define MESHOFFSETTING_TOLERANCE_H
#include "base_type.h"
namespace meshoffsetting
{
    inline constexpr Scalarm kScalarTolerance = std::numeric_limits<Scalarm>::epsilon() * Scalarm(64);
    inline constexpr Scalarm kScalarToleranceSquared = kScalarTolerance * kScalarTolerance;
    inline constexpr Scalarm Distance = Scalarm(1e-6);
    inline constexpr Scalarm SquaredDistance = Distance * Distance;
    inline constexpr Scalarm Parameter = Scalarm(1e-6);
}
#endif