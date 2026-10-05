// File       : rayleighPlessetCavitationMassTransfer.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "rayleighPlessetCavitationMassTransfer.h"

// std
#include <algorithm>
#include <cmath>

namespace accel
{

double
RayleighPlessetCavitationMassTransfer::rate(const massTransferState& s) const
{
    // guard against invalid densities (e.g. uninitialized fields)
    if (!validDensities(s))
        return 0.0;

    const double p = ratePressure(s.pressure);
    const double alphaV = clampedAlpha(s.alphaVapor);
    const double pv = saturationPressure();

    if (p < pv)
    {
        const double dp = std::max(0.0, pv - p);
        return cfg_.vaporizationCoefficient_ * 3.0 *
               cfg_.nucleationSiteVolumeFraction_ * (1.0 - alphaV) *
               s.rhoVapor / cfg_.nucleationSiteRadius_ *
               std::sqrt(2.0 / 3.0 * dp / s.rhoLiquid);
    }
    if (p > pv)
    {
        const double dp = std::max(0.0, p - pv);
        return -cfg_.condensationCoefficient_ * 3.0 * alphaV * s.rhoVapor /
               cfg_.nucleationSiteRadius_ *
               std::sqrt(2.0 / 3.0 * dp / s.rhoLiquid);
    }
    return 0.0;
}

} /* namespace accel */
