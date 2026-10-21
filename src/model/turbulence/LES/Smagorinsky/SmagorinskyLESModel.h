// File       : SmagorinskyLESModel.h
// Created    : Fri Jul 31 2026
// Author     : GitHub Copilot
// Description: LES Smagorinsky turbulence model
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SMAGORINSKYLESMODEL_H
#define SMAGORINSKYLESMODEL_H

#include "turbulenceModel.h"

namespace accel
{

// Smagorinsky SGS model with Van Driest near-wall damping:
//
//   mut = rho * (Cs * delta * fD)^2 * |S|,  fD = 1 - exp(-yPlus / APlus)
//
// The coefficients Cs and APlus are owned per domain by `struct turbulence`
// in domain.h (keys `smagorinsky_constant` / `van_driest_a_plus`); they are
// deliberately not duplicated here.
class SmagorinskyLESModel : public turbulenceModel
{
public:
    SmagorinskyLESModel(realm* realm);

    void updateTurbulentDynamicViscosity(
        const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // SMAGORINSKYLESMODEL_H
