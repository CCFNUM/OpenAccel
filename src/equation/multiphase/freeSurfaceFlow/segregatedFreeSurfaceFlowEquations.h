// File       : segregatedFreeSurfaceFlowEquations.h
// Created    : Sun Jan 26 2025 22:53:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Segregated equation system for multiphase free-surface flow
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEGREGATEDFREESURFACEFLOWEQUATIONS_H
#define SEGREGATEDFREESURFACEFLOWEQUATIONS_H

// code
#include "bulkNavierStokesEquation.h"
#include "bulkPressureCorrectionEquation.h"
#include "freeSurfaceFlowModel.h"
#include "turbulenceModel.h"
#include "volumeFractionEquation.h"

namespace accel
{

class segregatedFreeSurfaceFlowEquations : public equation,
                                           public freeSurfaceFlowModel
{
public:
    static constexpr equationID ID = equationID::segregatedFreeSurface;

    segregatedFreeSurfaceFlowEquations(realm* realm);

    void addDomain(std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override
    {
    }

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

private:
    // instances for segregated equations
    std::unique_ptr<bulkNavierStokesEquation> U_eq_;

    std::unique_ptr<bulkPressureCorrectionEquation> pCorr_eq_;

    std::vector<std::unique_ptr<volumeFractionEquation>> alpha_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDFREESURFACEFLOWEQUATIONS_H
