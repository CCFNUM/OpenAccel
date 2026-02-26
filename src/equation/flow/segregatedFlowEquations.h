// File       : segregatedFlowEquations.h
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Fabian Wermelinger
// Description: Navier-Stokes convenience wrapper to solve coupled momentum and
// pressure correction separately. Alternatively, these two
// equations can be solved by defining them explicitly in the YAML
// config file of the run. This class is just a wrapper for that
// since the two equations should be solved together.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEGREGATEDFLOWEQUATIONS_H
#define SEGREGATEDFLOWEQUATIONS_H

// code
#include "flowModel.h"
#include "navierStokesEquation.h"
#include "pressureCorrectionEquation.h"

namespace accel
{

class segregatedFlowEquations : public equation, public flowModel
{
public:
    static constexpr equationID ID = equationID::segregatedFlow;

    segregatedFlowEquations(realm* realm);

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
    std::unique_ptr<navierStokesEquation> U_eq_;

    std::unique_ptr<pressureCorrectionEquation> pCorr_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDFLOWEQUATIONS_H
