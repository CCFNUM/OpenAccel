// File       : turbulentIntermittencyTransitionSSTEquation.h
// Created    : Tue Jan 14 2025
// Author     : Adam Fares
// Description: Turbulent intermittency equation for the transition SST model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTINTERMITTENCYTRANSITIONSSTEQUATION_H
#define TURBULENTINTERMITTENCYTRANSITIONSSTEQUATION_H

#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "transitionShearStressTransportModel.h"
#include "turbulentIntermittencyTransitionSSTAssembler.h"
#include "types.h"

namespace accel
{

class turbulentIntermittencyTransitionSSTEquation : public equation,
                                                    public linearSystem<1>
{
private:
    transitionShearStressTransportModel* model_;

    using Assembler = turbulentIntermittencyTransitionSSTAssembler;

public:
    turbulentIntermittencyTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void preTimeStep() override;

    void printScales() override;

protected:
    void setResidualScales_() override;

    void applyDependencyUpdates_(const domain* domain,
                                 const stk::mesh::EntityRank entityRank,
                                 STKScalarField& stk_dst) override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTINTERMITTENCYTRANSITIONSSTEQUATION_H
