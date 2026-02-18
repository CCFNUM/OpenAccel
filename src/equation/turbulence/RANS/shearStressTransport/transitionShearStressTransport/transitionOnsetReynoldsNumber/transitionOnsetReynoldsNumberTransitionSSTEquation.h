// File : transitionOnsetReynoldsNumberTransitionSSTEquation.h
// Created : Tue Jan 14 2025
// Author : Adam Fares
// Description: Transition onset Reynolds number equation for the transition SST
// model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef TRANSITIONONSETREYNOLDSNUMBERTRANSITIONSSTEQUATION_H
#define TRANSITIONONSETREYNOLDSNUMBERTRANSITIONSSTEQUATION_H

#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "transitionOnsetReynoldsNumberTransitionSSTAssembler.h"
#include "transitionShearStressTransportModel.h"
#include "types.h"

namespace accel
{

class transitionOnsetReynoldsNumberTransitionSSTEquation
    : public equation,
      public linearSystem<1>
{
private:
    transitionShearStressTransportModel* model_;

    using Assembler = transitionOnsetReynoldsNumberTransitionSSTAssembler;

public:
    transitionOnsetReynoldsNumberTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    virtual void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void preTimeStep() override;

    void printScales() override;

protected:
    void setResidualScales_() override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TRANSITIONONSETREYNOLDSNUMBERTRANSITIONSSTEQUATION_H
