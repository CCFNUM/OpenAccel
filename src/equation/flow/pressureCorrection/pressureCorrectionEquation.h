// File : pressureCorrectionEquation.h
// Created : Thu Mar 14 2024 12:50:04 (+0100)
// Author : Fabian Wermelinger
// Description: Pressure correction (continuity) equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESSURECORRECTIONEQUATION_H
#define PRESSURECORRECTIONEQUATION_H

// code
#include "equation.h"
#include "flowModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "pressureCorrectionAssembler.h"

namespace accel
{

class pressureCorrectionEquation : public equation, public linearSystem<1>
{
private:
    flowModel* model_;

    using Assembler = pressureCorrectionAssembler;

public:
    pressureCorrectionEquation(realm* realm, flowModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

protected:
    void setResidualScales_() override;

    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // PRESSURECORRECTIONEQUATION_H
