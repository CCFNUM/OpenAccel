// File : bulkPressureCorrectionEquation.h
// Created : Sun Feb 02 2025 20:44:56 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Pressure correction equation for bulk free-surface flow
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BULKPRESSURECORRECTIONEQUATION_H
#define BULKPRESSURECORRECTIONEQUATION_H

// code
#include "bulkPressureCorrectionAssembler.h"
#include "equation.h"
#include "freeSurfaceFlowModel.h"
#include "linearSystem.h"

namespace accel
{

class bulkPressureCorrectionEquation : public equation, public linearSystem<1>
{
private:
    freeSurfaceFlowModel* model_;

    using Assembler = bulkPressureCorrectionAssembler;

public:
    bulkPressureCorrectionEquation(realm* realm, freeSurfaceFlowModel* model);

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

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // BULKPRESSURECORRECTIONEQUATION_H
