// File : turbulentEddyFrequencyEquation.h
// Created : Thu Mar 14 2024 12:50:04 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Turbulent eddy frequency equation for the SST model
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTEDDYFREQUENCYSSTEQUATION_H
#define TURBULENTEDDYFREQUENCYSSTEQUATION_H

#include "shearStressTransportModel.h"
#include "turbulentEddyFrequencyEquation.h"
#include "turbulentEddyFrequencySSTAssembler.h"

namespace accel
{

class turbulentEddyFrequencySSTEquation : public turbulentEddyFrequencyEquation
{
private:
    shearStressTransportModel* model_;

    using Assembler = turbulentEddyFrequencySSTAssembler;

public:
    turbulentEddyFrequencySSTEquation(realm* realm,
                                      shearStressTransportModel* model);

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTEDDYFREQUENCYSSTEQUATION_H
