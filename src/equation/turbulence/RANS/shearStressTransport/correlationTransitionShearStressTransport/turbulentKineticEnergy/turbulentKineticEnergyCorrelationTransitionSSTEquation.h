// File       : turbulentKineticEnergyCorrelationTransitionSSTEquation.h
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: TKE equation for correlation-based transition SST (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTEQUATION_H
#define TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTEQUATION_H

#include "correlationTransitionShearStressTransportModel.h"
#include "turbulentKineticEnergyCorrelationTransitionSSTAssembler.h"
#include "turbulentKineticEnergyEquation.h"

namespace accel
{

class turbulentKineticEnergyCorrelationTransitionSSTEquation
    : public turbulentKineticEnergyEquation
{
private:
    correlationTransitionShearStressTransportModel* model_;

    using Assembler = turbulentKineticEnergyCorrelationTransitionSSTAssembler;

public:
    turbulentKineticEnergyCorrelationTransitionSSTEquation(
        realm* realm,
        correlationTransitionShearStressTransportModel* model);

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void solve() override;

protected:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTEQUATION_H
