// File       : turbulentKineticEnergyTransitionSSTEquation.h
// Created    : Mon Jan 14 2025
// Author     : Adam Fares
// Description: Turbulent kinetic energy equation for the transition SST model
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTKINETICENERGYTRANSITIONSSTEQUATION_H
#define TURBULENTKINETICENERGYTRANSITIONSSTEQUATION_H

#include "transitionShearStressTransportModel.h"
#include "turbulentKineticEnergyEquation.h"
#include "turbulentKineticEnergyTransitionSSTAssembler.h"

namespace accel
{

class turbulentKineticEnergyTransitionSSTEquation
    : public turbulentKineticEnergyEquation
{
private:
    transitionShearStressTransportModel* model_;

    using Assembler = turbulentKineticEnergyTransitionSSTAssembler;

public:
    turbulentKineticEnergyTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model);

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void solve() override;

    void applyOverrides() override;

protected:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYTRANSITIONSSTEQUATION_H
