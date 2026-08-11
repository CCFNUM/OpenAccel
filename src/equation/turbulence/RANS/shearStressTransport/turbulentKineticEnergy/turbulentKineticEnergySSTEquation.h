// File       : turbulentKineticEnergyEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Turbulent kinetic energy equation for the SST model
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTKINETICENERGYSSTEQUATION_H
#define TURBULENTKINETICENERGYSSTEQUATION_H

#include "shearStressTransportModel.h"
#include "turbulentKineticEnergyEquation.h"
#include "turbulentKineticEnergySSTAssembler.h"

namespace accel
{

class turbulentKineticEnergySSTEquation : public turbulentKineticEnergyEquation
{
private:
    shearStressTransportModel* model_;

    using Assembler = turbulentKineticEnergySSTAssembler;

public:
    turbulentKineticEnergySSTEquation(realm* realm,
                                      shearStressTransportModel* model);

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void solve() override;

    void applyOverrides() override;

protected:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYSSTEQUATION_H
