// File       : turbulentKineticEnergySSTAssembler.h
// Created    : Thu Feb 22 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the turbulent kinetic energy equation in SST model
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTKINETICENERGYSSTASSEMBLER_H
#define TURBULENTKINETICENERGYSSTASSEMBLER_H

#include "shearStressTransportModel.h"
#include "turbulentKineticEnergyAssembler.h"

namespace accel
{

class turbulentKineticEnergySSTAssembler
    : public turbulentKineticEnergyAssembler
{
private:
    shearStressTransportModel* model_;

public:
    turbulentKineticEnergySSTAssembler(shearStressTransportModel* model)
        : turbulentKineticEnergyAssembler(model), model_(model)
    {
    }

protected:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;
    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;
    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    // Boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYSSTASSEMBLER_H
