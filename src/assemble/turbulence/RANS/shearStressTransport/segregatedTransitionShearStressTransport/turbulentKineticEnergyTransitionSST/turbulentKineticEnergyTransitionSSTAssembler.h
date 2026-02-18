// File : turbulentKineticEnergyTransitionSSTAssembler.h
// Created : Wed Jan 15 2025
// Author : Adam Fares
// Description: Assembler for the turbulent kinetic energy in transition SST
// model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYTRANSITIONSSTASSEMBLER_H
#define TURBULENTKINETICENERGYTRANSITIONSSTASSEMBLER_H

#include "transitionShearStressTransportModel.h"
#include "turbulentKineticEnergyAssembler.h"

namespace accel
{

class turbulentKineticEnergyTransitionSSTAssembler
    : public turbulentKineticEnergyAssembler
{
private:
    transitionShearStressTransportModel* model_;

public:
    turbulentKineticEnergyTransitionSSTAssembler(
        transitionShearStressTransportModel* model)
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

#endif // TURBULENTKINETICENERGYTRANSITIONSSTASSEMBLER_H
