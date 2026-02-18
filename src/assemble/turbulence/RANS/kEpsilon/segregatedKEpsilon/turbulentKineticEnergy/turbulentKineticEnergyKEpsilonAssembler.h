// File : turbulentKineticEnergyKEpsilonAssembler.h
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : Achraf Nagihi
// Description: Assembler for the turbulent kinetic energy equation in k-epsilon
// model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYKEPSILONASSEMBLER_H
#define TURBULENTKINETICENERGYKEPSILONASSEMBLER_H

#include "kEpsilonModel.h"
#include "turbulentKineticEnergyAssembler.h"

namespace accel
{

class turbulentKineticEnergyKEpsilonAssembler
    : public turbulentKineticEnergyAssembler
{
private:
    kEpsilonModel* model_;

public:
    turbulentKineticEnergyKEpsilonAssembler(kEpsilonModel* model)
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

#endif // TURBULENTKINETICENERGYKEPSILONASSEMBLER_H
