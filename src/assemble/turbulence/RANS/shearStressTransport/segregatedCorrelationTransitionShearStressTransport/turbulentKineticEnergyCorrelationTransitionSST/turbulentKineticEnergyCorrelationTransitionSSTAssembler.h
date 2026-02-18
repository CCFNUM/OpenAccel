// File : turbulentKineticEnergyCorrelationTransitionSSTAssembler.h
// Created : Sun Dec 29 2024
// Author : Adam Fares
// Description: Assembler for TKE equation in correlation-based
// transition SST model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTASSEMBLER_H
#define TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTASSEMBLER_H

#include "correlationTransitionShearStressTransportModel.h"
#include "turbulentKineticEnergyAssembler.h"

namespace accel
{

class turbulentKineticEnergyCorrelationTransitionSSTAssembler
    : public turbulentKineticEnergyAssembler
{
private:
    correlationTransitionShearStressTransportModel* model_;

public:
    turbulentKineticEnergyCorrelationTransitionSSTAssembler(
        correlationTransitionShearStressTransportModel* model)
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

#endif // TURBULENTKINETICENERGYCORRELATIONTRANSITIONSSTASSEMBLER_H
