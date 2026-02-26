// File       : turbulentKineticEnergyKEpsilonEquation.h
// Created    : Thu Feb 22 2025 13:38:51 (+0100)
// Author     : Achraf Nagihi
// Description: Turbulent kinetic energy equation for the k-epsilon model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYKEPSILONEQUATION_H
#define TURBULENTKINETICENERGYKEPSILONEQUATION_H

#include "kEpsilonModel.h"
#include "turbulentKineticEnergyEquation.h"
#include "turbulentKineticEnergyKEpsilonAssembler.h"

namespace accel
{

class turbulentKineticEnergyKEpsilonEquation
    : public turbulentKineticEnergyEquation
{
private:
    kEpsilonModel* model_;

    using Assembler = turbulentKineticEnergyKEpsilonAssembler;

public:
    turbulentKineticEnergyKEpsilonEquation(realm* realm, kEpsilonModel* model);

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void solve() override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYKEPSILONEQUATION_H
