// File       : thermalEnergyEquation.cpp
// Created    : Mon Dec 01 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Thermal energy equation implementation details
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALENERGYEQUATION_H
#define THERMALENERGYEQUATION_H

#ifndef WITH_THERMAL_TEMPERATURE

// code
#include "equation.h"
#include "heatTransferModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "thermalEnergyAssembler.h"
#include "types.h"

namespace accel
{

class thermalEnergyEquation : public equation,
                              public heatTransferModel,
                              public linearSystem<1>
{
    using Assembler = thermalEnergyAssembler;

public:
    static constexpr equationID ID = equationID::thermalEnergy;

    thermalEnergyEquation(realm* realm);

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif

#endif // THERMALENERGYEQUATION_H
