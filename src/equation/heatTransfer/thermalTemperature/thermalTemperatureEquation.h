// File       : thermalTemperatureEquation.h
// Created    : Thu Feb 22 2024 08:42:10 (+0100)
// Author     : Fabian Wermelinger
// Description: System to solve scalar temperature equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALTEMPERATUREEQUATION_H
#define THERMALTEMPERATUREEQUATION_H

#ifdef WITH_THERMAL_TEMPERATURE

// code
#include "equation.h"
#include "heatTransferModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "thermalTemperatureAssembler.h"
#include "types.h"

namespace accel
{

class thermalTemperatureEquation : public equation,
                                   public heatTransferModel,
                                   public linearSystem<1>
{
    using Assembler = thermalTemperatureAssembler;

public:
    static constexpr equationID ID = equationID::thermalTemperature;

    thermalTemperatureEquation(realm* realm);

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

#endif // THERMALTEMPERATUREEQUATION_H
