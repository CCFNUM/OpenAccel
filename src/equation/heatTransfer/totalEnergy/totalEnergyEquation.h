// File : totalEnergyEquation.h
// Created : Thu Mar 27 2025 10:42:10 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Total energy conservation equation for heat transfer
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TOTALENERGYEQUATION_H
#define TOTALENERGYEQUATION_H

// code
#include "equation.h"
#include "heatTransferModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "totalEnergyAssembler.h"
#include "types.h"

namespace accel
{

class totalEnergyEquation : public equation,
                            public heatTransferModel,
                            public linearSystem<1>
{
    using Assembler = totalEnergyAssembler;

public:
    static constexpr equationID ID = equationID::totalEnergy;

    totalEnergyEquation(realm* realm);

    void checkDomain(const std::shared_ptr<domain> domain) override;

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

#endif // TOTALENERGYEQUATION_H
