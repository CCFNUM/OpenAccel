// File : turbulentKineticEnergyAssembler.h
// Created : Thu Feb 22 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Abstract base assembler for the turbulent kinetic energy
// equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTKINETICENERGYASSEMBLER_H
#define TURBULENTKINETICENERGYASSEMBLER_H

#include "phiAssembler.h"
#include "turbulenceModel.h"

namespace accel
{

class turbulentKineticEnergyAssembler : public phiAssembler<1>
{
public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    turbulentKineticEnergyAssembler(turbulenceModel* model) : Base(model)
    {
    }

protected:
    virtual void assembleNodeTermsFusedSteady_(const domain* domain,
                                               Context* ctx) override = 0;
    virtual void
    assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                              Context* ctx) override = 0;
    virtual void
    assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                               Context* ctx) override = 0;

    void postAssemble_(const domain* domain, Context* ctx) override
    {
        Base::postAssemble_(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    }
};

} /* namespace accel */

#endif // TURBULENTKINETICENERGYASSEMBLER_H
