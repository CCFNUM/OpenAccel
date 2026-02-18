// File : turbulentIntermittencyTransitionSSTAssembler.h
// Created : Wed Jan 15 2025
// Author : Adam Fares
// Description: Assembler for the intermittency equation in the transition SST
// model
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTINTERMITTENCYTRANSITIONSSTASSEMBLER_H
#define TURBULENTINTERMITTENCYTRANSITIONSSTASSEMBLER_H

// code
#include "phiAssembler.h"
#include "transitionShearStressTransportModel.h"

#include <cassert>

namespace accel
{

class turbulentIntermittencyTransitionSSTAssembler : public phiAssembler<1>
{
private:
    transitionShearStressTransportModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    turbulentIntermittencyTransitionSSTAssembler(
        transitionShearStressTransportModel* model)
        : Base(model), model_(model)
    {
    }

protected:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;
    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;
    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    void postAssemble_(const domain* domain, Context* ctx) override
    {
        Base::postAssemble_(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    }
};

} /* namespace accel */

#endif // TURBULENTINTERMITTENCYTRANSITIONSSTASSEMBLER_H
