// File       : turbulentIntermittencyCorrelationTransitionSSTAssembler.h
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: Assembler for gamma equation in correlation-based
// transition SST model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTASSEMBLER_H
#define TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTASSEMBLER_H

// code
#include "correlationTransitionShearStressTransportModel.h"
#include "phiAssembler.h"

namespace accel
{

class turbulentIntermittencyCorrelationTransitionSSTAssembler
    : public phiAssembler<1>
{
private:
    correlationTransitionShearStressTransportModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    turbulentIntermittencyCorrelationTransitionSSTAssembler(
        correlationTransitionShearStressTransportModel* model)
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

#endif // TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTASSEMBLER_H
