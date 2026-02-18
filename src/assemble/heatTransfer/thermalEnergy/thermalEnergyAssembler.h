// File : thermalEnergyAssembler.h
// Created : Mon Apr 14 2025 10:36:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Assembler for the thermal energy transport equation
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALENERGYASSEMBLER_H
#define THERMALENERGYASSEMBLER_H

#ifndef WITH_THERMAL_TEMPERATURE

#include "heatTransferModel.h"
#include "phiAssembler.h"

namespace accel
{

class thermalEnergyAssembler : public phiAssembler<1>
{
private:
    heatTransferModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    thermalEnergyAssembler(heatTransferModel* model)
        : Base(model), model_(model)
    {
    }

protected:
    void postAssemble_(const domain* domain, Context* ctx) override
    {
        Base::postAssemble_(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    }

    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;

    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;

    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
};

} /* namespace accel */

#endif /* WITH_THERMAL_TEMPERATURE */

#endif // THERMALENERGYASSEMBLER_H
