// File : wallScaleDiffusionAssembler.h
// Created : Thu Feb 22 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Assembles node and element terms for wall scale diffusion
// equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef WALLSCALEDIFFUSIONASSEMBLER_H
#define WALLSCALEDIFFUSIONASSEMBLER_H

#include "phiAssembler.h"
#include "wallScaleDiffusionModel.h"

namespace accel
{

class wallScaleDiffusionAssembler : public phiAssembler<1>
{
private:
    wallScaleDiffusionModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    wallScaleDiffusionAssembler(wallScaleDiffusionModel* model)
        : Base(model), model_(model)
    {
    }

protected:
    // kernel drivers
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;
    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;
    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;
    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;

    // Boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
};

} /* namespace accel */

#endif // WALLSCALEDIFFUSIONASSEMBLER_H
