// File : displacementDiffusionEquation.h
// Created : Tue Nov 26 2024
// Author : Mhamad Mahdi Alloush
// Description: Displacement diffusion equation for mesh deformation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DISPLACEMENTDIFFUSIONEQUATION_H
#define DISPLACEMENTDIFFUSIONEQUATION_H

// code
#include "displacementDiffusionAssembler.h"
#include "displacementDiffusionModel.h"
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "types.h"

namespace accel
{

class displacementDiffusionEquation : public equation,
                                      public displacementDiffusionModel,
                                      public linearSystem<SPATIAL_DIM>
{
    using Assembler = displacementDiffusionAssembler;

public:
    static constexpr equationID ID = equationID::displacementDiffusion;

    displacementDiffusionEquation(realm* realm);

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

    stk::mesh::PartVector collectDirichletBoundaryParts_();

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // DISPLACEMENTDIFFUSIONEQUATION_H
