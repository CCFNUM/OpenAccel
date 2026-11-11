// File       : displacementDiffusionEquation.h
// Created    : Tue Nov 26 2024
// Author     : Mhamad Mahdi Alloush
// Description: Displacement diffusion equation for mesh deformation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

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
#include <set>

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

    // domains whose stiffness is evaluated on the initial mesh (computed once)
    std::set<label> stiffnessFrozen_;

    // reads the fixed stiffness from the restart file
    bool restoreStiffness_(STKScalarField& Gamma);
};

} /* namespace accel */

#endif // DISPLACEMENTDIFFUSIONEQUATION_H
