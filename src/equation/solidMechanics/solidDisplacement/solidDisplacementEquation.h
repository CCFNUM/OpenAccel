// File       : solidDisplacementEquation.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Solid displacement equation with Aitken acceleration for
// structural mechanics
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOLIDDISPLACEMENTEQUATION_H
#define SOLIDDISPLACEMENTEQUATION_H

// code
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "solidDisplacementAssembler.h"
#include "solidMechanicsModel.h"
#include "types.h"

namespace accel
{

class solidDisplacementEquation : public equation,
                                  public solidMechanicsModel,
                                  public linearSystem<SPATIAL_DIM>
{
    using Assembler = solidDisplacementAssembler;

public:
    static constexpr equationID ID = equationID::solidDisplacement;

    solidDisplacementEquation(realm* realm);

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

    void applyDependencyUpdates_(const domain* domain,
                                 const stk::mesh::EntityRank entityRank,
                                 STKScalarField& stk_dst) override;

    stk::mesh::PartVector collectDirichletBoundaryParts_();

private:
    std::unique_ptr<Assembler> assembler_;

    // ALE-FSI
};

} /* namespace accel */

#endif // SOLIDDISPLACEMENTEQUATION_H
