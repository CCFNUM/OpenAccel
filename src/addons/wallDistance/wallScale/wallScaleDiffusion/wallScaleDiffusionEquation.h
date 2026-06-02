// File       : wallScaleDiffusionEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Equation class for solving the wall scale diffusion problem
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef WALLSCALEDIFFUSIONEQUATION_H
#define WALLSCALEDIFFUSIONEQUATION_H

// code
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "types.h"
#include "wallScaleDiffusionAssembler.h"
#include "wallScaleDiffusionModel.h"

namespace accel
{

class wallScaleDiffusionEquation : public equation,
                                   public wallScaleDiffusionModel,
                                   public linearSystem<1>
{
private:
    using Assembler = wallScaleDiffusionAssembler;

public:
    static constexpr equationID ID = equationID::wallScaleDiffusion;

    wallScaleDiffusionEquation(realm* realm);

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

    stk::mesh::PartVector collectDirichletBoundaryParts_();

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // WALLSCALEDIFFUSIONEQUATION_H
