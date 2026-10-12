// File       : navierStokesEquation.h
// Created    : Wed Mar 13 2024 10:02:56 (+0100)
// Author     : Fabian Wermelinger
// Description: Coupled momentum equations
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef NAVIERSTOKESEQUATION_H
#define NAVIERSTOKESEQUATION_H

// code
#include "equation.h"
#include "flowModel.h"
#include "linearSystem.h"
#include "navierStokesAssembler.h"

namespace accel
{

class navierStokesEquation : public equation, public linearSystem<SPATIAL_DIM>
{
    flowModel* model_;

    using Assembler = navierStokesAssembler;

public:
    static constexpr equationID ID = equationID::coupledNavierStokes;

    navierStokesEquation(realm* realm, flowModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

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

#endif // NAVIERSTOKESEQUATION_H
