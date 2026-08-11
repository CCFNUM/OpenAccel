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
#ifdef HAS_UTOPIA
#include "ContextUtopia.h"
#endif /* HAS_UTOPIA */

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

    void applyOverrides() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

#ifdef HAS_UTOPIA
    // This API is only defined for segregated FSI via Utopia
    void utopiaUpdateSolution();
    void utopiaResetFluid();
    void utopiaReportNormIncrement();

protected:
    void utopiaSegregatedVelocityIncrement_();
#else
protected:
#endif /* HAS_UTOPIA */
    void setResidualScales_() override;

#ifdef HAS_UTOPIA
    std::shared_ptr<fsi::ContextUtopia<linearSystem::BLOCKSIZE>> ctx_fsi_;
#endif /* HAS_UTOPIA */

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // NAVIERSTOKESEQUATION_H
