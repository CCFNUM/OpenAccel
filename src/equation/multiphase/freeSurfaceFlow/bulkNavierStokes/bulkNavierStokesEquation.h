// File       : bulkNavierStokesEquation.h
// Created    : Thu Jan 30 2025 10:02:56 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Bulk Navier-Stokes momentum equation for free-surface flow
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef BULKNAVIERSTOKESEQUATION_H
#define BULKNAVIERSTOKESEQUATION_H

// code
#include "freeSurfaceFlowModel.h"
#include "navierStokesEquation.h"

namespace accel
{

class bulkNavierStokesEquation : public navierStokesEquation
{
    freeSurfaceFlowModel* model_;

    using Assembler = navierStokesAssembler;

public:
    bulkNavierStokesEquation(realm* realm, freeSurfaceFlowModel* model);

    void setup() override;

    void postInitialize() override;

    void preSolve() override;

    void preTimeStep() override;
};

} /* namespace accel */

#endif // BULKNAVIERSTOKESEQUATION_H
