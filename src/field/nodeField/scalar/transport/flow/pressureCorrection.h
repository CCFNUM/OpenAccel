// File       : pressureCorrection.h
// Created    : Mon May 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar working field for the pressure correction (p')
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef PRESSURECORRECTION_H
#define PRESSURECORRECTION_H

#include "nodeScalarField.h"

namespace accel
{

// Dedicated storage for the pressure correction p' produced by the segregated
// (SIMPLE/SIMPLEC) pressure-correction equation. It exists solely so that the
// segregated velocity correction can use a properly relaxed/limited pressure
// gradient while reconstructing the true increment p' from its own gradient
// field. The field carries a gradient (computeGradient = true) but no boundary
// conditions, no restart state and no time states: it is overwritten every
// outer iteration from the raw linear-solve correction.
class pressureCorrection : public nodeScalarField
{
public:
    // Constructors

    pressureCorrection(realm* realmPtr,
                       const std::string name,
                       unsigned numberOfStates,
                       bool highResolution);
};

} // namespace accel

#endif // PRESSURECORRECTION_H
