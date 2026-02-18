// File : volumeFraction.h
// Created : Sun Jan 26 2025 22:09:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Volume fraction field for multiphase flows with compressive
// advection
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef VOLUMEFRACTION_H
#define VOLUMEFRACTION_H

#include "nodeScalarField.h"

namespace accel
{

class volumeFraction : public nodeScalarField
{
public:
    // Constructors

    volumeFraction(realm* realmPtr,
                   const std::string name,
                   unsigned numberOfStates,
                   bool highResolution);

    // Override blending factor computation to use one-sided Barth-Jespersen
    // bound, enabling compressive advection (beta > 1)
    void updateBlendingFactorField(label iZone);
};

} // namespace accel

#endif // VOLUMEFRACTION_H
