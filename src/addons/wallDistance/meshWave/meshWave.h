// File       : meshWave.h
// Created    : Sat Mar 01 2026
// Author     : Mhamad Mahdi Alloush
// Description: Wall distance computation via direct Euclidean distance to
// nearest wall node (meshWave method)
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MESHWAVE_H
#define MESHWAVE_H

// code
#include "types.h"

namespace accel
{

class wallDistance;

class meshWave
{
private:
    wallDistance* wallDistancePtr_;

    void computeWallDistance_();

public:
    // Constructors
    meshWave(wallDistance* wallDistancePtr);

    // Methods

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // MESHWAVE_H
