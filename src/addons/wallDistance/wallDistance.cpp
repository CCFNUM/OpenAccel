// File       : wallDistance.cpp
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "wallDistance.h"

namespace accel
{

wallDistance::wallDistance(realm* realm) : fieldBroker(realm)
{
    yMinRef();

    // create deformation
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.geometricWallDistanceCalculation_)
    {
        errorMsg("Geometric wall distance calculator not implemented yet");
    }
    else
    {
        wallScalePtr_ = std::make_unique<wallScale>(this);
    }
}

void wallDistance::setup()
{
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.geometricWallDistanceCalculation_)
    {
        errorMsg("Geometric wall distance calculator not implemented yet");
    }
    else
    {
        assert(wallScalePtr_);
        wallScalePtr_->setup();
    }
}

void wallDistance::reset()
{
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.geometricWallDistanceCalculation_)
    {
        errorMsg("Geometric wall distance calculator not implemented yet");
    }
    else
    {
        assert(wallScalePtr_);
    }
}

void wallDistance::initialize()
{
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.geometricWallDistanceCalculation_)
    {
        errorMsg("Geometric wall distance calculator not implemented yet");
    }
    else
    {
        assert(wallScalePtr_);
        wallScalePtr_->initialize();
    }
}

void wallDistance::update()
{
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.geometricWallDistanceCalculation_)
    {
        errorMsg("Geometric wall distance calculator not implemented yet");
    }
    else
    {
        assert(wallScalePtr_);
        wallScalePtr_->update();
    }
}

} // namespace accel
