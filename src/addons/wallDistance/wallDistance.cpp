// File       : wallDistance.cpp
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "wallDistance.h"

namespace accel
{

wallDistance::wallDistance(realm* realm) : fieldBroker(realm)
{
    yMinRef();

    const auto method =
        this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.wallDistanceMethod_;

    switch (method)
    {
        case wallDistanceMethod::poisson:
            wallScalePtr_ = std::make_unique<wallScale>(this);
            break;
        case wallDistanceMethod::meshWave:
            meshWavePtr_ = std::make_unique<meshWave>(this);
            break;
        case wallDistanceMethod::signedDistanceFunction:
            errorMsg("Signed distance function wall distance method not "
                     "implemented yet");
            break;
    }
}

void wallDistance::setup()
{
    const auto method =
        this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.wallDistanceMethod_;

    switch (method)
    {
        case wallDistanceMethod::poisson:
            assert(wallScalePtr_);
            wallScalePtr_->setup();
            break;
        case wallDistanceMethod::meshWave:
            assert(meshWavePtr_);
            meshWavePtr_->setup();
            break;
        case wallDistanceMethod::signedDistanceFunction:
            errorMsg("Signed distance function wall distance method not "
                     "implemented yet");
            break;
    }
}

void wallDistance::reset()
{
    const auto method =
        this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.wallDistanceMethod_;

    switch (method)
    {
        case wallDistanceMethod::poisson:
            assert(wallScalePtr_);
            break;
        case wallDistanceMethod::meshWave:
            assert(meshWavePtr_);
            break;
        case wallDistanceMethod::signedDistanceFunction:
            errorMsg("Signed distance function wall distance method not "
                     "implemented yet");
            break;
    }
}

void wallDistance::initialize()
{
    const auto method =
        this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.wallDistanceMethod_;

    switch (method)
    {
        case wallDistanceMethod::poisson:
            assert(wallScalePtr_);
            wallScalePtr_->initialize();
            break;
        case wallDistanceMethod::meshWave:
            assert(meshWavePtr_);
            meshWavePtr_->initialize();
            break;
        case wallDistanceMethod::signedDistanceFunction:
            errorMsg("Signed distance function wall distance method not "
                     "implemented yet");
            break;
    }
}

void wallDistance::update()
{
    const auto method =
        this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.wallDistanceMethod_;

    switch (method)
    {
        case wallDistanceMethod::poisson:
            assert(wallScalePtr_);
            wallScalePtr_->update();
            break;
        case wallDistanceMethod::meshWave:
            assert(meshWavePtr_);
            meshWavePtr_->update();
            break;
        case wallDistanceMethod::signedDistanceFunction:
            errorMsg("Signed distance function wall distance method not "
                     "implemented yet");
            break;
    }
}

} // namespace accel
