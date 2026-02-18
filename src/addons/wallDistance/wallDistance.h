// File : wallDistance.h
// Created : Tue Feb 10 2026 12:50:04 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Computes and manages wall distance fields using wall scale
// methods
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WALLDISTANCE_H
#define WALLDISTANCE_H

// code
#include "equation.h"
#include "fieldBroker.h"
#include "wallScale.h"

namespace accel
{

class wallDistance : public fieldBroker
{
private:
    std::unique_ptr<wallScale> wallScalePtr_ = nullptr;

public:
    // Constructors
    wallDistance(realm* realm);

    using fieldBroker::yMinRef;
    using fieldBroker::yScaleRef;

    // Methods

    void reset();

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // WALLDISTANCE_H
