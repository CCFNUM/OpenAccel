// File       : wallScaleDiffusionModel.h
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Model providing wall scale field access for diffusion-based wall
// distance
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WALLSCALEDIFFUSIONMODEL_H
#define WALLSCALEDIFFUSIONMODEL_H

// code
#include "model.h"

namespace accel
{

class wallScaleDiffusionModel : public model
{
public:
    wallScaleDiffusionModel(realm* realm);

    using fieldBroker::yScaleRef;
};

} /* namespace accel */

#endif // WALLSCALEDIFFUSIONMODEL_H
