// File       : wallScaleDiffusionModel.h
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Model providing wall scale field access for diffusion-based wall
//              distance
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

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
