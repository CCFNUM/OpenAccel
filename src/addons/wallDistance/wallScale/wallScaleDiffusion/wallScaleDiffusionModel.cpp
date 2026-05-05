// File       : wallScaleDiffusionModel.cpp
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "wallScaleDiffusionModel.h"
#include "simulation.h"

namespace accel
{

wallScaleDiffusionModel::wallScaleDiffusionModel(realm* realm) : model(realm)
{
    // create field instances
    yScaleRef();
}

} // namespace accel
