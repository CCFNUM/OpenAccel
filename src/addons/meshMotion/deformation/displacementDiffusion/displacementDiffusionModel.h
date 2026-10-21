// File       : displacementDiffusionModel.h
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Model of the displacement diffusion equation
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DISPLACEMENTDIFFUSIONMODEL_H
#define DISPLACEMENTDIFFUSIONMODEL_H

// code
#include "meshDisplacementModel.h"

namespace accel
{

// the wall/interface motion and the Dirichlet parts come from the shared model
class displacementDiffusionModel : public meshDisplacementModel
{
public:
    displacementDiffusionModel(realm* realm);
};

} /* namespace accel */

#endif // DISPLACEMENTDIFFUSIONMODEL_H
