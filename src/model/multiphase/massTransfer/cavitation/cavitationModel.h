// File       : cavitationModel.h
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Validation of the cavitation parameters of a fluid pair
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef CAVITATIONMODEL_H
#define CAVITATIONMODEL_H

#include "domain.h"

namespace accel
{

// errorMsg with `path` (YAML path) prefixed on an invalid parameter
void validateCavitationModel(const fluidPairModel::massTransfer& cfg,
                             const std::string& path);

} /* namespace accel */

#endif // CAVITATIONMODEL_H
