// File       : thermoModel.cpp
// Created    : Thu Apr 03 2025 17:05:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "thermoModel.h"

namespace accel
{

thermoModel::thermoModel(realm* realm) : model(realm)
{
}

const scalar thermoModel::universalGasConstant_ = 8314.4598;

} /* namespace accel */
