// File       : thermoModel.cpp
// Created    : Thu Apr 03 2025 17:05:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "thermoModel.h"

namespace accel
{

thermoModel::thermoModel(realm* realm) : model(realm)
{
}

const scalar thermoModel::universalGasConstant_ = 8314.4598;

} /* namespace accel */
