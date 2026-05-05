// File       : volumeFractionAssembler.cpp
// Created    : Mon Jan 27 2025
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "volumeFractionAssembler.h"

namespace accel
{

volumeFractionAssembler::volumeFractionAssembler(freeSurfaceFlowModel* model,
                                                 label phaseIndex)
    : phiAssembler<1>(reinterpret_cast<fieldBroker*>(model)), model_(model),
      phaseIndex_(phaseIndex)
{
}

} /* namespace accel */
