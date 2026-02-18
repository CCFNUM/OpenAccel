// File : massFlowRate.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "massFlowRate.h"
#include "boundary.h"
#include "domain.h"
#include "mesh.h"
#include "realm.h" // required for initialization
#include "simulation.h"
#include "zone.h"

namespace accel
{

massFlowRate::massFlowRate(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates)
    : elementScalarField(realmPtr, name, numberOfStates)
{
    // create divergence field
    std::string divergenceFieldName = "divergence";

    // Find the position of the dot
    size_t dot_pos = name.find('.');

    // Append the dot and the remaining part of name1 (if found)
    if (dot_pos != std::string::npos)
    {
        divergenceFieldName += name.substr(dot_pos);
    }

    // Instantiate divergence field
    divFieldPtr_ = std::make_unique<nodeField<1>>(
        this->meshPtr(), divergenceFieldName, 1, false);

    // set size of massFlowRateFraction
    sideMassFlowRateFraction_.resize(this->meshPtr()->nZones());
    for (label iZone = 0; iZone < this->meshPtr()->nZones(); iZone++)
    {
        sideMassFlowRateFraction_[iZone].resize(
            this->meshPtr()->zonePtr(iZone)->nBoundaries(), 0.0);
    }
}

// Access

nodeField<1>& massFlowRate::divRef()
{
    return *divFieldPtr_;
}

const nodeField<1>& massFlowRate::divRef() const
{
    return *divFieldPtr_;
}

} // namespace accel
