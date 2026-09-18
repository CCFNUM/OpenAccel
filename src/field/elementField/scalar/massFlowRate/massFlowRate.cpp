// File       : massFlowRate.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

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

    // The mass flux and its divergence carry iteration history through the
    // Rhie-Chow interpolation; persist them so a restart resumes with the
    // correct (consistently assembled) state instead of a raw reconstruction.
    realmPtr->registerRestartField(name);
    realmPtr->registerRestartField(divergenceFieldName);

    // The boundary side mass flux carries the Rhie-Chow history as well.
    // Registered here unconditionally (not at side field instantiation) so the
    // restart set does not depend on registration order: the interface path
    // instantiates the side field without restart registration, which would
    // otherwise shadow the registration in elementField::registerSideField.
    // Interface-side values are not yet restored on a restart, only boundary
    // ones (elementField::restoreSideField); they will be eventually.
    realmPtr->registerRestartField(name + "_side");

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
