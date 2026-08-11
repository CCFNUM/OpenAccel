// File       : overrides.cpp
// Created    : Mon Feb 09 2026 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "overrides.h"
#include "controls.h"
#include "macros.h"
#include "messager.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

overrides::overrides(realm* realm) : realmPtr_(realm)
{
}

void overrides::read(const YAML::Node& inputNode)
{
    const auto& physicalAnalysisNode = inputNode["physical_analysis"];

    if (!physicalAnalysisNode["overrides"])
    {
        return;
    }

    const auto& overridesNode = physicalAnalysisNode["overrides"];

    if (overridesNode["fsi"])
    {
        fsi_ = true;
    }

    if (overridesNode["masking"])
    {
        maskingPtr_ = std::make_unique<masking>(realmPtr_);
        maskingPtr_->read(overridesNode["masking"]);

        // the fields have to exist before the mesh is populated
        maskingPtr_->setup();
    }
}

void overrides::initialize()
{
    if (maskingPtr_)
    {
        maskingPtr_->initialize();
    }
}

YAML::Node overrides::getFSINode() const
{
    assert(fsi_);

    return realmPtr_->simulationRef()
        .getYAMLPhysicalAnalysisNode()["overrides"]["fsi"];
}

} // namespace accel
