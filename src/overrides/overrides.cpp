// File       : overrides.cpp
// Created    : Mon Feb 09 2026 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "overrides.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

overrides::overrides(realm* realm) : realmPtr_(realm)
{
}

void overrides::read(const YAML::Node& inputNode)
{
#ifdef HAS_UTOPIA
    const auto& physicalAnalysisNode = inputNode["physical_analysis"];

    if (physicalAnalysisNode["overrides"])
    {
        const auto& overridesNode = physicalAnalysisNode["overrides"];

        if (overridesNode["fsi"])
        {
            fsi_ = true;
        }
    }
#endif /* HAS_UTOPIA */
}

#ifdef HAS_UTOPIA
YAML::Node overrides::getFSINode() const
{
    assert(fsi_);

    return realmPtr_->simulationRef()
        .getYAMLPhysicalAnalysisNode()["overrides"]["fsi"];
}
#endif /* HAS_UTOPIA */

} // namespace accel
