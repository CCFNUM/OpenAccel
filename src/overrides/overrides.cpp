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
#include "immersedDistance.h"

#include <mpi.h>
#include <iomanip>
#include <iostream>

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
        declareImmersedDistance(*realmPtr_);
    }
}

void overrides::initialize()
{
    if (maskingPtr_)
    {

        // Time our immersed signed-distance calculation.
        const auto comm = messager::comm();

        MPI_Barrier(comm);
        const double start = MPI_Wtime();

        maskingPtr_->initialize();

        computeImmersedDistance(*realmPtr_, *maskingPtr_);

        const double localSeconds = MPI_Wtime() - start;

        double elapsedSeconds = 0.0;
        MPI_Reduce(
            &localSeconds,
            &elapsedSeconds,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0,
            comm);

        if (messager::master())
        {
            std::cout << "SSDF_TOTAL_SECONDS="
                      << std::setprecision(12)
                      << elapsedSeconds
                      << std::endl;
        }
    }
}

YAML::Node overrides::getFSINode() const
{
    assert(fsi_);

    return realmPtr_->simulationRef()
        .getYAMLPhysicalAnalysisNode()["overrides"]["fsi"];
}

} // namespace accel
