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

#ifdef HAS_UTOPIA
    if (overridesNode["fsi"])
    {
        fsi_ = true;
    }
#endif /* HAS_UTOPIA */

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

        // Instrument existing operations only; no additional geometry calculations.
        const auto comm = messager::comm();
        MaskingTimings maskingTimings;
        ImmersedDistanceTimings distanceTimings;
        MPI_Barrier(comm);
        const double start = MPI_Wtime();

        maskingPtr_->initialize(&maskingTimings);
        const double maskingSeconds = MPI_Wtime() - start;

        const double distanceStart = MPI_Wtime();
        computeImmersedDistance(*realmPtr_, *maskingPtr_, &distanceTimings);
        const double immersedSeconds = MPI_Wtime() - distanceStart;
        const double totalSeconds = MPI_Wtime() - start;

        // Each value is the maximum rank time. Components add approximately to the
        // total in serial; maxima may come from different ranks under MPI.
        const char* names[] = {
            "SSDF_TOTAL_SECONDS", "SSDF_MASKING_INIT_SECONDS",
            "SSDF_SURFACE_SETUP_SECONDS", "SSDF_GRAPH_SECONDS",
            "SSDF_COVERED_CLASSIFICATION_SECONDS", "SSDF_LAYER_PROPAGATION_SECONDS",
            "SSDF_REFERENCE_PAIRING_SECONDS", "SSDF_MASKING_OTHER_SECONDS",
            "SSDF_DISTANCE_SECONDS", "SSDF_SIGN_SECONDS",
            "SSDF_DISTANCE_SIGN_2D_SECONDS", "SSDF_FIELD_SECONDS",
            "SSDF_PREPARATION_OTHER_SECONDS"};
        const double local[] = {
            totalSeconds, maskingSeconds, maskingTimings.surfaceSetup,
            maskingTimings.graph, maskingTimings.coveredClassification,
            maskingTimings.layers, maskingTimings.referencePairing,
            maskingSeconds - maskingTimings.surfaceSetup - maskingTimings.graph
                - maskingTimings.coveredClassification - maskingTimings.layers
                - maskingTimings.referencePairing,
            distanceTimings.distance, distanceTimings.sign,
            distanceTimings.combined2D, distanceTimings.field,
            immersedSeconds - distanceTimings.distance - distanceTimings.sign
                - distanceTimings.combined2D - distanceTimings.field};
        constexpr int count = sizeof(local) / sizeof(local[0]);
        double maximum[count] = {};
        MPI_Reduce(local, maximum, count, MPI_DOUBLE, MPI_MAX, 0, comm);
        if (messager::master())
        {
            for (int i = 0; i < count; ++i)
                std::cout << names[i] << "=" << std::setprecision(12)
                          << maximum[i] << std::endl;
        }
    }
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
