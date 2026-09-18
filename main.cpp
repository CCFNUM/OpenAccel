// File       : main.cpp
// Created    : Thu Mar 26 2026 10:17:26 (+0100)
// Author     : Fabian Wermelinger
// Description: Main application entry
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include <Kokkos_Core.hpp>
#include <csignal>
#include <mpi.h>
#ifdef HAS_PETSC
#include <petscsys.h>
#if PETSC_VERSION_LT(3, 18, 1)
#define ErrorWrapPetscCall(c) CHKERRQ(c)
#else
#define ErrorWrapPetscCall(c) PetscCall(c)
#endif
#endif /* HAS_PETSC */
#ifdef HAS_HYPRE
#include <HYPRE_utilities.h>
#endif /* HAS_HYPRE */

// code libraries
#include "macros.h"
#include "simulation.h"

volatile sig_atomic_t g_signalSent = 0;
volatile sig_atomic_t g_signalID = 0;

extern "C" void handleSignal(int signal)
{
    g_signalSent = 1;
    g_signalID = signal;
    std::signal(signal, SIG_DFL); // next one kills immediately (no cleanup)
}

int main(int argc, char* argv[])
{
    // The package tackles real world physics in 2D or 3D
    assert(SPATIAL_DIM >= 2);

    using Sim = ::accel::simulation;

    // register handlers
    std::signal(SIGINT, handleSignal);
    std::signal(SIGTERM, handleSignal);

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
        accel::errorMsg("Provided MPI thread-level support is not sufficient");
    }
#ifdef HAS_PETSC
    // Initialize the Petsc environment
    ErrorWrapPetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
#endif /* HAS_PETSC */

#ifdef HAS_HYPRE
    HYPRE_Initialize();
#endif /* HAS_HYPRE */
    Kokkos::initialize(argc, argv);
    {
        // create and run simulation
        Sim* realm = new Sim(argc, const_cast<const char**>(argv));

        realm->run();

        delete realm;
    }
    Kokkos::finalize();

#ifdef HAS_HYPRE
    HYPRE_Finalize();
#endif /* HAS_HYPRE */

#ifdef HAS_PETSC
    // Finalize the Petsc environment.
    ErrorWrapPetscCall(PetscFinalize());
#endif /* HAS_PETSC */
    MPI_Finalize();

    if (g_signalSent != 0)
    {
        std::signal(g_signalID, SIG_DFL);
        std::raise(g_signalID);
    }

    return 0;
}
