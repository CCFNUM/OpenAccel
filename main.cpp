// File       : main.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <csignal>
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

namespace
{
// Global pointer used by the signal handler to trigger a clean shutdown.
::accel::simulation* g_realm = nullptr;

void handleSignal(int sig)
{
    delete g_realm;
    g_realm = nullptr;
    // Restore the default handler and re-raise so the process exits with the
    // correct status (e.g. SIGINT still shows as interrupted to the shell).
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}
} // namespace

int main(int argc, char* argv[])
{
    // The package tackles real world physics in 2D or 3D
    assert(SPATIAL_DIM >= 2);

    using Sim = ::accel::simulation;

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

        // Register signal handlers so that SIGINT/SIGTERM trigger the
        // simulation destructor (which closes any open gnuplot windows).
        g_realm = realm;
        std::signal(SIGINT, handleSignal);
        std::signal(SIGTERM, handleSignal);

        realm->run();

        g_realm = nullptr;
        std::signal(SIGINT, SIG_DFL);
        std::signal(SIGTERM, SIG_DFL);
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

    return 0;
}
