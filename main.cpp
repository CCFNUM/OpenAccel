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

#ifdef _WIN32
// Included last so its macros do not pollute the STK/Trilinos headers above.
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

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

#ifdef _WIN32
    // The MSYS2 netcdf DLL bundles the AWS C++ SDK (S3/NCZarr support), whose
    // process-exit handler deadlocks joining its event-loop threads during the
    // Windows shutdown sequence (LdrShutdownProcess -> netcdf onexit table ->
    // aws-c-common SleepConditionVariableSRW). All simulation output has
    // already been written and closed by this point, and Kokkos/MPI have been
    // finalized above, so terminate immediately to skip the deadlocking
    // per-DLL teardown. TerminateProcess is the only exit path that bypasses
    // the onexit tables entirely. Re-raising a caught signal below would go
    // through that same teardown (the CRT's default SIGINT/SIGTERM action is
    // _exit(3)), so report an interrupted run with that status here instead.
    std::cout.flush();
    std::cerr.flush();
    ::TerminateProcess(::GetCurrentProcess(), g_signalSent != 0 ? 3 : 0);
#endif

    if (g_signalSent != 0)
    {
        std::signal(g_signalID, SIG_DFL);
        std::raise(g_signalID);
    }

    return 0;
}
