// File       : main.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <Kokkos_Core.hpp>
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

    return 0;
}
