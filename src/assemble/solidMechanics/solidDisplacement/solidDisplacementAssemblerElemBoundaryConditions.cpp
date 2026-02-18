// File : solidDisplacementAssemblerElemBoundaryConditions.cpp
// Created : Sun Feb 01 2026 02:30:10 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidDisplacementAssembler.h"

namespace accel
{

void solidDisplacementAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    phiAssembler<SPATIAL_DIM>::assembleElemTermsBoundary_(domain, ctx);
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} /* namespace accel */
