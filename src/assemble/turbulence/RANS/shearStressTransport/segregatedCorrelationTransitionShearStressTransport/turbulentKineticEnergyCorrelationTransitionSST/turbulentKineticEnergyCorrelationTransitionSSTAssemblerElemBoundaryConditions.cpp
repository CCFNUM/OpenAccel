// File       :
// turbulentKineticEnergyCorrelationTransitionSSTAssemblerElemBoundaryConditions.cpp
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: Boundary conditions for TKE equation in correlation-based
// transition SST model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyCorrelationTransitionSSTAssembler.h"

namespace accel
{

void turbulentKineticEnergyCorrelationTransitionSSTAssembler::
    assembleElemTermsBoundary_(const domain* domain, Context* ctx)
{
    assert(phi_);

    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            phi_->boundaryConditionRef(domain->index(), iBoundary).type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                        case boundaryConditionType::intensityAndLengthScale:
                        case boundaryConditionType::
                            intensityAndEddyViscosityRatio:
                            assembleElemTermsBoundaryInletFixedValue_(
                                domain, boundary, ctx);
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                assembleElemTermsBoundaryOutletZeroGradient_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    assembleElemTermsBoundaryOpening_(domain, boundary, ctx);
                }
                break;
        }
    }
}

} /* namespace accel */
