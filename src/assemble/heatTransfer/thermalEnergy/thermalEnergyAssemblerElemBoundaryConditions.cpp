// File       : thermalEnergyAssemblerElemBoundaryConditions.cpp
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WITH_THERMAL_TEMPERATURE

#include "thermalEnergyAssembler.h"

namespace accel
{

void thermalEnergyAssembler::assembleElemTermsBoundary_(const domain* domain,
                                                        Context* ctx)
{
    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            model_->TRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                assembleElemTermsBoundaryWallFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::specifiedFlux:
                            {
                                assembleElemTermsBoundaryWallSpecifiedFlux_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            break;

                        case boundaryConditionType::mixed:
                            {
                                assembleElemTermsBoundaryWallMixed_(
                                    domain, boundary, ctx);
                            }
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
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                assembleElemTermsBoundaryInletFixedValue_(
                                    domain, boundary, ctx);
                            }
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

            default:
                break;
        }
    }
}

} // namespace accel

#endif
