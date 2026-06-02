// File       : displacementDiffusionAssemblerElemBoundaryConditions.cpp
// Created    : Tue Nov 26 2024
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "displacementDiffusionAssembler.h"

namespace accel
{

void displacementDiffusionAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
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
            default:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                        case boundaryConditionType::periodicDisplacement:
                        case boundaryConditionType::rigidBodySolution:
                            assembleElemTermsBoundaryWallFixedValue_(
                                domain, boundary, ctx);
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::symmetry:
                break;
        }
    }
}

} /* namespace accel */
