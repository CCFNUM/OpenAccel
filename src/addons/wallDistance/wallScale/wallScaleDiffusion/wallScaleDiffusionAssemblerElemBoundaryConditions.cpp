// File : wallScaleDiffusionAssemblerElemBoundaryConditions.cpp
// Created : Thu Jun 19 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "wallScaleDiffusionAssembler.h"

namespace accel
{

void wallScaleDiffusionAssembler::assembleElemTermsBoundary_(
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
            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            assembleElemTermsBoundaryWallFixedValue_(
                                domain, boundary, ctx);
                            break;

                        case boundaryConditionType::zeroGradient:
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            default:
                break;
        }
    }
}

} /* namespace accel */
