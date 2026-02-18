// File : thermalTemperatureAssembler.h
// Created : Thu Apr 14 2024 8:36:38 (+0100)
// Author : Fabian Wermelinger
// Description: Assembler for the temperature-based heat transfer equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALTEMPERATUREASSEMBLER_H
#define THERMALTEMPERATUREASSEMBLER_H

#ifdef WITH_THERMAL_TEMPERATURE

#include "heatTransferModel.h"
#include "phiAssembler.h"

namespace accel
{

class thermalTemperatureAssembler : public phiAssembler<1>
{
private:
    heatTransferModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    thermalTemperatureAssembler(heatTransferModel* model)
        : Base(model), model_(model)
    {
    }

protected:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;

    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;

    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;
    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;

#ifdef HAS_INTERFACE
    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;
    void assembleElemTermsInterfaceSideHTC_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);
#endif /* HAS_INTERFACE */

    // boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    void assembleElemTermsBoundaryWallFixedValue_(const domain* domain,
                                                  const boundary* boundary,
                                                  Context* ctx) override;
    void assembleElemTermsBoundaryWallSpecifiedFlux_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx) override;
    void assembleElemTermsBoundaryWallMixed_(const domain* domain,
                                             const boundary* boundary,
                                             Context* ctx) override;
    void assembleElemTermsBoundaryInletFixedValue_(const domain* domain,
                                                   const boundary* boundary,
                                                   Context* ctx) override;
    void assembleElemTermsBoundaryOutletZeroGradient_(const domain* domain,
                                                      const boundary* boundary,
                                                      Context* ctx) override;
    void assembleElemTermsBoundaryOpening_(const domain* domain,
                                           const boundary* boundary,
                                           Context* ctx) override;

    void postAssemble_(const domain* domain, Context* ctx) override
    {
        Base::postAssemble_(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    }
};

} /* namespace accel */

#endif /* WITH_THERMAL_TEMPERATURE */

#endif // THERMALTEMPERATUREASSEMBLER_H
