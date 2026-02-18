// File : totalEnergyAssembler.h
// Created : Thu Apr 02 2025 15:40:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Assembler for the total energy conservation equation
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TOTALENERGYASSEMBLER_H
#define TOTALENERGYASSEMBLER_H

#include "heatTransferModel.h"
#include "phiAssembler.h"

namespace accel
{

class totalEnergyAssembler : public phiAssembler<1>
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
    totalEnergyAssembler(heatTransferModel* model) : Base(model), model_(model)
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

    void assembleElemTermsInterfaceSideHybrid_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);
#endif /* HAS_INTERFACE */

    // boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
    void assembleElemTermsBoundarySymmetry_(const domain* domain,
                                            const boundary* boundary,
                                            Context* ctx);
    void assembleElemTermsBoundaryWallZeroGradient_(const domain* domain,
                                                    const boundary* boundary,
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

#endif // TOTALENERGYASSEMBLER_H
