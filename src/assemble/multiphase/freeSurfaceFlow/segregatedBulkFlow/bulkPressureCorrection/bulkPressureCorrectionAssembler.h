// File       : bulkPressureCorrectionAssembler.h
// Created    : Wed Feb 05 2025 01:09:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the bulk pressure correction in free-surface flows
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BULKPRESSURECORRECTIONASSEMBLER_H
#define BULKPRESSURECORRECTIONASSEMBLER_H

// code
#include "phiAssembler.h"

namespace accel
{

// Forward declarations
class freeSurfaceFlowModel;
class domain;

// TODO: Some specializations for momentumAssembler based on phiAssembler result
// in some degree of code duplication in this class. If possible, phiAssembler
// should be generic enough to handle momentum assembly without the need for
// implementing a new class.
class bulkPressureCorrectionAssembler : public phiAssembler<1>
{
public:
    using Base = phiAssembler<1>;
    using Base::BLOCKSIZE;
    using Context = typename Base::Context;

protected:
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;
    using Base::field_broker_;

    label phaseIndex_ = -1;

public:
    // Constructors
    bulkPressureCorrectionAssembler(freeSurfaceFlowModel* model);

    void adjustMatrixForPressureReference(const domain* domain, Context* ctx);

    void setPhaseIndex(label phaseIndex)
    {
        phaseIndex_ = phaseIndex;
    }

protected:
    void preAssemble_(const domain*, Context*) override
    {
    }

    void postAssemble_(const domain*, Context*) override
    {
    }

private:
    freeSurfaceFlowModel* model_;

    // kernel drivers
    void assembleNodeTermsFused_(const domain* domain, Context* ctx) override;
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;
    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;
    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;
    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;

#ifdef HAS_INTERFACE
    void assembleElemTermsInterfaces_(const domain* domain,
                                      Context* ctx) override;
    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;
    void assembleElemTermsInterfaceSideNoSlipWall_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);
#endif /* HAS_INTERFACE */

    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    // boundary conditions
    void assembleElemTermsBoundarySymmetry_(const domain* domain,
                                            const boundary* boundary,
                                            Context* ctx);
    void assembleElemTermsBoundaryWallNoSlip_(const domain* domain,
                                              const boundary* boundary,
                                              Context* ctx);
    void
    assembleElemTermsBoundaryInletSpecifiedVelocity_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx);
    void assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
    void assembleElemTermsBoundaryInletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
    void
    assembleElemTermsBoundaryOutletSpecifiedPressure_(const domain* domain,
                                                      const boundary* boundary,
                                                      Context* ctx);
    void assembleElemTermsBoundaryOutletOutflow_(const domain* domain,
                                                 const boundary* boundary,
                                                 Context* ctx);
    void assembleElemTermsBoundaryOpening_(const domain* domain,
                                           const boundary* boundary,
                                           Context* ctx) override;
    void assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
    void
    assembleElemTermsBoundaryInletSpecifiedPressure_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx);
};

} // namespace accel

#endif // BULKPRESSURECORRECTIONASSEMBLER_H
