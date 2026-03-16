// File       : navierStokesAssembler.h
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the Navier-Stokes momentum equations in segregated
// flow
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NAVIERSTOKESASSEMBLER_H
#define NAVIERSTOKESASSEMBLER_H

// code
#include "phiAssembler.h"

namespace accel
{

// Forward declarations
class flowModel;
class domain;
class cnavierStokesAssembler;

// TODO: Some specializations for navierStokesAssembler based on phiAssembler
// result in some degree of code duplication in this class. If possible,
// phiAssembler should be generic enough to handle momentum assembly without the
// need for implementing a new class.
class navierStokesAssembler : public phiAssembler<SPATIAL_DIM>
{
public:
    using Base = phiAssembler<SPATIAL_DIM>;
    using Base::BLOCKSIZE;
    using Context = typename Base::Context;

protected:
    friend cnavierStokesAssembler;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;
    using Base::field_broker_;

public:
    // Constructors
    navierStokesAssembler(flowModel* model);

    void setupDUCoefficients(const domain* domain);

    void computeDUCoefficients(const domain* domain, Context* ctx);

    void assembleNormalRelaxation(const domain* domain,
                                  Context* ctx,
                                  const scalar urf);

protected:
    void preAssemble_(const domain*, Context*) override
    {
    }

    void postAssemble_(const domain*, Context*) override;

    void assembleBoundaryRelaxation_(const domain* domain,
                                     Vector& b,
                                     const scalar urf) override;

    void applySymmetryConditions_(const domain* domain, Context* ctx) override;

private:
    flowModel* model_;

    // kernel drivers
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
    void assembleElemTermsBoundarySlipWall_(const domain* domain,
                                            const boundary* boundary,
                                            Context* ctx);
    void
    assembleElemTermsBoundaryInletSpecifiedVelocity_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx);
    void
    assembleElemTermsBoundaryOutletSpecifiedPressure_(const domain* domain,
                                                      const boundary* boundary,
                                                      Context* ctx);
    void assembleElemTermsBoundaryOutletOutflow_(const domain* domain,
                                                 const boundary* boundary,
                                                 Context* ctx);
    void assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
    void assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
    void assembleElemTermsBoundaryOpening_(const domain* domain,
                                           const boundary* boundary,
                                           Context* ctx) override;
    void
    assembleElemTermsBoundaryInletSpecifiedPressure_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx);
    void assembleElemTermsBoundaryInletSpecifiedTotalPressure_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);
};

} // namespace accel

#endif // NAVIERSTOKESASSEMBLER_H
