// File       : navierStokesAssembler.h
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the Navier-Stokes momentum equations in segregated
//              flow
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef NAVIERSTOKESASSEMBLER_H
#define NAVIERSTOKESASSEMBLER_H

// code
#include "phiAssembler.h"

namespace accel
{

// Forward declarations
class flowModel;
class domain;

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
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;
    using Base::field_broker_;

public:
    // Constructors
    navierStokesAssembler(flowModel* model);

    void preAssemble(const domain* domain, Context* ctx) override
    {
    }

    void postAssemble(const domain*, Context*) override;

    void setupDUCoefficients(const domain* domain);

    void computeDUCoefficients(const domain* domain, Context* ctx);

protected:
    void assembleBoundaryRelaxation_(const domain* domain,
                                     Context* ctx,
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

    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    // boundary conditions
    void assembleElemTermsBoundarySymmetry_(const domain* domain,
                                            const boundary* boundary,
                                            Context* ctx) override;
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
