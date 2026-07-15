// File       : totalEnergyAssembler.h
// Created    : Thu Apr 02 2025 15:40:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the total energy conservation equation
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

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

    void preAssemble(const domain* domain, Context* ctx) override
    {
        Base::preAssemble(domain, ctx);

        // Zero qDot on interface IP slots before scatter.
        auto& qDotSideSTKFieldRef =
            model_->qDotRef().sideFieldRef().stkFieldRef();
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isConformalTreatment())
                continue;
            const bool zeroBothSides = interf->isInternal();
            const bool zeroOwnSide = !interf->isInternal();

            if (!zeroBothSides && !zeroOwnSide)
                continue; // nothing to zero from this side

            auto zeroSide = [&](const interfaceSideInfo* sidePtr)
            {
                for (const auto& faceVec : sidePtr->ipInfoVec())
                {
                    for (const ipInfo* ip : faceVec)
                    {
                        if (ip == nullptr || ip->isExposed_)
                            continue;
                        scalar* q = stk::mesh::field_data(qDotSideSTKFieldRef,
                                                          ip->currentFace_);
                        q[ip->currentGaussPointId_] = 0.0;
                    }
                }
            };
            if (zeroBothSides)
            {
                zeroSide(interf->masterInfoPtr());
                zeroSide(interf->slaveInfoPtr());
            }
            else
            {
                zeroSide(interf->interfaceSideInfoPtr(domain->index()));
            }
        }
    }

    void postAssemble(const domain* domain, Context* ctx) override
    {
        Base::postAssemble(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx, 0.75);
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

    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;
    void assembleElemTermsInterfaceSideHTC_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);

    // boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
    void assembleElemTermsBoundarySymmetry_(const domain* domain,
                                            const boundary* boundary,
                                            Context* ctx) override;
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
};

} /* namespace accel */

#endif // TOTALENERGYASSEMBLER_H
