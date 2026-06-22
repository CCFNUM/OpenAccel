// File       : thermalEnergyAssembler.h
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the thermal energy transport equation
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef THERMALENERGYASSEMBLER_H
#define THERMALENERGYASSEMBLER_H

#ifndef WITH_THERMAL_TEMPERATURE

#include "heatTransferModel.h"
#include "phiAssembler.h"

namespace accel
{

class thermalEnergyAssembler : public phiAssembler<1>
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
    thermalEnergyAssembler(heatTransferModel* model)
        : Base(model), model_(model)
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
                continue; // slave domain of a cross-domain GGI iface

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
        assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    }

protected:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;

    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;

    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;

    void assembleElemTermsInterfaceSideHTC_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);

    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    void assembleElemTermsBoundaryWallFixedValue_(const domain* domain,
                                                  const boundary* boundary,
                                                  Context* ctx) override;

    void assembleElemTermsBoundaryWallMixed_(const domain* domain,
                                             const boundary* boundary,
                                             Context* ctx) override;

    void assembleElemTermsBoundaryWallSpecifiedFlux_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx) override;
};

} /* namespace accel */

#endif /* WITH_THERMAL_TEMPERATURE */

#endif // THERMALENERGYASSEMBLER_H
