// File       : volumeFractionAssembler.h
// Created    : Mon Jan 27 2025
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the volume fraction transport in free-surface
//              flows
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef VOLUMEFRACTIONASSEMBLER_H
#define VOLUMEFRACTIONASSEMBLER_H

#include "freeSurfaceFlowModel.h"
#include "phiAssembler.h"

namespace accel
{

class volumeFractionAssembler : public phiAssembler<1>
{
private:
    freeSurfaceFlowModel* model_;

    label phaseIndex_ = -1;

public:
    using Base = phiAssembler<1>;

    volumeFractionAssembler(freeSurfaceFlowModel* model, label phaseIndex);

protected:
    // Assembly

    // Node terms: same as the base class (divergence correction and the
    // transient/false transient terms) plus the interphase mass transfer
    // (phase change) source, assembled in the same node loop. The source of
    // this phase is, in the advective form rho_k (d alpha_k/dt + u.grad
    // alpha_k) used by the equation,
    //   S_k = s_k mdot_lv - alpha_k rho_k D,  D = mdot_lv (1/rho_v - 1/rho_l)
    // with s_v = +1, s_l = -1 (D is dropped if the continuity source is off).
    // The -alpha rho D part is implicit when D > 0.
    void assembleNodeTermsFused_(const domain* domain, Context* ctx) override;
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

    // Auxiliary field access

    nodeField<1, SPATIAL_DIM>& rhoRef() override
    {
        return model_->rhoRef(phaseIndex_);
    }

    const nodeField<1, SPATIAL_DIM>& rhoRef() const override
    {
        return model_->rhoRef(phaseIndex_);
    }

    elementField<scalar, 1>& mDotRef() override
    {
        return model_->mDotRef(phaseIndex_);
    }

    const elementField<scalar, 1>& mDotRef() const override
    {
        return model_->mDotRef(phaseIndex_);
    }

    virtual nodeField<1>& divRef() override
    {
        return model_->mDotRef(phaseIndex_).divRef();
    }

    const virtual nodeField<1>& divRef() const override
    {
        return model_->mDotRef(phaseIndex_).divRef();
    }
};

} /* namespace accel */

#endif // VOLUMEFRACTIONASSEMBLER_H
