// File       : volumeFractionAssembler.h
// Created    : Mon Jan 27 2025
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the volume fraction transport in free-surface
// flows
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

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
    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;

#ifdef HAS_INTERFACE
    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;
#endif /* HAS_INTERFACE */

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
