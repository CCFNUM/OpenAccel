// File       : phasicContinuityAssembler.h
// Description: Assembler for conservative Euler-Euler phasic continuity

#ifndef PHASICCONTINUITYASSEMBLER_H
#define PHASICCONTINUITYASSEMBLER_H

#include "eulerEulerModel.h"
#include "phiAssembler.h"

namespace accel
{

class phasicContinuityAssembler : public phiAssembler<1>
{
public:
    phasicContinuityAssembler(eulerEulerModel* model, label phaseIndex);

protected:
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
        return model_->intrinsicMDotRef(phaseIndex_);
    }

    const elementField<scalar, 1>& mDotRef() const override
    {
        return model_->intrinsicMDotRef(phaseIndex_);
    }

    nodeField<1>& divRef() override
    {
        return model_->intrinsicMDotRef(phaseIndex_).divRef();
    }

    const nodeField<1>& divRef() const override
    {
        return model_->intrinsicMDotRef(phaseIndex_).divRef();
    }

    velocity& advectionVelocityRef() override
    {
        return model_->URef(phaseIndex_);
    }

    const velocity& advectionVelocityRef() const override
    {
        return model_->URef(phaseIndex_);
    }

private:
    eulerEulerModel* model_;
    label phaseIndex_;
};

} // namespace accel

#endif // PHASICCONTINUITYASSEMBLER_H
