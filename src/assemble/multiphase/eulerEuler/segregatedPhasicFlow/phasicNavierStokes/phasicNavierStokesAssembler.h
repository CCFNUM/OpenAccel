// File       : phasicNavierStokesAssembler.h
// Description: Assembler for one Euler-Euler phasic momentum equation

#ifndef PHASICNAVIERSTOKESASSEMBLER_H
#define PHASICNAVIERSTOKESASSEMBLER_H

#include "eulerEulerModel.h"
#include "phiAssembler.h"

namespace accel
{

class phasicNavierStokesAssembler : public phiAssembler<SPATIAL_DIM>
{
public:
    using Base = phiAssembler<SPATIAL_DIM>;
    using Context = typename Base::Context;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;

    phasicNavierStokesAssembler(eulerEulerModel* model, label phaseIndex);

    void setup(const std::vector<std::shared_ptr<domain>>& domains);
    void postAssemble(const domain* domain, Context* ctx) override;

protected:
    nodeField<1, SPATIAL_DIM>& rhoRef() override
    {
        return model_->phaseMassCoefficientRef(phaseIndex_);
    }

    const nodeField<1, SPATIAL_DIM>& rhoRef() const override
    {
        return model_->phaseMassCoefficientRef(phaseIndex_);
    }

    elementField<scalar, 1>& mDotRef() override
    {
        return model_->mDotRef(phaseIndex_);
    }

    const elementField<scalar, 1>& mDotRef() const override
    {
        return model_->mDotRef(phaseIndex_);
    }

    nodeField<1>& divRef() override
    {
        return model_->mDotRef(phaseIndex_).divRef();
    }

    const nodeField<1>& divRef() const override
    {
        return model_->mDotRef(phaseIndex_).divRef();
    }

    velocity& advectionVelocityRef() override
    {
        return model_->URef(phaseIndex_);
    }

    const velocity& advectionVelocityRef() const override
    {
        return model_->URef(phaseIndex_);
    }

    void assembleNodeTermsFused_(const domain* domain,
                                 Context* ctx) override;

    // The phase momentum equation is a vector transport equation: the generic
    // phiAssembler switch cannot dispatch the `noSlip`/`slip` velocity types,
    // so the boundary dispatch and the vector symmetry treatment are
    // specialized here (see
    // phasicNavierStokesAssemblerElemBoundaryConditions.cpp).
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    void applySymmetryConditions_(const domain* domain, Context* ctx) override;

    // The generic phiAssembler opening kernel is a scalar-transport kernel: on
    // reversed flow it advects in `phi_`'s boundary side value and pins the
    // face nodes towards it. For momentum that is both wrong (the opening
    // entrains along the specified flow direction, not towards a stored
    // velocity) and unsafe, because the phasic velocity registers only the
    // side flow-direction fields at an opening -- never the side/node-side
    // value fields the generic kernel dereferences. This specialization
    // mirrors navierStokesAssembler::assembleElemTermsBoundaryOpening_ with
    // the phase mass flux and the alpha_k-weighted effective viscosity.
    void assembleElemTermsBoundaryOpening_(const domain* domain,
                                           const boundary* boundary,
                                           Context* ctx) override;

private:
    void assembleSources_(const domain* domain, Context* ctx);
    void computeDUCoefficients_(const domain* domain, Context* ctx);

    eulerEulerModel* model_;
    label phaseIndex_;
};

} // namespace accel

#endif // PHASICNAVIERSTOKESASSEMBLER_H
