// File       : eulerEulerModel.h
// Description: Inhomogeneous Euler-Euler multiphase model infrastructure

#ifndef EULEREULERMODEL_H
#define EULEREULERMODEL_H

#include "multiphaseModel.h"

namespace accel
{

class eulerEulerModel : public multiphaseModel
{
public:
    explicit eulerEulerModel(realm* realm);

    using fieldBroker::duRef;
    using fieldBroker::duTildeRef;
    using fieldBroker::dragDiagonalRef;
    using fieldBroker::interphaseMomentumSourceRef;
    using fieldBroker::intrinsicMDotRef;
    using fieldBroker::muRef;
    using fieldBroker::phaseMassCoefficientRef;
    using fieldBroker::psiRef;
    using fieldBroker::updateCompressibility;
    using fieldBroker::updateDensity;
    using fieldBroker::updateDynamicViscosity;
    using fieldBroker::updateVelocity;

    std::vector<label>
    activePhaseIndices(const std::shared_ptr<domain> domain) const;

    std::vector<label> activePhaseIndices(const domain* domain) const;

    bool isPhaseActive(const std::shared_ptr<domain> domain,
                       label globalPhaseIndex) const;

    bool isPrimaryPhase(const std::shared_ptr<domain> domain,
                        label globalPhaseIndex) const;

    const fluidPairModel* findFluidPairModel(
        const std::shared_ptr<domain> domain,
        label globalPhaseIndexA,
        label globalPhaseIndexB) const;

    void setupPhaseFields(const std::shared_ptr<domain> domain);

    void initializePhaseFields(const std::shared_ptr<domain> domain);

    // How far redistributeBodyForces may move a phase's body force away from
    // its raw nodal value, as a fraction of that value. Binds only across a
    // free surface, where alpha jumps by orders of magnitude between
    // neighbouring nodes; elsewhere the correction is far smaller than this.
    static constexpr scalar redistributionLimit_ = 1.0;

    void applyVolumeFractionClosure(const std::shared_ptr<domain> domain);

    void updatePhaseMassCoefficient(const std::shared_ptr<domain> domain,
                                    label phaseIndex);

    void updatePhaseMassCoefficientPrevTime(
        const std::shared_ptr<domain> domain,
        label phaseIndex);

    // Build the per-phase body force alpha_k (rho_k - rho_ref) g and put it on
    // the same discrete footing as the compact pressure gradient by reusing
    // flowModel::redistributeBodyForces, staged one phase at a time through
    // the inherited `body_forces` scratch field.
    void updatePhaseBodyForces(const std::shared_ptr<domain> domain);

    void updateInterphaseMomentumSources(
        const std::shared_ptr<domain> domain);

    void updatePhaseMassFlux(const std::shared_ptr<domain> domain,
                             label phaseIndex,
                             bool includeRhieChow = true);

    void updatePhaseFluxDivergence(const std::shared_ptr<domain> domain,
                                   label phaseIndex);

    void fillPhaseEffectiveViscosity(const domain* domain,
                                     label phaseIndex,
                                     STKScalarField& field) const;

    void correctPhaseVelocity(const std::shared_ptr<domain> domain,
                              label phaseIndex);

    // Boundaries where this phase's velocity is prescribed: a specified-value
    // inlet, or a no-slip / specified-value wall.
    //
    // The phasic momentum equation delegates its boundaries to the generic
    // phiAssembler kernels, which impose the value weakly through the boundary
    // flux and leave the boundary node's own row free. Nothing else pins it,
    // so the prescribed value is simply not honoured -- measured on the bubble
    // column at t=0.2: |U.air| on the no-slip walls averaged 0.32 (max 0.93)
    // against a prescribed 0, and the inlet ran at 0.22 against a prescribed
    // 0.02. Those unconstrained boundary nodes are where the run destabilises.
    //
    // Used in three places that must agree: the momentum rows are pinned, the
    // nodal values are written from the boundary condition, and the pressure
    // step's velocity correction skips these nodes (a prescribed velocity has
    // zero correction -- correcting it is what broke earlier attempts).
    stk::mesh::PartVector velocityDirichletParts(
        const std::shared_ptr<domain> domain,
        label phaseIndex) const;

    void applyPhaseVelocityDirichlet(const std::shared_ptr<domain> domain,
                                     label phaseIndex);

    // Euler-Euler pressure update.
    //
    // This intentionally hides flowModel::updatePressure. The inherited
    // version routes an `opening` boundary through
    // flowModel::updatePressureBoundarySideFieldOpening_, which builds the
    // opening total-pressure side field from the *mixture* URef()/rhoRef()/
    // mDotRef(). An Euler-Euler realm never creates those bulk fields -- only
    // the per-phase ones exist -- so the inherited path dereferences null
    // fields. The quantities are taken from a configured reference phase, or
    // reconstructed from the phases when none is configured. Call sites hold
    // an eulerEulerModel*, so this override is selected statically.
    void updatePressure(const std::shared_ptr<domain> domain);

private:
    // Opening total-pressure side field. A configured reference phase supplies
    // rho, U, and mDot directly; otherwise reconstruct mixture quantities:
    //   rho_m   = sum_k alpha_k rho_k
    //   U_m     = sum_k alpha_k rho_k U_k / rho_m   (mass averaged)
    //   mDot_m  = sum_k mDot_k                      (sign selects in/outflow)
    void updatePhasicPressureSideFieldOpening_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updatePhasicPressureSideFields_(
        const std::shared_ptr<domain> domain);

    void setupEulerEulerAuxiliaryFields_(
        const std::shared_ptr<domain> domain,
        label phaseIndex);

    void updatePhaseMassFluxInterior_(const std::shared_ptr<domain> domain,
                                      label phaseIndex,
                                      bool includeRhieChow);

    void updatePhaseMassFluxSideParts_(
        const std::shared_ptr<domain> domain,
        label phaseIndex,
        const stk::mesh::PartVector& parts);

    // Keep each phase's boundary reversal flag consistent with its own mass
    // flux. In particular, a zero-gradient outlet must not assemble a
    // negative implicit advective coefficient when that phase reverses: the
    // established bulk-flow path turns such a face into a temporary slip
    // wall and clips its outgoing flux to zero. Openings retain the signed
    // flux and use the flag to select their entrainment values.
    void updatePhaseFlowReversalFlags_(
        const std::shared_ptr<domain> domain,
        label phaseIndex);
};

} // namespace accel

#endif // EULEREULERMODEL_H
