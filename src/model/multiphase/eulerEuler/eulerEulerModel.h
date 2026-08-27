// File       : eulerEulerModel.h
// Description: Inhomogeneous Euler-Euler multiphase model infrastructure

#ifndef EULEREULERMODEL_H
#define EULEREULERMODEL_H

#include "elementScalarField.h"
#include "multiphaseModel.h"

#include <map>
#include <memory>
#include <utility>

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

    // Pressure-response coefficient after local phase/drag coupling.  This
    // replaces the uncoupled du coefficient in the Euler-Euler pressure and
    // flux corrections; free-surface fields do not use it.
    nodeField<1>& coupledDuRef(label phaseIndex);
    const nodeField<1>& coupledDuRef(label phaseIndex) const;

    // Build the local OpenFOAM-style phase/drag pressure response matrix.
    void updateCoupledPressureResponse(const std::shared_ptr<domain> domain);

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

    // Volume fraction passed through the same redistribution and the same
    // limiter as the body force, for the drag law to use.
    //
    // The terminal slip follows K u_r = F_buoy. updatePhaseBodyForces builds
    // F_buoy from the *redistributed* alpha_k, while K is proportional to the
    // *pointwise* alpha_k, so u_r^2 ends up proportional to the ratio of the
    // two and inherits alpha's grid-scale content -- although the terminal
    // slip is physically independent of alpha altogether. On the laminar
    // bubble column alpha's 2*dx component is 30% of its mean, which predicts
    // a slip mode of 0.044 against a mean slip of 0.295; the measured value is
    // 0.046, while OpenFOAM, which uses one alpha for both, has 2e-5.
    //
    // Feeding the same smoothed alpha to the drag makes that ratio exactly one
    // (the body-force limiter clips the force by +/-|raw|, which for a
    // constant (rho_k-rho_ref)g is the same relative clip as limiting alpha
    // by +/-|alpha|), so the slip stops responding to alpha's noise while its
    // mean value is unchanged.
    void updatePhaseSmoothedVolumeFraction(
        const std::shared_ptr<domain> domain);

    nodeField<1>& phaseSmoothedAlphaRef(label phaseIndex);

    void updateInterphaseMomentumSources(
        const std::shared_ptr<domain> domain);

    void updatePhaseMassFlux(const std::shared_ptr<domain> domain,
                             label phaseIndex,
                             bool includeRhieChow = true);

    // Refresh alpha weighting without reconstructing the intrinsic flux that
    // pressure correction has already made conservative.
    void updatePhaseMassFluxWeighting(
        const std::shared_ptr<domain> domain,
        label phaseIndex);

    // Add the pressure-correction contribution to the predictor phase flux.
    void correctPhaseMassFlux(const std::shared_ptr<domain> domain,
                              label phaseIndex);

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

    // Flux-corrected (Zalesak/MULES-style) alpha update: replaces the
    // implicit predictor's diagonal boundedness stabilizer and the global
    // clip-and-redistribute step with a locally bounded, exactly
    // conservative correction. See phasicContinuityEquation::solve() for
    // the call site and domain->multiphase_.eulerEulerFCT_.fluxCorrectedTransport_
    // for the switch that enables it.
    //
    // Ported from freeSurfaceFlowModel's correctFCT -- same algorithm, same
    // CVFEM SCS-face structure (MasterElementRepo/adjacentNodes()).
    //
    // The high-order flux IS limited (van Leer TVD, see computeFctFH_).
    // An earlier version used a raw shape-function/central FH on the theory
    // that the MULES layer alone would bound it; that was wrong and is the
    // single thing that made this path destabilising. Measured: ~90% of
    // faces ran with lambda == 1, so the limiter never engaged and alpha
    // was advected by pure central differencing, growing grid-scale
    // oscillations that drove max|U.air| from 0.4 to 6 m/s by t=3.5 and
    // then to overflow. OpenFOAM's alpha equation is `Gauss vanLeer` with
    // MULES as a second layer; freeSurfaceFlowModel blends with its NVD
    // factor. A bounded high-order flux is a requirement, not a refinement.
    //
    // Still deliberately not ported: the VOF interface-compression term
    // (needs an interface-normal field Euler-Euler has no equivalent of).
    // FL is a genuine upwind reconstruction (computeFctFL_) -- mDotRef is
    // NOT the low-order upwind flux it was first assumed to be (verified:
    // reusing it made A == 0 everywhere).
    void correctVolumeFractionFCT(const std::shared_ptr<domain> domain,
                                  label phaseIndex);

private:
    struct pairDragField
    {
        label phaseA;
        label phaseB;
        std::unique_ptr<nodeField<1>> coefficient;
    };

    nodeField<1>& pairDragRef(label phaseA, label phaseB);

    std::vector<std::unique_ptr<nodeField<1>>> phaseSmoothedAlpha_;
    std::vector<std::unique_ptr<nodeField<1>>> coupledDu_;
    std::vector<pairDragField> pairDragFields_;

    // Flux-corrected transport scratch fields, one set per phase.
    std::vector<std::unique_ptr<elementScalarField>> fctFL_;
    std::vector<std::unique_ptr<elementScalarField>> fctFH_;
    std::vector<std::unique_ptr<elementScalarField>> fctA_;
    std::vector<std::unique_ptr<elementScalarField>> fctLambda_;
    std::vector<std::unique_ptr<nodeField<1>>> fctQPlus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctQMinus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctPPlus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctPMinus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctSumAPlus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctSumAMinus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctLimiterPlus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctLimiterMinus_;
    std::vector<std::unique_ptr<nodeField<1>>> fctLocalAlphaMax_;
    std::vector<std::unique_ptr<nodeField<1>>> fctLocalAlphaMin_;

    // Index of phaseIndex within phases_/the fct*_ vectors above.
    label fctLocalPhaseIndex_(label phaseIndex) const;

    // Low-order flux: intrinsicMDot (alpha-independent) times the upwind
    // nodal alpha at each SCS integration point (alpha[left] if the
    // intrinsic flux is positive, alpha[right] otherwise) -- the actual
    // low-order/bounded candidate FCT needs. Not mDotRef(phaseIndex):
    // that field turned out to already be a shape-function-weighted
    // (central-like) reconstruction, not upwind -- see class-level comment.
    void computeFctFL_(const std::shared_ptr<domain> domain,
                       label phaseIndex);

    // High-order candidate flux: intrinsicMDot times a van Leer TVD
    // reconstruction of alpha at each SCS integration point (NOT raw
    // central). A bounded high-order flux is a requirement here, not a
    // refinement: with raw central, ~90% of faces ran with limiter
    // lambda == 1, so the MULES layer never engaged and alpha was advected
    // by pure central differencing -- unbounded, and it grew grid-scale
    // oscillations that blew the phase velocities up. OpenFOAM uses
    // `div(phi,alpha) Gauss vanLeer` for exactly this reason, with MULES
    // only as a second safety layer; freeSurfaceFlowModel blends with its
    // NVD factor. See the .cpp for the reconstruction itself.
    void computeFctFH_(const std::shared_ptr<domain> domain,
                       label phaseIndex);

    // Antidiffusive correction A = FH - FL.
    void computeFctA_(const std::shared_ptr<domain> domain,
                      label phaseIndex);

    // Per-node max allowable change Q+/Q-, from local (neighbour) alpha
    // extrema -- the max/min of each node's directly-connected neighbours
    // via SCS faces, clipped to the global [residualVolumeFraction,
    // 1-residualVolumeFraction] bounds -- following the corrector-MULES
    // form: alpha here is already the implicit predictor's (bounded)
    // result. This mirrors nodeField::updateMinMaxFields's own neighbour-
    // scan, reimplemented directly rather than called: that accessor
    // depends on setupMinMaxFields() having been invoked, which only
    // happens automatically for the highResolution advection scheme (see
    // the .cpp history/comments) -- Euler-Euler's alpha uses upwind, so
    // those fields are never constructed there. A first version used fixed
    // global bounds only, which was NOT enough on its own: without a local
    // cap, alpha near the inlet could jump straight to the global ceiling
    // in one step and then cascade outward, saturating the whole domain --
    // exactly the failure a local neighbour bound is supposed to prevent.
    //
    // The neighbour scan also folds in boundary/opening faces, using the
    // same inflow/outflow alpha selection updatePhaseMassFluxSideParts_
    // already uses -- a boundary-adjacent node's bound needs to see what's
    // actually entering there, not just its interior neighbours, or its
    // own correction stays capped by a ceiling that ignores the one real
    // "neighbour" driving it (measured cause of a large, unphysical
    // domain-mean drift when this was interior-only).
    void computeFctQ_(const std::shared_ptr<domain> domain,
                      label phaseIndex);

    // Per-node unlimited pending antidiffusive flux P+/P- (computed once,
    // from the unlimited A).
    void computeFctP_(const std::shared_ptr<domain> domain,
                      label phaseIndex);

    // One Zalesak limiter iteration: lambda+/- = Q/P on iter 0, or
    // (sumA+Q)/P on later iterations (already-limited flux frees up
    // capacity on the opposite sign). Returns the max change, for an
    // optional early exit.
    scalar computeFctLambdaNode_(const std::shared_ptr<domain> domain,
                                 label phaseIndex,
                                 label iter);

    // Interpolate the per-node limiter to each SCS face:
    // lambda_ip = min(lambda-[donor], lambda+[receiver]).
    void computeFctLambdaIP_(const std::shared_ptr<domain> domain,
                             label phaseIndex);

    // Accumulate the limited flux actually used this iteration, so the next
    // iteration's lambda computation can see how much capacity it consumed.
    void computeFctSumA_(const std::shared_ptr<domain> domain,
                         label phaseIndex);

    // Apply the fully-limited correction:
    // alpha -= dt/(rho*V) * sum_faces(lambda*A).
    void applyFctCorrection_(const std::shared_ptr<domain> domain,
                             label phaseIndex);

    // Final global safety net after the locally-conservative FCT
    // correction above. applyFctCorrection_ is proven exactly conservative
    // *for whatever alpha it's handed* (checked directly: total domain
    // alpha*volume before/after a single call matches to float precision)
    // -- but that says nothing about whether the alpha it's handed, straight
    // out of the implicit predictor solve, already tracks the true physical
    // total. The old clip-and-redistribute step this whole path replaces
    // also forced the domain total to exactly match
    // oldMass - dt*boundaryRate (the value boundary in/outflow implies),
    // independent of what the assembled matrix's own boundary handling
    // produced -- correcting for whatever approximation error is in that,
    // which otherwise compounds over thousands of steps. Reuses that same
    // computation and weighted redistribution here, as a residual correction
    // on top of the FCT result rather than the primary mechanism (measured
    // necessary: FCT alone left the domain drifting to ~95% water by t=8 on
    // bubbleColumnLaminar, vs OpenFOAM's flat ~69%).
    void applyGlobalMassSafetyNet_(const std::shared_ptr<domain> domain,
                                   label phaseIndex);

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
