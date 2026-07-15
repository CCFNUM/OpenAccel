// File       : segregatedFlowEquations.cpp
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Fabian Wermelinger
// Description: Segregated Navier-Stokes equations implementation details
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "segregatedFlowEquations.h"
#include "realm.h"

namespace accel
{

segregatedFlowEquations::segregatedFlowEquations(realm* realm)
    : equation("Segregated Flow"), flowModel(realm)
{
    U_eq_ = std::make_unique<navierStokesEquation>(realm, this);
    pCorr_eq_ = std::make_unique<pressureCorrectionEquation>(realm, this);
    U_eq_->setFallbackName(canonicalEquationIDName(this->equation::name_));
    pCorr_eq_->setFallbackName(canonicalEquationIDName(this->equation::name_));

    // set relaxation factor for mass flux field to 0.75 for steady-state cases,
    // if and only if, not specified by user
    if (controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.relaxMass_ < 0.9999)
    {
        this->mDotRef().setURF(
            this->controlsRef()
                .solverRef()
                .solverControl_.basicSettings_.convergenceControl_
                .relaxationParameters_.relaxMass_);
    }

    // set sub-iterations
    subIters_ = this->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.segregatedFlow_;
}

void segregatedFlowEquations::addDomain(std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    U_eq_->addDomain(domain);
    pCorr_eq_->addDomain(domain);
}

bool segregatedFlowEquations::isConverged() const
{
    return U_eq_->isConverged() && pCorr_eq_->isConverged();
}

void segregatedFlowEquations::setup()
{
    U_eq_->setup();
    pCorr_eq_->setup();

    equation::isCreated_ = U_eq_->isCreated() && pCorr_eq_->isCreated();
}

void segregatedFlowEquations::initialize()
{
    U_eq_->initialize();
    pCorr_eq_->initialize();
}

void segregatedFlowEquations::postInitialize()
{
    U_eq_->postInitialize();
    pCorr_eq_->postInitialize();

    equation::isInitialized_ =
        U_eq_->isInitialized() && pCorr_eq_->isInitialized();
}

void segregatedFlowEquations::solve()
{
    // frozen flow ... do not assemble/solve for flow equations
    if (controlsRef().solverRef().solverControl_.expertParameters_.freezeFlow_)
    {
        return;
    }

    // predictor step: solve momentum
    {
        U_eq_->preSolve();
        U_eq_->solve();
        U_eq_->postSolve();
    }

    // freeze_pressure: momentum predictor only, skip the pressure corrector.
    if (controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.freezePressure_)
    {
        this->URef().updateScale();
        FOREACH_DOMAIN(updateVelocityGradientField);
        FOREACH_DOMAIN(updateVelocityBlendingFactorField);
        return;
    }

    // corrrector step: solve pressure correction
    {
        if (messager::master() && pCorr_eq_->subIters() > 1)
        {
            std::cout << std::endl
                      << "Sub-iterating " + pCorr_eq_->name() << ":\n";
        }

        for (label subIter = 1; subIter <= pCorr_eq_->subIters(); subIter++)
        {
            // only print sub-iter if it is active (> 1)
            if (messager::master() && pCorr_eq_->subIters() > 1)
            {
                std::cout << std::endl << " sub-iter: " << subIter << "\n";
            }

            pCorr_eq_->preSolve();
            pCorr_eq_->solve();
            pCorr_eq_->postSolve();

            // correct density (only if compressible)
            FOREACH_DOMAIN_IF(updateDensity, domain->isMaterialCompressible());

            // update shock-sensor damping of the high-resolution blending
            // factor
            FOREACH_DOMAIN(updateBlendingDampingFactorField_);

            // update density-related fields
            FOREACH_DOMAIN_IF(updateDensityGradientField,
                              domain->isMaterialCompressible());
            FOREACH_DOMAIN_IF(updateDensityBlendingFactorField,
                              domain->isMaterialCompressible());

            // correct mass flux and do updates:
            // 1) mass divergence
            // 2) flow reversal flag side field
            FOREACH_DOMAIN(updateSideMassFlowRateFraction);
            FOREACH_DOMAIN(transformMassFlowRateToAbsolute);
            FOREACH_DOMAIN(updateMassFlowRate);
            FOREACH_DOMAIN(transformMassFlowRateToRelative);
            FOREACH_DOMAIN(updateFlowReversalFlag);
            FOREACH_DOMAIN(updateMassDivergenceField);

            // velocity Correction Step (Pressure-Velocity Coupling)
            //
            // Adjusts velocity to satisfy the continuity equation at
            // sub-control volume faces: v_corr = -D * grad(p'), where D is the
            // momentum influence coefficient (du, or duTilde for SIMPLEC) and
            // grad(p') is the gradient of the dedicated pressure correction
            // field. That gradient was built during the pressure correction
            // step with gradient relaxation forced to 1, so it carries the
            // full (un-relaxed) increment regardless of how the pressure field
            // gradient itself is relaxed.
            FOREACH_DOMAIN_RAW({
                stk::mesh::MetaData& metaData = meshRef().metaDataRef();

                // SIMPLE-Consistent
                const bool consistent =
                    controlsRef()
                        .solverRef()
                        .solverControl_.expertParameters_.consistent_;

                const auto* duSTKFieldPtr =
                    consistent
                        ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::duTilde_ID)
                        : metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::du_ID);

                const auto& gradPCorrSTKFieldRef =
                    pCorrRef().gradRef().stkFieldRef();
                auto& USTKFieldRef = URef().stkFieldRef();

                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
                stk::mesh::BucketVector const& nodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, nodeBucket);
                    const scalar* dub =
                        stk::mesh::field_data(*duSTKFieldPtr, nodeBucket);
                    const scalar* dpCorrdxb =
                        stk::mesh::field_data(gradPCorrSTKFieldRef, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            Ub[SPATIAL_DIM * iNode + i] -=
                                dub[SPATIAL_DIM * iNode + i] *
                                dpCorrdxb[SPATIAL_DIM * iNode + i];
                        }
                    }
                }
            });

            // update pressure gradient for the next momentum solve (now free to
            // use the default relaxation factor)
            FOREACH_DOMAIN(updatePressureGradientField);

            // update velocity scale: raw
            this->URef().updateScale();
            this->updatePressureScale();

            // update velocity gradient
            FOREACH_DOMAIN(updateVelocityGradientField);

            // update velocity high-res fields
            FOREACH_DOMAIN(updateVelocityBlendingFactorField);

            // if converged .. break sub-iter loop
            if (pCorr_eq_->isConverged())
            {
                break;
            }
        }
    }
}

void segregatedFlowEquations::postSolve()
{
    // calculate post-processed quantities

    FOREACH_DOMAIN(updateCourantNumberField_);
    FOREACH_DOMAIN(updateMachNumberField_);
    FOREACH_DOMAIN(updateTotalPressureField_);
    FOREACH_DOMAIN(updateRelativeVelocityField_);
    FOREACH_DOMAIN(updateUWallCoeffs);     // laminar
    FOREACH_DOMAIN(updateWallShearStress); // laminar
    FOREACH_DOMAIN(updateMassImbalance_);
    FOREACH_DOMAIN(updateInterfaceMassImbalance_);
    FOREACH_DOMAIN(updateInterfaceMomentumImbalance_);

    this->reportFlowData_();
}

void segregatedFlowEquations::preTimeStep()
{
    resetCourantNumber();

    U_eq_->preTimeStep();
    pCorr_eq_->preTimeStep();
}

void segregatedFlowEquations::printScales()
{
    U_eq_->printScales();
    pCorr_eq_->printScales();
}

void segregatedFlowEquations::updateMassFlowRate(
    const std::shared_ptr<domain> domain)
{
    // update according to rhie-chow
    fieldBroker::updateMassFlowRate(domain);
}

} /* namespace accel */
