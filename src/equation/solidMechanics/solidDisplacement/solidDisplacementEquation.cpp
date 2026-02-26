// File       : solidDisplacementEquation.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidDisplacementEquation.h"

namespace accel
{

solidDisplacementEquation::solidDisplacementEquation(realm* realm)
    : equation("Solid Displacement"), solidMechanicsModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    // Collect equation names
    std::vector<std::string> eqNames;

#if SPATIAL_DIM == 2
    eqNames.push_back("x-Disp");
    eqNames.push_back("y-Disp");
#elif SPATIAL_DIM == 3
    eqNames.push_back("x-Disp");
    eqNames.push_back("y-Disp");
    eqNames.push_back("z-Disp");
#endif

    // set
    this->setEquationName(eqNames);

    // set relaxation factor for solid displacement field. Default is 1
    this->DRef().setURF(
        controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.solidDisplacementRelaxationFactor_);

    // set sub-iterations
    subIters_ = this->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.solidDisplacement_;

    // Initialize Aitken acceleration parameters from config
    const auto& accel = controlsRef()
                            .solverRef()
                            .solverControl_.advancedOptions_.equationControls_
                            .acceleration_.solidDisplacement_;
    useAitken_ = accel.aitkenEnabled_;
    aitkenOmegaInit_ = accel.aitkenInitialOmega_;
    aitkenOmegaMin_ = accel.aitkenOmegaMin_;
    aitkenOmegaMax_ = accel.aitkenOmegaMax_;
    aitkenOmega_ = aitkenOmegaInit_;
}

void solidDisplacementEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::solid);
}

bool solidDisplacementEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void solidDisplacementEquation::setup()
{
    // setup of fields on defined domains
    FOREACH_DOMAIN(setupDensity);
    FOREACH_DOMAIN(setupYoungModulus);
    FOREACH_DOMAIN(setupPoissonRatio);
    FOREACH_DOMAIN(setupDisplacement);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler: Gamma is dummy
    assembler_->setup(&DRef(), diffusion, domainVector_, 0.0);

    // setup linear solver
    // FIXME: Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), fieldBroker::meshRef());

    equation::isCreated_ = true;
}

void solidDisplacementEquation::initialize()
{
    // raw initialization of solid displacement
    FOREACH_DOMAIN(initializeDisplacement);

    // raw initialization of properties
    FOREACH_DOMAIN(initializeDensity);
    FOREACH_DOMAIN(initializeYoungModulus);
    FOREACH_DOMAIN(initializePoissonRatio);

    // update displacement gradient
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // calculate initial stress and strain
    FOREACH_DOMAIN(updateStressAndStrain_);

    // update enthalpy scale
    DRef().updateScale();
    rhoRef().updateScale();
    ERef().updateScale();
    nuRef().updateScale();
}

void solidDisplacementEquation::postInitialize()
{
    equation::isInitialized_ = true;
}

void solidDisplacementEquation::preSolve()
{
    // solid displacement updates
    FOREACH_DOMAIN(updateDisplacementPrevIterField);
    FOREACH_DOMAIN(updateDisplacement);

    // other property updates
    FOREACH_DOMAIN(updateDensity);
    FOREACH_DOMAIN(updateYoungModulus);
    FOREACH_DOMAIN(updatePoissonRatio);
}

void solidDisplacementEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    linearSystem::simulationRef().getProfiler().push("linear_system_assembly");

    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get());

    // fix system on all dirichlet boundaries
    assembler_->fix(this->collectDirichletBoundaryParts_(), {}, ctx.get());

    // fix system on all mixed boundaries: pass ignored dof's
    for (const auto& domain : domainVector_)
    {
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            switch (bcType)
            {
                case boundaryConditionType::mixed:
                    {
                        const scalar* fixedValueFlag =
                            (this->DRef()
                                 .boundaryConditionRef(domain->index(),
                                                       iBoundary)
                                 .template data<SPATIAL_DIM>(
                                     "fixed_value_flag"))
                                .value();

                        std::vector<label> ignoreDofs;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            if (fixedValueFlag[i] < 0.5)
                            {
                                ignoreDofs.push_back(i);
                            }
                        }
                        assembler_->fix(
                            domain->zonePtr()->boundaryPtr(iBoundary)->parts(),
                            {},
                            ctx.get(),
                            ignoreDofs);
                    }
                    break;

                default:
                    break;
            }
        }
    }

    linearSystem::simulationRef().getProfiler().pop();

    // solve linear system
    linearSystem::solve();

    // Compute relaxation factor (Aitken or static URF)
    const Vector& correction = ctx->getXVector();
    scalar relaxationFactor =
        useAitken_ ? computeAitkenOmega_(correction) : DRef().urf();

    // correction
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
            domain.get(),
            correction,
            stk::topology::NODE_RANK,
            DRef().stkFieldRef(),
            relaxationFactor);

        // synchronize
        DRef().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) update gradient
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // 2) update scales
    DRef().updateScale();
    rhoRef().updateScale();
    ERef().updateScale();
    nuRef().updateScale();
}

void solidDisplacementEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateDisplacementPrevTimeField);

    // Reset Aitken state for new timestep
    if (useAitken_)
    {
        aitkenIter_ = 0;
        aitkenOmega_ = aitkenOmegaInit_;
        aitkenResidualPrev_.clear();
    }

    // Reset FSI-level Aitken iteration state for new timestep.
    // fsiDfluidPrev_ is intentionally kept: it holds the converged interface
    // position from the last timestep and must serve as the starting point for
    // the first outer iteration, otherwise the fluid mesh snaps back to zero.
    fsiResidualPrev_
        .clear(); // cleared → signals "first outer iter of timestep"
    for (auto& [idx, omega] : fsiAitkenOmega_)
    {
        omega = fsiOmegaInit_;
    }
    if (fsiActive_)
    {
        fsiResidualNorm_ = 1.0;
    }
}

void solidDisplacementEquation::postSolve()
{
    // calculate post-processed quantities
    FOREACH_DOMAIN(updateStressAndStrain_);
}

stk::mesh::PartVector
solidDisplacementEquation::collectDirichletBoundaryParts_()
{
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            switch (bcType)
            {
                case boundaryConditionType::specifiedValue:
                    {
                        const auto& boundaryRef =
                            domain->zonePtr()->boundaryRef(iBoundary);

                        for (auto part : boundaryRef.parts())
                        {
                            incPartVec.push_back(part);
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    return incPartVec;
}

void solidDisplacementEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    const scalar disp_scale_inv = 1.0 / (this->DRef().scale() + ::accel::SMALL);
    for (int i = 0; i < SPATIAL_DIM; i++)
    {
        residual_scales[i] = disp_scale_inv;
    }
}

void solidDisplacementEquation::applyDependencyUpdates_(
    const domain* domain,
    const stk::mesh::EntityRank entityRank,
    STKScalarField& /*stk_dst*/)
{
    assert(entityRank == stk::mesh::EntityRank::NODE_RANK);

    // if the zone deforms, and a fluid-solid interface exists, transfer
    // displacement from solid side to fluid side
    if (domain->zonePtr()->meshDeforming())
    {
#ifdef HAS_INTERFACE
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                DRef().transfer(interf->index(),
                                dataTransferType::copy,
                                !interf->isMasterZone(domain->index()));

                // --- FSI-level Aitken relaxation ---
                // Mark FSI as active so convergence check is enabled
                fsiActive_ = true;

                // Identify fluid zone and its interface side info
                const label fluidZoneIndex =
                    interf->isMasterZone(domain->index())
                        ? interf->slaveZoneIndex()
                        : interf->masterZoneIndex();

                const interfaceSideInfo* fluidSide =
                    interf->interfaceSideInfoPtr(fluidZoneIndex);

                auto& DSTKFieldRef = DRef().stkFieldRef();
                const label interfIdx = interf->index();

                // Select fluid interface nodes
                stk::mesh::Selector selFluidNodes =
                    fieldBroker::meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(fluidSide->currentPartVec_);
                stk::mesh::BucketVector const& fluidNodeBuckets =
                    fieldBroker::meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selFluidNodes);

                // Count total local nodes on fluid interface
                size_t nTotal = 0;
                for (auto ib = fluidNodeBuckets.begin();
                     ib != fluidNodeBuckets.end();
                     ++ib)
                {
                    nTotal += (*ib)->size();
                }
                const size_t vecSize = nTotal * SPATIAL_DIM;

                // First EVER call for this interface (simulation start, t=0):
                // initialise the accumulated fluid position to zero.
                // On all subsequent timesteps fsiDfluidPrev_ already holds the
                // converged interface position from the previous timestep.
                const bool isVeryFirstCall =
                    fsiDfluidPrev_.find(interfIdx) == fsiDfluidPrev_.end();
                if (isVeryFirstCall)
                {
                    fsiDfluidPrev_[interfIdx].assign(vecSize, 0.0);
                    fsiAitkenOmega_[interfIdx] = fsiOmegaInit_;
                    fsiResidualNormMax_[interfIdx] = 0.0;
                }

                // First outer iteration of the current timestep:
                // fsiResidualPrev_ was cleared in preTimeStep(), so it has no
                // entry yet.  This flag is consistent across all MPI ranks.
                const bool isFirstIterOfTimestep =
                    fsiResidualPrev_.find(interfIdx) == fsiResidualPrev_.end();

                // Read D_solid values just transferred to fluid interface nodes
                std::vector<scalar> DSolid(vecSize);
                {
                    size_t offset = 0;
                    for (auto ib = fluidNodeBuckets.begin();
                         ib != fluidNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& b = **ib;
                        const scalar* Db =
                            stk::mesh::field_data(DSTKFieldRef, b);
                        for (size_t iNode = 0; iNode < b.size(); ++iNode)
                        {
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                DSolid[offset++] = Db[SPATIAL_DIM * iNode + i];
                            }
                        }
                    }
                }

                // Compute residual: r^k = D_solid - D_fluid_prev
                std::vector<scalar> rk(vecSize);
                for (size_t j = 0; j < vecSize; ++j)
                {
                    rk[j] = DSolid[j] - fsiDfluidPrev_[interfIdx][j];
                }

                // Update Aitken omega from the second outer iteration onward.
                // isFirstIterOfTimestep is consistent across all MPI ranks
                // (fsiResidualPrev_ cleared globally in preTimeStep).
                scalar& omega = fsiAitkenOmega_[interfIdx];
                if (!isFirstIterOfTimestep)
                {
                    // Accumulate locally; ranks with no fluid interface nodes
                    // contribute zero — all ranks still call the collective.
                    scalar denom = 0.0;
                    scalar dot = 0.0;
                    if (vecSize > 0)
                    {
                        const auto& rPrev = fsiResidualPrev_[interfIdx];
                        for (size_t j = 0; j < vecSize; ++j)
                        {
                            const scalar dr = rk[j] - rPrev[j];
                            denom += dr * dr;
                            dot += rPrev[j] * dr;
                        }
                    }
                    messager::sumReduce(denom);
                    messager::sumReduce(dot);
                    if (denom > SMALL)
                    {
                        scalar omegaNew = std::abs(-omega * dot / denom);
                        omega =
                            std::clamp(omegaNew, fsiOmegaMin_, fsiOmegaMax_);
                    }
                }

                // Compute relaxed interface position:
                // D_fluid_new = D_fluid_prev + omega * r^k
                std::vector<scalar> DFluidNew(vecSize);
                for (size_t j = 0; j < vecSize; ++j)
                {
                    DFluidNew[j] = fsiDfluidPrev_[interfIdx][j] + omega * rk[j];
                }

                // Overwrite D field on fluid interface nodes with relaxed value
                {
                    size_t offset = 0;
                    for (auto ib = fluidNodeBuckets.begin();
                         ib != fluidNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& b = **ib;
                        scalar* Db = stk::mesh::field_data(DSTKFieldRef, b);
                        for (size_t iNode = 0; iNode < b.size(); ++iNode)
                        {
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                Db[SPATIAL_DIM * iNode + i] =
                                    DFluidNew[offset++];
                            }
                        }
                    }
                }

                // Store history for next outer iteration
                fsiResidualPrev_[interfIdx] = rk;
                fsiDfluidPrev_[interfIdx] = DFluidNew;

                // Compute normalized FSI residual norm
                scalar normSq = 0.0;
                for (const scalar v : rk)
                {
                    normSq += v * v;
                }
                messager::sumReduce(normSq);
                const scalar norm = std::sqrt(normSq);
                fsiResidualNormMax_[interfIdx] =
                    std::max(fsiResidualNormMax_[interfIdx], norm);
                fsiResidualNorm_ =
                    norm / (fsiResidualNormMax_[interfIdx] + SMALL);

                // Open residual file on first ever call for this interface
                if (fsiResidualStreams_.find(interfIdx) ==
                    fsiResidualStreams_.end())
                {
                    initializeFsiResidualFile_(interfIdx, interf->name());
                }
                writeFsiResidualLine_(interfIdx, omega, fsiResidualNorm_);

                if (messager::master())
                {
                    std::cout << "  FSI Aitken [" << interf->name() << "]"
                              << "  omega=" << std::scientific
                              << std::setprecision(4) << omega
                              << "  |r|_norm=" << std::scientific
                              << std::setprecision(4) << fsiResidualNorm_
                              << std::endl;
                }
            }
        }
#endif /* HAS_INTERFACE */
    }
}

void solidDisplacementEquation::printScales()
{
    if (messager::master())
    {
        std::cout << this->DRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << DRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << DRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << DRef().scale() << std::endl
                  << std::endl;

        std::cout << rhoRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << rhoRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << rhoRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << rhoRef().scale() << std::endl
                  << std::endl;

        std::cout << ERef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << ERef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << ERef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << ERef().scale() << std::endl
                  << std::endl;

        std::cout << nuRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << nuRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << nuRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << nuRef().scale() << std::endl
                  << std::endl;
    }
}

void solidDisplacementEquation::initializeFsiResidualFile_(
    label interfIdx,
    const std::string& interfName)
{
    if (!messager::master())
        return;

    // Build base path following the linearSystem naming convention:
    // postProcessing/0/residuals/fsi_<interfName>_###.out
    std::string baseName = "fsi_" + interfName;
    std::replace(baseName.begin(), baseName.end(), ' ', '_');

    const fs::path baseFilePath =
        linearSystem::simulationRef().getResidualDirectory() / baseName;
    const size_t len = baseFilePath.string().size();

    int nfiles = 0;
    for (const auto& entry : fs::directory_iterator{
             linearSystem::simulationRef().getResidualDirectory()})
    {
        if (entry.path().string().compare(
                0, len, baseFilePath.string(), 0, len) == 0)
        {
            ++nfiles;
        }
    }

    std::ostringstream oss;
    oss << baseFilePath.string() << "_" << std::setfill('0') << std::setw(3)
        << ++nfiles << ".out";

    auto stream = std::make_shared<std::ofstream>(oss.str());
    assert(stream->is_open());

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    auto& fout = *stream;
    fout << "# Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << "# Git revision: " << accel::git_revision << '\n';
    fout << "# FSI coupling residual history — interface: " << interfName
         << '\n';
    fout << "# \n";
    fout << "# " << "global_iterations" << '\t' << "inner_iterations" << '\t'
         << "sim_time[s]" << '\t' << "aitken_omega" << '\t'
         << "fsi_residual_norm" << '\n';

    fsiResidualStreams_[interfIdx] = stream;
}

void solidDisplacementEquation::writeFsiResidualLine_(label interfIdx,
                                                      scalar omega,
                                                      scalar residualNorm)
{
    if (!messager::master())
        return;

    auto it = fsiResidualStreams_.find(interfIdx);
    if (it == fsiResidualStreams_.end())
        return;

    auto& fout = *(it->second);
    fout << linearSystem::simulationRef().getGlobalIterationCount() << '\t'
         << linearSystem::simulationRef().getIterationCount() << '\t'
         << std::setprecision(3) << std::scientific
         << linearSystem::simulationRef().getSimulationTime() << '\t'
         << std::setprecision(6) << std::scientific << omega << '\t'
         << std::setprecision(6) << std::scientific << residualNorm
         << std::endl;
}

scalar solidDisplacementEquation::computeAitkenOmega_(const Vector& correction)
{
    const size_t n = correction.size();

    if (aitkenIter_ == 0)
    {
        // First iteration: store correction, use initial omega
        aitkenResidualPrev_ = correction;
        aitkenIter_++;
        return aitkenOmegaInit_;
    }

    // Compute norms and dot products for Aitken formula
    // r_prev = previous correction, r_curr = current correction
    scalar normPrevSq = 0.0;  // |r_{k-1}|²
    scalar normCurrSq = 0.0;  // |r_k|²
    scalar dotPrevCurr = 0.0; // r_{k-1} · r_k

    for (size_t i = 0; i < n; ++i)
    {
        normPrevSq += aitkenResidualPrev_[i] * aitkenResidualPrev_[i];
        normCurrSq += correction[i] * correction[i];
        dotPrevCurr += aitkenResidualPrev_[i] * correction[i];
    }

    // MPI reduction for parallel runs
    messager::sumReduce(normPrevSq);
    messager::sumReduce(normCurrSq);
    messager::sumReduce(dotPrevCurr);

    // Compute Aitken relaxation factor:
    // ω_k = ω_{k-1} * (|r_{k-1}|² - r_{k-1}·r_k) / |Δr|²
    // where |Δr|² = |r_k - r_{k-1}|² = |r_k|² - 2*r_{k-1}·r_k + |r_{k-1}|²
    scalar normDrSq = normCurrSq - 2.0 * dotPrevCurr + normPrevSq;
    scalar numerator = normPrevSq - dotPrevCurr;

    if (normDrSq > SMALL)
    {
        // Compute Aitken omega using standard formula (preCICE-style)
        scalar omegaNew = aitkenOmega_ * (numerator / normDrSq);
        aitkenOmega_ = std::clamp(omegaNew, aitkenOmegaMin_, aitkenOmegaMax_);
    }

    // Store current correction for next iteration
    aitkenResidualPrev_ = correction;
    aitkenIter_++;

#ifndef NDEBUG
    // Print Aitken omega for monitoring
    if (messager::master())
    {
        std::cout << "Aitken omega: " << std::scientific << std::setprecision(4)
                  << aitkenOmega_ << " (iter " << aitkenIter_ << ")"
                  << std::endl;
    }
#endif

    return aitkenOmega_;
}

} /* namespace accel */
