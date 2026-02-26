// File       :
// turbulentIntermittencyCorrelationTransitionSSTAssemblerNodeTerms.cpp
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: Node terms for gamma
// equation in correlation-based
// transition SST model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "turbulentIntermittencyCorrelationTransitionSSTAssembler.h"

namespace accel
{

void turbulentIntermittencyCorrelationTransitionSSTAssembler::
    assembleNodeTermsFusedSteady_(const domain* domain, Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants (Menter 2015)
    const scalar flength = model_->flength();
    const scalar caTwo = model_->caTwo();
    const scalar ceTwo = model_->ceTwo();
    const scalar Ctu1 = model_->Ctu1();
    const scalar Ctu2 = model_->Ctu2();
    const scalar Ctu3 = model_->Ctu3();
    const scalar CPG1 = model_->CPG1();
    const scalar CPG2 = model_->CPG2();
    const scalar CPG3 = model_->CPG3();
    const scalar CPG1_lim = model_->CPG1_lim();
    const scalar CPG2_lim = model_->CPG2_lim();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // other
    scalar dt = model_->controlsRef().getPhysicalTimescale();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // false transient
            scalar lhsfac = rho * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // compute strain rate magnitude
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // compute vorticity magnitude (for dissipation term)
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar vorticity = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                    dudx[SPATIAL_DIM * j + i]);
                    wijMag += vorticity * vorticity;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute turbulence intensity Tu
            scalar Tu = std::min(
                100.0 * std::sqrt(2.0 / 3.0 * k) / (omega * y + SMALL), 100.0);

            // Compute pressure gradient parameter lambda_0L
            scalar dvnn = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                dvnn += dudx[SPATIAL_DIM * i + i];
            }
            scalar lambda0L = -7.57e-3 * dvnn * y * y * rho / mu + 0.0128;
            lambda0L = std::max(std::min(lambda0L, 1.0), -1.0);

            // Compute FPG and apply to Tu (FPG multiplies Tu)
            scalar fpg = correlationTransitionShearStressTransportModel::FPG(
                lambda0L, CPG1, CPG2, CPG3, CPG1_lim, CPG2_lim);

            // Compute critical Reynolds number Re_theta_c
            scalar Re0c = Ctu1 + Ctu2 * std::exp(-Ctu3 * Tu * fpg);

            // Compute Rev and Fonset
            scalar Rev = rho * y * y * sijMag / mu;
            scalar fOnset1 = Rev / (2.2 * Re0c + SMALL);
            scalar fOnset2 = std::min(fOnset1, 2.0);
            scalar Rt = rho * k / (omega * mu + SMALL);
            scalar fOnset3 = std::max(1.0 - std::pow(Rt / 3.5, 3.0), 0.0);
            scalar fOnset = std::max(fOnset2 - fOnset3, 0.0);

            // Compute Fturb
            scalar fTurb = std::exp(-Rt * Rt * Rt * Rt / 16.0);

            // Production term (no sqrt variant)
            const scalar Pgamma =
                flength * rho * sijMag * fOnset * gamma * (1.0 - gamma);

            // Dissipation term
            const scalar Dgamma =
                caTwo * rho * wijMag * fTurb * gamma * (ceTwo * gamma - 1.0);

            // Linearization (positivity-preserving)
            const scalar PgammaDir =
                flength * rho * sijMag * fOnset * (1.0 - gamma);
            const scalar PgammaDirP = -flength * rho * sijMag * fOnset;
            const scalar DgammaDir =
                caTwo * rho * wijMag * fTurb * (ceTwo * gamma - 1.0);
            const scalar DgammaDirP = caTwo * rho * wijMag * fTurb * ceTwo;

            const scalar gamma_pos1 = std::max(DgammaDir - PgammaDir, 0.0);
            const scalar gamma_pos2 = std::max(DgammaDirP - PgammaDirP, 0.0);

            rhs[0] += (Pgamma - Dgamma) * vol;
            lhs[0] += (gamma_pos1 + gamma_pos2 * gamma) * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentIntermittencyCorrelationTransitionSSTAssembler::
    assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                              Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants (Menter 2015)
    const scalar flength = model_->flength();
    const scalar caTwo = model_->caTwo();
    const scalar ceTwo = model_->ceTwo();
    const scalar Ctu1 = model_->Ctu1();
    const scalar Ctu2 = model_->Ctu2();
    const scalar Ctu3 = model_->Ctu3();
    const scalar CPG1 = model_->CPG1();
    const scalar CPG2 = model_->CPG2();
    const scalar CPG3 = model_->CPG3();
    const scalar CPG1_lim = model_->CPG1_lim();
    const scalar CPG2_lim = model_->CPG2_lim();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->meshRef().controlsRef().getTimestep();
    const auto c = BDF1::coeff();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gammaOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOld, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * gamma + lhsfacOld * gammaOld);

            // compute strain rate magnitude
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // compute vorticity magnitude (for dissipation term)
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar vorticity = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                    dudx[SPATIAL_DIM * j + i]);
                    wijMag += vorticity * vorticity;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute turbulence intensity Tu
            scalar Tu = std::min(
                100.0 * std::sqrt(2.0 / 3.0 * k) / (omega * y + SMALL), 100.0);

            // Compute pressure gradient parameter lambda_0L
            scalar dvnn = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                dvnn += dudx[SPATIAL_DIM * i + i];
            }
            scalar lambda0L = -7.57e-3 * dvnn * y * y * rho / mu + 0.0128;
            lambda0L = std::max(std::min(lambda0L, 1.0), -1.0);

            // Compute FPG and apply to Tu (FPG multiplies Tu)
            scalar fpg = correlationTransitionShearStressTransportModel::FPG(
                lambda0L, CPG1, CPG2, CPG3, CPG1_lim, CPG2_lim);

            // Compute critical Reynolds number Re_theta_c
            scalar Re0c = Ctu1 + Ctu2 * std::exp(-Ctu3 * Tu * fpg);

            // Compute Rev and Fonset
            scalar Rev = rho * y * y * sijMag / mu;
            scalar fOnset1 = Rev / (2.2 * Re0c + SMALL);
            scalar fOnset2 = std::min(fOnset1, 2.0);
            scalar Rt = rho * k / (omega * mu + SMALL);
            scalar fOnset3 = std::max(1.0 - std::pow(Rt / 3.5, 3.0), 0.0);
            scalar fOnset = std::max(fOnset2 - fOnset3, 0.0);

            // Compute Fturb
            scalar fTurb = std::exp(-Rt * Rt * Rt * Rt / 16.0);

            // Production term (no sqrt variant)
            const scalar Pgamma =
                flength * rho * sijMag * fOnset * gamma * (1.0 - gamma);

            // Dissipation term
            const scalar Dgamma =
                caTwo * rho * wijMag * fTurb * gamma * (ceTwo * gamma - 1.0);

            // Linearization (positivity-preserving)
            const scalar PgammaDir =
                flength * rho * sijMag * fOnset * (1.0 - gamma);
            const scalar PgammaDirP = -flength * rho * sijMag * fOnset;
            const scalar DgammaDir =
                caTwo * rho * wijMag * fTurb * (ceTwo * gamma - 1.0);
            const scalar DgammaDirP = caTwo * rho * wijMag * fTurb * ceTwo;

            const scalar gamma_pos1 = std::max(DgammaDir - PgammaDir, 0.0);
            const scalar gamma_pos2 = std::max(DgammaDirP - PgammaDirP, 0.0);

            rhs[0] += (Pgamma - Dgamma) * vol;
            lhs[0] += (gamma_pos1 + gamma_pos2 * gamma) * vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * gamma * divUm * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentIntermittencyCorrelationTransitionSSTAssembler::
    assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                               Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants (Menter 2015)
    const scalar flength = model_->flength();
    const scalar caTwo = model_->caTwo();
    const scalar ceTwo = model_->ceTwo();
    const scalar Ctu1 = model_->Ctu1();
    const scalar Ctu2 = model_->Ctu2();
    const scalar Ctu3 = model_->Ctu3();
    const scalar CPG1 = model_->CPG1();
    const scalar CPG2 = model_->CPG2();
    const scalar CPG3 = model_->CPG3();
    const scalar CPG1_lim = model_->CPG1_lim();
    const scalar CPG2_lim = model_->CPG2_lim();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->meshRef().controlsRef().getTimestep();
    const auto c = BDF2::coeff(dt, mesh.controlsRef().getTimestep(-1));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);

            scalar gammaOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOld, node);
            scalar gammaOldOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOldOld, node);

            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);

            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * gamma + lhsfacOld * gammaOld +
                       lhsfacOldOld * gammaOldOld);

            // compute strain rate magnitude
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // compute vorticity magnitude (for dissipation term)
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar vorticity = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                    dudx[SPATIAL_DIM * j + i]);
                    wijMag += vorticity * vorticity;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute turbulence intensity Tu
            scalar Tu = std::min(
                100.0 * std::sqrt(2.0 / 3.0 * k) / (omega * y + SMALL), 100.0);

            // Compute pressure gradient parameter lambda_0L
            scalar dvnn = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                dvnn += dudx[SPATIAL_DIM * i + i];
            }
            scalar lambda0L = -7.57e-3 * dvnn * y * y * rho / mu + 0.0128;
            lambda0L = std::max(std::min(lambda0L, 1.0), -1.0);

            // Compute FPG and apply to Tu (FPG multiplies Tu)
            scalar fpg = correlationTransitionShearStressTransportModel::FPG(
                lambda0L, CPG1, CPG2, CPG3, CPG1_lim, CPG2_lim);

            // Compute critical Reynolds number Re_theta_c
            scalar Re0c = Ctu1 + Ctu2 * std::exp(-Ctu3 * Tu * fpg);

            // Compute Rev and Fonset
            scalar Rev = rho * y * y * sijMag / mu;
            scalar fOnset1 = Rev / (2.2 * Re0c + SMALL);
            scalar fOnset2 = std::min(fOnset1, 2.0);
            scalar Rt = rho * k / (omega * mu + SMALL);
            scalar fOnset3 = std::max(1.0 - std::pow(Rt / 3.5, 3.0), 0.0);
            scalar fOnset = std::max(fOnset2 - fOnset3, 0.0);

            // Compute Fturb
            scalar fTurb = std::exp(-Rt * Rt * Rt * Rt / 16.0);

            // Production term (no sqrt variant)
            const scalar Pgamma =
                flength * rho * sijMag * fOnset * gamma * (1.0 - gamma);

            // Dissipation term
            const scalar Dgamma =
                caTwo * rho * wijMag * fTurb * gamma * (ceTwo * gamma - 1.0);

            // Linearization (positivity-preserving)
            const scalar PgammaDir =
                flength * rho * sijMag * fOnset * (1.0 - gamma);
            const scalar PgammaDirP = -flength * rho * sijMag * fOnset;
            const scalar DgammaDir =
                caTwo * rho * wijMag * fTurb * (ceTwo * gamma - 1.0);
            const scalar DgammaDirP = caTwo * rho * wijMag * fTurb * ceTwo;

            const scalar gamma_pos1 = std::max(DgammaDir - PgammaDir, 0.0);
            const scalar gamma_pos2 = std::max(DgammaDirP - PgammaDirP, 0.0);

            rhs[0] += (Pgamma - Dgamma) * vol;
            lhs[0] += (gamma_pos1 + gamma_pos2 * gamma) * vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * gamma * divUm * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
