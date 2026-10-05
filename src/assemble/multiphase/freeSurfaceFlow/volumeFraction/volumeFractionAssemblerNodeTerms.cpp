// File       : volumeFractionAssemblerNodeTerms.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Node based terms of the volume fraction equation: divergence
//              correction, (false) transient and the interphase mass transfer
//              (phase change) source, assembled in one node loop
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "volumeFractionAssembler.h"

namespace accel
{

void volumeFractionAssembler::assembleNodeTermsFused_(const domain* domain,
                                                      Context* ctx)
{
    auto& mesh = field_broker_->meshRef();
    if (mesh.controlsRef().isTransient())
    {
        auto scheme = mesh.controlsRef()
                          .solverRef()
                          .solverControl_.basicSettings_.transientScheme_;
        switch (scheme)
        {
            case transientSchemeType::firstOrderBackwardEuler:
                assembleNodeTermsFusedFirstOrderUnsteady_(domain, ctx);
                break;

            case transientSchemeType::secondOrderBackwardEuler:
                assembleNodeTermsFusedSecondOrderUnsteady_(domain, ctx);
                break;

            default:
                break;
        }
    }
    else
    {
        assembleNodeTermsFusedSteady_(domain, ctx);
    }
}

void volumeFractionAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    std::vector<scalar> lhs(1);
    std::vector<scalar> rhs(1);
    std::vector<label> scratchIds(1);
    std::vector<scalar> scratchVals(1);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // other
    scalar dt = field_broker_->controlsRef().getPhysicalTimescale();

    // mass transfer pairs of this domain (phase change source)
    const auto mtPairs = model_->massTransferPairs(domain);

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

        // field chunks in bucket
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            lhs[0] = 0.0;
            rhs[0] = 0.0;

            // get values of current node
            scalar rho = rhob[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // false transient: added later
            scalar lhsfac = rho * vol / dt;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }

            // divergence correction
            scalar phii = phib[iNode];
            rhs[0] -= -div * phii;

            // false transient
            lhs[0] += lhsfac;

            // interphase mass transfer (phase change) source: residual
            // vol (s mdot - alpha rho D) with s = +1 for the vapor, -1 for the
            // liquid and 0 for any other phase, D = mdot (1/rho_v - 1/rho_l);
            // the -alpha rho D part (mixture expansion) acts on every phase
            // and is implicit (positive diagonal) when D > 0
            for (const auto* pair : mtPairs)
            {
                scalar sign = 0.0;
                if (pair->vaporIndex_ == phaseIndex_)
                    sign = 1.0;
                else if (pair->liquidIndex_ == phaseIndex_)
                    sign = -1.0;
                const scalar mdot =
                    *stk::mesh::field_data(*pair->mdotSTKFieldPtr_, node);

                scalar D = 0.0;
                if (pair->model_->includeContinuitySource())
                {
                    const scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(pair->liquidIndex_).stkFieldPtr(),
                        node);
                    const scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(pair->vaporIndex_).stkFieldPtr(),
                        node);
                    D = massTransferModel::volumeSource(mdot, rhoL, rhoV);
                }

                rhs[0] += vol * (sign * mdot - phii * rho * D);
                if (D > 0.0)
                {
                    lhs[0] += vol * rho * D;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void volumeFractionAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    std::vector<scalar> lhs(1);
    std::vector<scalar> rhs(1);
    std::vector<label> scratchIds(1);
    std::vector<scalar> scratchVals(1);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = mesh.controlsRef().getTimestep();
    const auto c = BDF1::coeff();

    // mass transfer pairs of this domain (phase change source)
    const auto mtPairs = model_->massTransferPairs(domain);

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

        // field chunks in bucket
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);
        scalar* phibOld = stk::mesh::field_data(*phiSTKFieldPtrOld, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            lhs[0] = 0.0;
            rhs[0] = 0.0;

            // get values of current node
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: added later
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }

            scalar phii = phib[iNode];
            scalar phiOldi = phibOld[iNode];

            // divergence correction
            rhs[0] -= -div * phii;

            // transient
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * phii + lhsfacOld * phiOldi);

            // interphase mass transfer (phase change) source: residual
            // vol (s mdot - alpha rho D) with s = +1 for the vapor, -1 for the
            // liquid and 0 for any other phase, D = mdot (1/rho_v - 1/rho_l);
            // the -alpha rho D part (mixture expansion) acts on every phase
            // and is implicit (positive diagonal) when D > 0
            for (const auto* pair : mtPairs)
            {
                scalar sign = 0.0;
                if (pair->vaporIndex_ == phaseIndex_)
                    sign = 1.0;
                else if (pair->liquidIndex_ == phaseIndex_)
                    sign = -1.0;
                const scalar mdot =
                    *stk::mesh::field_data(*pair->mdotSTKFieldPtr_, node);

                scalar D = 0.0;
                if (pair->model_->includeContinuitySource())
                {
                    const scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(pair->liquidIndex_).stkFieldPtr(),
                        node);
                    const scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(pair->vaporIndex_).stkFieldPtr(),
                        node);
                    D = massTransferModel::volumeSource(mdot, rhoL, rhoV);
                }

                rhs[0] += vol * (sign * mdot - phii * rho * D);
                if (D > 0.0)
                {
                    lhs[0] += vol * rho * D;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void volumeFractionAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    std::vector<scalar> lhs(1);
    std::vector<scalar> rhs(1);
    std::vector<label> scratchIds(1);
    std::vector<scalar> scratchVals(1);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        this->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = mesh.controlsRef().getTimestep();
    const auto c = BDF2::coeff(dt, mesh.controlsRef().getTimestep(-1));

    // mass transfer pairs of this domain (phase change source)
    const auto mtPairs = model_->massTransferPairs(domain);

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

        // field chunks in bucket
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* rhobOldOld =
            stk::mesh::field_data(*rhoSTKFieldPtrOldOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);
        scalar* phibOld = stk::mesh::field_data(*phiSTKFieldPtrOld, nodeBucket);
        scalar* phibOldOld =
            stk::mesh::field_data(*phiSTKFieldPtrOldOld, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            lhs[0] = 0.0;
            rhs[0] = 0.0;

            // get values of current node
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar rhoOldOld = rhobOldOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: added later
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }

            scalar phii = phib[iNode];
            scalar phiOldi = phibOld[iNode];
            scalar phiOldOldi = phibOldOld[iNode];

            // divergence correction
            rhs[0] -= -div * phii;

            // transient
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * phii + lhsfacOld * phiOldi +
                       lhsfacOldOld * phiOldOldi);

            // interphase mass transfer (phase change) source: residual
            // vol (s mdot - alpha rho D) with s = +1 for the vapor, -1 for the
            // liquid and 0 for any other phase, D = mdot (1/rho_v - 1/rho_l);
            // the -alpha rho D part (mixture expansion) acts on every phase
            // and is implicit (positive diagonal) when D > 0
            for (const auto* pair : mtPairs)
            {
                scalar sign = 0.0;
                if (pair->vaporIndex_ == phaseIndex_)
                    sign = 1.0;
                else if (pair->liquidIndex_ == phaseIndex_)
                    sign = -1.0;
                const scalar mdot =
                    *stk::mesh::field_data(*pair->mdotSTKFieldPtr_, node);

                scalar D = 0.0;
                if (pair->model_->includeContinuitySource())
                {
                    const scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(pair->liquidIndex_).stkFieldPtr(),
                        node);
                    const scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(pair->vaporIndex_).stkFieldPtr(),
                        node);
                    D = massTransferModel::volumeSource(mdot, rhoL, rhoV);
                }

                rhs[0] += vol * (sign * mdot - phii * rho * D);
                if (D > 0.0)
                {
                    lhs[0] += vol * rho * D;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
