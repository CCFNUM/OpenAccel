// File       : volumeFractionAssemblerNodeTerms.cpp
// Created    : Fri Oct 09 2026
// Author     : Mhamad Mahdi Alloush
// Description: Node terms of the volume fraction equation incl. mass transfer
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "volumeFractionAssembler.h"

namespace accel
{

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
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // other
    scalar dt = field_broker_->controlsRef().getPhysicalTimescale();

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

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

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

            // interphase mass transfer: s mdot - alpha rho D, D = mdot
            // (1/rho_v - 1/rho_l), s = +1 vapor, -1 liquid, 0 others
            for (const auto& fpm : domain->fluidPairModels())
            {
                const auto& mt = fpm.massTransfer_;
                if (mt.option_ != massTransferModelOption::none)
                {
                    scalar sign =
                        mt.vaporIndex_ == phaseIndex_
                            ? 1.0
                            : (mt.liquidIndex_ == phaseIndex_ ? -1.0 : 0.0);

                    scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(mt.liquidIndex_).stkFieldPtr(), node);
                    scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(mt.vaporIndex_).stkFieldPtr(), node);
                    scalar mdot = *stk::mesh::field_data(
                        *model_->massTransferRateSTKFieldPtr(fpm), node);
                    scalar dmdotdalpha = *stk::mesh::field_data(
                        *model_->massTransferRateAlphaCoeffSTKFieldPtr(fpm),
                        node);

                    scalar rhoDfac = mt.includeContinuitySource_
                                         ? rho * (1.0 / rhoV - 1.0 / rhoL)
                                         : 0.0;
                    scalar expansion = phii * rhoDfac;

                    rhs[0] += vol * mdot * (sign - expansion);

                    // -dS/dalpha where positive (d mdot/d alpha = s d
                    // mdot/d alpha_v)
                    lhs[0] +=
                        vol *
                        (std::max(mdot * rhoDfac, 0.0) +
                         std::max(-sign * dmdotdalpha * (sign - expansion),
                                  0.0));
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

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

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

            // interphase mass transfer: s mdot - alpha rho D, D = mdot
            // (1/rho_v - 1/rho_l), s = +1 vapor, -1 liquid, 0 others
            for (const auto& fpm : domain->fluidPairModels())
            {
                const auto& mt = fpm.massTransfer_;
                if (mt.option_ != massTransferModelOption::none)
                {
                    scalar sign =
                        mt.vaporIndex_ == phaseIndex_
                            ? 1.0
                            : (mt.liquidIndex_ == phaseIndex_ ? -1.0 : 0.0);

                    scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(mt.liquidIndex_).stkFieldPtr(), node);
                    scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(mt.vaporIndex_).stkFieldPtr(), node);
                    scalar mdot = *stk::mesh::field_data(
                        *model_->massTransferRateSTKFieldPtr(fpm), node);
                    scalar dmdotdalpha = *stk::mesh::field_data(
                        *model_->massTransferRateAlphaCoeffSTKFieldPtr(fpm),
                        node);

                    scalar rhoDfac = mt.includeContinuitySource_
                                         ? rho * (1.0 / rhoV - 1.0 / rhoL)
                                         : 0.0;
                    scalar expansion = phii * rhoDfac;

                    rhs[0] += vol * mdot * (sign - expansion);

                    // -dS/dalpha where positive (d mdot/d alpha = s d
                    // mdot/d alpha_v)
                    lhs[0] +=
                        vol *
                        (std::max(mdot * rhoDfac, 0.0) +
                         std::max(-sign * dmdotdalpha * (sign - expansion),
                                  0.0));
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

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

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

            // interphase mass transfer: s mdot - alpha rho D, D = mdot
            // (1/rho_v - 1/rho_l), s = +1 vapor, -1 liquid, 0 others
            for (const auto& fpm : domain->fluidPairModels())
            {
                const auto& mt = fpm.massTransfer_;
                if (mt.option_ != massTransferModelOption::none)
                {
                    scalar sign =
                        mt.vaporIndex_ == phaseIndex_
                            ? 1.0
                            : (mt.liquidIndex_ == phaseIndex_ ? -1.0 : 0.0);

                    scalar rhoL = *stk::mesh::field_data(
                        *model_->rhoRef(mt.liquidIndex_).stkFieldPtr(), node);
                    scalar rhoV = *stk::mesh::field_data(
                        *model_->rhoRef(mt.vaporIndex_).stkFieldPtr(), node);
                    scalar mdot = *stk::mesh::field_data(
                        *model_->massTransferRateSTKFieldPtr(fpm), node);
                    scalar dmdotdalpha = *stk::mesh::field_data(
                        *model_->massTransferRateAlphaCoeffSTKFieldPtr(fpm),
                        node);

                    scalar rhoDfac = mt.includeContinuitySource_
                                         ? rho * (1.0 / rhoV - 1.0 / rhoL)
                                         : 0.0;
                    scalar expansion = phii * rhoDfac;

                    rhs[0] += vol * mdot * (sign - expansion);

                    // -dS/dalpha where positive (d mdot/d alpha = s d
                    // mdot/d alpha_v)
                    lhs[0] +=
                        vol *
                        (std::max(mdot * rhoDfac, 0.0) +
                         std::max(-sign * dmdotdalpha * (sign - expansion),
                                  0.0));
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
