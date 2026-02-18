// File : transitionOnsetReynoldsNumberTransitionSSTNodeTerms.cpp
// Created : Wed Jan 15 2025
// Author : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "transitionOnsetReynoldsNumberTransitionSSTAssembler.h"

namespace accel
{

void transitionOnsetReynoldsNumberTransitionSSTAssembler::
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
    const STKScalarField* ReThetaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar ce2 = model_->ce2();
    scalar cThetat = model_->cThetat();

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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
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
            rhs[0] -= -div * ReTheta;

            // Compute vorticity wij
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute Velocity Mag
            scalar Umag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                Umag += U[i] * U[i];
            }
            Umag = std::max(std::sqrt(std::abs(Umag)), SMALL);

            scalar dUids = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    dUids += dudx[SPATIAL_DIM * i + j] * U[i] * U[j];
                }
            }
            scalar dUdsMag = dUids / std::pow(std::max(Umag, SMALL), 2.0);

            scalar Tfactor = (500.0 * mu) / (rho * Umag * Umag + SMALL);

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1.0 * std::pow((ReOmega / 1e5), 2.0));

            scalar delta =
                (375.0 * wijMag * mu * ReTheta * y) / (rho * Umag * Umag);
            if (delta == 0.0)
                delta = SMALL;

            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar Tu = std::max(
                100.0 * std::sqrt((2.0 * k) / 3.0) / (Umag + SMALL), 0.027);

            // Reynolds Theta Equivalent Iteration
            scalar ReTheta_eq = 0.0;

            const scalar tol = 1e-12;
            const label maxIter = 150;

            scalar lambda = 0.0;
            scalar ReThetatOld, ReThetatNew, ReThetatTol;
            label iteration = 0;

            // Tu is in percent (%), NOT FRACTIONAL!!
            scalar FTu = 0.0;

            if (Tu > 1.3)
            {
                FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
            }
            else
            {
                FTu = 1173.51 - 589.428 * Tu + 0.2196 / std::pow(Tu, 2.0);
            }

            if (lambda > 0.0)
            {
                ReThetatNew =
                    FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                     std::exp(-2.0 * Tu));
            }
            else
            {
                ReThetatNew =
                    FTu *
                    (1.0 + (12.986 * lambda + 123.66 * std::pow(lambda, 2.0) +
                            405.689 * std::pow(lambda, 3.0)) *
                               std::exp(-std::pow(Tu / 1.5, 1.5)));
            }
            ReThetatNew = std::max(ReThetatNew, 20.0);

            ReThetatTol = ReThetatNew * tol;

            do
            {
                ReThetatOld = ReThetatNew;
                lambda = std::max(
                    std::min(std::pow(ReThetatOld, 2.0) * mu * dUdsMag /
                                 (std::pow(Umag, 2.0) + SMALL),
                             0.1),
                    -0.1);

                if (Tu > 1.3)
                {
                    FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
                }
                else
                {
                    FTu = 1173.51 - 589.428 * Tu +
                          0.2196 / (std::pow(Tu, 2.0) + SMALL);
                }

                if (lambda > 0.0)
                {
                    ReThetatNew =
                        FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                         std::exp(-2.0 * Tu));
                }
                else
                {
                    ReThetatNew =
                        FTu * (1.0 + (12.986 * lambda +
                                      123.66 * std::pow(lambda, 2.0) +
                                      405.689 * std::pow(lambda, 3.0)) *
                                         std::exp(-std::pow(Tu / 1.5, 1.5)));
                }
                ReThetatNew = std::max(ReThetatNew, 20.0);

                if (++iteration > maxIter &&
                    std::abs(ReThetatNew - ReThetatOld) > tol)
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << std::abs(ReThetatNew - ReThetatOld)
                              << "\033[0m" << std::endl;
                }
                if (std::isnan(ReThetatNew))
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << ReThetatNew << "\033[0m" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } while (std::abs(ReThetatNew - ReThetatOld) > ReThetatTol);

            ReTheta_eq = ReThetatNew;

            // Calculate Production Term
            const scalar Ptheta = cThetat * (rho / (Tfactor + SMALL)) *
                                  (ReTheta_eq - ReTheta) * (1.0 - fThetat);
            const scalar linearizedProductionFlux =
                (cThetat * (rho / (Tfactor + SMALL)) * (1.0 - fThetat));

            // Adding RHS - LHS
            rhs[0] += (Ptheta)*vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol);

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void transitionOnsetReynoldsNumberTransitionSSTAssembler::
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
    const STKScalarField* ReThetaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar ce2 = model_->ce2();
    scalar cThetat = model_->cThetat();

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
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar ReThetaOld =
                *stk::mesh::field_data(*ReThetaSTKFieldPtrOld, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * ReTheta;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * ReTheta + lhsfacOld * ReThetaOld);

            // Compute Vorticity Mag Wij
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute Velocity Mag
            scalar Umag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                Umag += U[i] * U[i];
            }
            Umag = std::max(std::sqrt(std::abs(Umag)), SMALL);

            scalar dUids = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    dUids += dudx[SPATIAL_DIM * i + j] * U[i] * U[j];
                }
            }
            scalar dUdsMag = dUids / std::pow(std::max(Umag, SMALL), 2.0);

            scalar Tfactor = (500.0 * mu) / (rho * Umag * Umag + SMALL);

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1.0 * std::pow((ReOmega / 1e5), 2.0));

            scalar delta =
                (375.0 * wijMag * mu * ReTheta * y) / (rho * Umag * Umag);
            if (delta == 0.0)
                delta = SMALL;

            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar Tu = std::max(
                100.0 * std::sqrt((2.0 * k) / 3.0) / (Umag + SMALL), 0.027);

            // Reynolds Theta Equivalent Iteration
            scalar ReTheta_eq = 0.0;

            const scalar tol = 1e-12;
            const label maxIter = 150;

            scalar lambda = 0.0;
            scalar ReThetatOld, ReThetatNew, ReThetatTol;
            label iteration = 0;

            // Tu is in percent (%), NOT FRACTIONAL!!
            scalar FTu = 0.0;

            if (Tu > 1.3)
            {
                FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
            }
            else
            {
                FTu = 1173.51 - 589.428 * Tu + 0.2196 / std::pow(Tu, 2.0);
            }

            if (lambda > 0.0)
            {
                ReThetatNew =
                    FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                     std::exp(-2.0 * Tu));
            }
            else
            {
                ReThetatNew =
                    FTu *
                    (1.0 + (12.986 * lambda + 123.66 * std::pow(lambda, 2.0) +
                            405.689 * std::pow(lambda, 3.0)) *
                               std::exp(-std::pow(Tu / 1.5, 1.5)));
            }
            ReThetatNew = std::max(ReThetatNew, 20.0);

            ReThetatTol = ReThetatNew * tol;

            do
            {
                ReThetatOld = ReThetatNew;
                lambda = std::max(
                    std::min(std::pow(ReThetatOld, 2.0) * mu * dUdsMag /
                                 (std::pow(Umag, 2.0) + SMALL),
                             0.1),
                    -0.1);

                if (Tu > 1.3)
                {
                    FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
                }
                else
                {
                    FTu = 1173.51 - 589.428 * Tu +
                          0.2196 / (std::pow(Tu, 2.0) + SMALL);
                }

                if (lambda > 0.0)
                {
                    ReThetatNew =
                        FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                         std::exp(-2.0 * Tu));
                }
                else
                {
                    ReThetatNew =
                        FTu * (1.0 + (12.986 * lambda +
                                      123.66 * std::pow(lambda, 2.0) +
                                      405.689 * std::pow(lambda, 3.0)) *
                                         std::exp(-std::pow(Tu / 1.5, 1.5)));
                }
                ReThetatNew = std::max(ReThetatNew, 20.0);

                if (++iteration > maxIter &&
                    std::abs(ReThetatNew - ReThetatOld) > tol)
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << std::abs(ReThetatNew - ReThetatOld)
                              << "\033[0m" << std::endl;
                }
                if (std::isnan(ReThetatNew))
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << ReThetatNew << "\033[0m" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } while (std::abs(ReThetatNew - ReThetatOld) > ReThetatTol);

            ReTheta_eq = ReThetatNew;

            // Calculate Production Term
            const scalar Ptheta = cThetat * (rho / (Tfactor + SMALL)) *
                                  (ReTheta_eq - ReTheta) * (1.0 - fThetat);
            const scalar linearizedProductionFlux =
                (cThetat * (rho / (Tfactor + SMALL)) * (1.0 - fThetat));

            // Adding RHS - LHS
            rhs[0] += (Ptheta)*vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol);

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * ReTheta * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void transitionOnsetReynoldsNumberTransitionSSTAssembler::
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
    const STKScalarField* ReThetaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar ce2 = model_->ce2();
    scalar cThetat = model_->cThetat();

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
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar ReThetaOld =
                *stk::mesh::field_data(*ReThetaSTKFieldPtrOld, node);
            scalar ReThetaOldOld =
                *stk::mesh::field_data(*ReThetaSTKFieldPtrOldOld, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * ReTheta;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * ReTheta + lhsfacOld * ReThetaOld +
                       lhsfacOldOld * ReThetaOldOld);

            // Compute Vorticity Mag Wij
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            // Compute Velocity Mag
            scalar Umag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                Umag += U[i] * U[i];
            }
            Umag = std::max(std::sqrt(std::abs(Umag)), SMALL);

            scalar dUids = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    dUids += dudx[SPATIAL_DIM * i + j] * U[i] * U[j];
                }
            }
            scalar dUdsMag = dUids / std::pow(std::max(Umag, SMALL), 2.0);

            scalar Tfactor = (500.0 * mu) / (rho * Umag * Umag + SMALL);

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1.0 * std::pow((ReOmega / 1e5), 2.0));

            scalar delta =
                (375.0 * wijMag * mu * ReTheta * y) / (rho * Umag * Umag);
            if (delta == 0.0)
                delta = SMALL;

            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar Tu = std::max(
                100.0 * std::sqrt((2.0 * k) / 3.0) / (Umag + SMALL), 0.027);

            // Reynolds Theta Equivalent Iteration
            scalar ReTheta_eq = 0.0;

            const scalar tol = 1e-12;
            const label maxIter = 150;

            scalar lambda = 0.0;
            scalar ReThetatOld, ReThetatNew, ReThetatTol;
            label iteration = 0;

            // Tu is in percent (%), NOT FRACTIONAL!!
            scalar FTu = 0.0;

            if (Tu > 1.3)
            {
                FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
            }
            else
            {
                FTu = 1173.51 - 589.428 * Tu + 0.2196 / std::pow(Tu, 2.0);
            }

            if (lambda > 0.0)
            {
                ReThetatNew =
                    FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                     std::exp(-2.0 * Tu));
            }
            else
            {
                ReThetatNew =
                    FTu *
                    (1.0 + (12.986 * lambda + 123.66 * std::pow(lambda, 2.0) +
                            405.689 * std::pow(lambda, 3.0)) *
                               std::exp(-std::pow(Tu / 1.5, 1.5)));
            }
            ReThetatNew = std::max(ReThetatNew, 20.0);

            ReThetatTol = ReThetatNew * tol;

            do
            {
                ReThetatOld = ReThetatNew;
                lambda = std::max(
                    std::min(std::pow(ReThetatOld, 2.0) * mu * dUdsMag /
                                 (std::pow(Umag, 2.0) + SMALL),
                             0.1),
                    -0.1);

                if (Tu > 1.3)
                {
                    FTu = 331.5 * std::pow((Tu - 0.5658), -0.671);
                }
                else
                {
                    FTu = 1173.51 - 589.428 * Tu +
                          0.2196 / (std::pow(Tu, 2.0) + SMALL);
                }

                if (lambda > 0.0)
                {
                    ReThetatNew =
                        FTu * (1.0 + 0.275 * (1.0 - std::exp(-35.0 * lambda)) *
                                         std::exp(-2.0 * Tu));
                }
                else
                {
                    ReThetatNew =
                        FTu * (1.0 + (12.986 * lambda +
                                      123.66 * std::pow(lambda, 2.0) +
                                      405.689 * std::pow(lambda, 3.0)) *
                                         std::exp(-std::pow(Tu / 1.5, 1.5)));
                }
                ReThetatNew = std::max(ReThetatNew, 20.0);

                if (++iteration > maxIter &&
                    std::abs(ReThetatNew - ReThetatOld) > tol)
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << std::abs(ReThetatNew - ReThetatOld)
                              << "\033[0m" << std::endl;
                }
                if (std::isnan(ReThetatNew))
                {
                    std::cout << "\033[1;31mERROR! - Reynold Theta Iterations "
                                 "Diverged - err = "
                              << ReThetatNew << "\033[0m" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } while (std::abs(ReThetatNew - ReThetatOld) > ReThetatTol);

            ReTheta_eq = ReThetatNew;

            // Calculate Production Term
            const scalar Ptheta = cThetat * (rho / (Tfactor + SMALL)) *
                                  (ReTheta_eq - ReTheta) * (1.0 - fThetat);
            const scalar linearizedProductionFlux =
                (cThetat * (rho / (Tfactor + SMALL)) * (1.0 - fThetat));

            // Adding RHS - LHS
            rhs[0] += (Ptheta)*vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol);

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * ReTheta * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
