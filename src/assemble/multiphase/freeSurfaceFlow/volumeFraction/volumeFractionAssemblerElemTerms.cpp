// File : volumeFractionAssemblerElemTerms.cpp
// Created : Mon Jan 27 2025
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "volumeFractionAssembler.h"

namespace accel
{

void volumeFractionAssembler::assembleElemTermsInterior_(const domain* domain,
                                                         Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const scalar fct =
        domain->multiphase_.freeSurfaceModel_.fluxCorrectedTransport_ ? 0.0
                                                                      : 1.0;

    // TODO: [2024-02-29] Account for BLOCKSIZE in
    // space for LHS/RHS; nodesPerElem*nodesPerElem* and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_beta;
    std::vector<scalar> ws_gradAlpha;
    std::vector<scalar> ws_Gamma;
    std::vector<scalar> ws_magU;
    std::vector<scalar> ws_scv_volume;
    std::vector<scalar> ws_nHat;
    std::vector<scalar> ws_rho;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;
    std::vector<scalar> ws_dndx_scv;
    std::vector<scalar> ws_deriv_scv;
    std::vector<scalar> ws_det_j_scv;

    // ip values
    std::vector<scalar> coordIp(SPATIAL_DIM);

    // pointers
    scalar* p_coordIp = &coordIp[0];

    // Get transport fields/side fields
    const auto& alphaSTKFieldRef = phi_->stkFieldRef();
    const auto& gradAlphaSTKFieldRef = phi_->gradRef().stkFieldRef();
    const auto& blendingFactorSTKFieldRef =
        phi_->blendingFactorRef().stkFieldRef();
    const auto& nHatSTKFieldRef = model_->nHatRef(phaseIndex_).stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const scalar gamma = model_->gamma(domain);

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for field
    const bool isShifted = phi_->isShifted();
    const bool isGradientShifted = phi_->isGradientShifted();

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master elements
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(
                elementBucket.topology());
        MasterElement* meSCV =
            accel::MasterElementRepo::get_volume_master_element(
                elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        const label numScvIp = meSCV->numIntPoints_;
        const label* ipNodeMap = meSCV->ipNodeMap();

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_alpha.resize(nodesPerElement);
        ws_beta.resize(nodesPerElement);
        ws_gradAlpha.resize(nodesPerElement * SPATIAL_DIM);
        ws_Gamma.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);

        ws_magU.resize(nodesPerElement);
        ws_rho.resize(nodesPerElement);
        ws_scv_volume.resize(numScvIp);
        ws_nHat.resize(nodesPerElement * SPATIAL_DIM);
        ws_dndx_scv.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_deriv_scv.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_det_j_scv.resize(numScvIp);

        // pointer to lhs/rhs
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_beta = &ws_beta[0];
        scalar* p_gradAlpha = &ws_gradAlpha[0];
        scalar* p_Gamma = &ws_Gamma[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];
        scalar* p_scv_volume = &ws_scv_volume[0];
        scalar* p_nHat = &ws_nHat[0];
        scalar* p_rho = &ws_rho[0];

        // extract shape function
        if (isShifted)
        {
            meSCS->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meSCS->shape_fcn(&p_coordinate_shape_function[0]);

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            // get elem
            stk::mesh::Entity elem = elementBucket[iElement];

            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(elem);
            label numNodes = bulkData.num_nodes(elem);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // pointers to real data
                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);

                // gather scalars
                p_Gamma[ni] = *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                // gather vectors
                const label niNdim = ni * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_coordinates[niNdim + i] = coords[i];
                }

                // gather 1-dim fields
                const scalar* phi =
                    stk::mesh::field_data(alphaSTKFieldRef, node);
                const scalar* beta =
                    stk::mesh::field_data(blendingFactorSTKFieldRef, node);
                const scalar* gradPhi =
                    stk::mesh::field_data(gradAlphaSTKFieldRef, node);
                const scalar* nHat =
                    stk::mesh::field_data(nHatSTKFieldRef, node);
                const scalar* rho = stk::mesh::field_data(rhoSTKFieldRef, node);

                p_alpha[ni] = phi[0];
                p_beta[ni] = beta[0];
                p_rho[ni] = rho[0];

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // gradient of volume fraction
                    p_gradAlpha[ni * SPATIAL_DIM + j] = gradPhi[j];

                    // interface normal
                    p_nHat[ni * SPATIAL_DIM + j] = nHat[j];
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

            scalar scv_error = 0.0;
            meSCV->determinant(
                1, &p_coordinates[0], &p_scv_volume[0], &scv_error);

            meSCV->grad_op(1,
                           &ws_coordinates[0],
                           &ws_dndx_scv[0],
                           &ws_deriv_scv[0],
                           &ws_det_j_scv[0],
                           &scv_error);

            // compute dndx
            if (isGradientShifted)
            {
                meSCS->shifted_grad_op(1,
                                       &ws_coordinates[0],
                                       &ws_dndx[0],
                                       &ws_deriv[0],
                                       &ws_det_j[0],
                                       &scs_error);
            }
            else
            {
                meSCS->grad_op(1,
                               &ws_coordinates[0],
                               &ws_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scs_error);
            }

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // save off mDot
                const scalar tmDot =
                    (stk::mesh::field_data(*mDotSTKFieldPtr_, elem))[ip];

                // zero out values of interest for this ip
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordIp[j] = 0.0;
                }

                // save off ip values; offset to Shape Function
                scalar GammaIp = 0.0;
                scalar phiIp = 0.0;
                const label offSetSF = ip * nodesPerElement;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];

                    GammaIp += r * p_Gamma[ic];

                    // compute scs ip value
                    phiIp += r * p_alpha[ic];
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordIp[i] +=
                            r_coord * p_coordinates[ic * SPATIAL_DIM + i];
                    }
                }

                //================================
                // Advection
                //================================

                // assemble advection; rhs and upwind contributions
                scalar phiUpwind;
                scalar dcorr = 0;
                if (tmDot > 0)
                {
                    phiUpwind = p_alpha[il];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[il * SPATIAL_DIM + j];
                        dcorr += p_beta[il] * dxj *
                                 p_gradAlpha[il * SPATIAL_DIM + j];
                    }
                }
                else
                {
                    phiUpwind = p_alpha[ir];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[ir * SPATIAL_DIM + j];
                        dcorr += p_beta[ir] * dxj *
                                 p_gradAlpha[ir * SPATIAL_DIM + j];
                    }
                }

                // total upwind advection
                const scalar aflux = tmDot * (phiUpwind + fct * dcorr);

                const label indexL = il;
                const label indexR = ir;

                const label rowL = indexL * nodesPerElement;
                const label rowR = indexR * nodesPerElement;

                const label rLiL_i = rowL + il;
                const label rLiR_i = rowL + ir;
                const label rRiR_i = rowR + ir;
                const label rRiL_i = rowR + il;

                // right hand side; L and R
                p_rhs[indexL] -= aflux;
                p_rhs[indexR] += aflux;

                // upwind advection left node
                const scalar alhsfacL = 0.5 * (tmDot + std::abs(tmDot));
                p_lhs[rLiL_i] += alhsfacL;
                p_lhs[rRiL_i] -= alhsfacL;

                // upwind advection right node
                const scalar alhsfacR = 0.5 * (tmDot - std::abs(tmDot));
                p_lhs[rRiR_i] -= alhsfacR;
                p_lhs[rLiR_i] += alhsfacR;
            }

            //================================
            // Interface Compression
            //================================

            // determine nodal velocity magnitude
            for (label ni = 0; ni < numNodes; ni++)
            {
                stk::mesh::Entity node = nodeRels[ni];

                // pointers to real data
                const scalar* vel = stk::mesh::field_data(USTKFieldRef, node);

                scalar magU = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    magU += vel[i] * vel[i];
                }

                ws_magU[ni] = stk::math::sqrt(magU);
            }

            // start loop over sub-control volumes
            for (label ip = 0; ip < numScvIp; ip++)
            {
                // nearest node to ip
                const label nn = ipNodeMap[ip];
                const scalar scv = p_scv_volume[ip];

                scalar sharpen = 0.0;
                for (label ni = 0; ni < numNodes; ni++)
                {
                    scalar sum = 0.0;
                    scalar lhsOne = 0.0;
                    scalar lhsTwo = 0.0;
                    const scalar nodalVof = p_alpha[ni];
                    const scalar nodalRho = p_rho[ni];
                    const label offsetDnDx =
                        SPATIAL_DIM * numNodes * ip + ni * SPATIAL_DIM;

                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        // compressive velocity
                        const scalar Ucomp = gamma * ws_magU[ni] *
                                             p_nHat[ni * SPATIAL_DIM + j] * fct;

                        lhsOne +=
                            nodalRho * Ucomp * ws_dndx_scv[offsetDnDx + j];
                        lhsTwo -= nodalRho * 2.0 * Ucomp * nodalVof *
                                  ws_dndx_scv[offsetDnDx + j];
                        sum += nodalRho * Ucomp * ws_dndx_scv[offsetDnDx + j];
                    }

                    sharpen += sum * (nodalVof * (1 - nodalVof));
                    p_lhs[nn * numNodes + ni] +=
                        (lhsOne * 1.0 + lhsTwo * 0.0) * scv;
                }

                // assemble rhs
                p_rhs[nn] -= sharpen * scv;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
