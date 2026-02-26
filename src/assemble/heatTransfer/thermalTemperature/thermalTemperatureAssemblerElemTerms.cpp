// File       : thermalTemperatureAssemblerElemTerms.cpp
// Created    : Thu Apr 14 2024 8:36:38 (+0100)
// Author     : Fabian Wermelinger
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef WITH_THERMAL_TEMPERATURE

#include "thermalTemperatureAssembler.h"

namespace accel
{

void thermalTemperatureAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // TODO: Account for BLOCKSIZE in space for LHS/RHS;
    // nodesPerElem*nodesPerElem* and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_h0;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_beta;
    std::vector<scalar> ws_gradT;
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;

    // ip values
    std::vector<scalar> coordIp(SPATIAL_DIM);

    // pointers
    scalar* p_coordIp = &coordIp[0];

    // Get transport fields/side fields
    const auto& gradTSTKFieldRef = phi_->gradRef().stkFieldRef();
    const auto& betaTSTKFieldRef = phi_->blendingFactorRef().stkFieldRef();
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for field?
    const bool isShifted = phi_->isShifted();

    // shifted ip's for gradients?
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

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

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
        ws_T.resize(nodesPerElement);
        ws_beta.resize(nodesPerElement);
        ws_gradT.resize(nodesPerElement * SPATIAL_DIM);
        ws_lambdaEff.resize(nodesPerElement);
        ws_cp.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);

        // pointer to lhs/rhs
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_T = &ws_T[0];
        scalar* p_beta = &ws_beta[0];
        scalar* p_gradT = &ws_gradT[0];
        scalar* p_lambdaEff = &ws_lambdaEff[0];
        scalar* p_cp = &ws_cp[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];

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
                p_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                p_beta[ni] = *stk::mesh::field_data(betaTSTKFieldRef, node);

                // gather N-dim fields
                const scalar* gradT =
                    stk::mesh::field_data(gradTSTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_gradT[ni * SPATIAL_DIM + j] = gradT[j];
                    p_coordinates[ni * SPATIAL_DIM + j] = coords[j];
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

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

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordIp[j] = 0.0;
                }

                // save off ip values; offset to Shape Function
                scalar lambdaEffIp = 0.0;
                scalar cpIp = 0.0;
                scalar h0Ip = 0.0;
                const label offSetSF = ip * nodesPerElement;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];

                    lambdaEffIp += r * p_lambdaEff[ic];
                    cpIp += r * p_cp[ic];

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
                scalar TUpwind;
                scalar dcorr = 0;
                if (tmDot > 0)
                {
                    TUpwind = p_T[il];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[il * SPATIAL_DIM + j];
                        dcorr +=
                            p_beta[il] * dxj * p_gradT[il * SPATIAL_DIM + j];
                    }
                }
                else
                {
                    TUpwind = p_T[ir];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[ir * SPATIAL_DIM + j];
                        dcorr +=
                            p_beta[ir] * dxj * p_gradT[ir * SPATIAL_DIM + j];
                    }
                }

                // total upwind advection
                const scalar aflux = tmDot * cpIp * (TUpwind + dcorr);

                const label rowL = il * nodesPerElement;
                const label rowR = ir * nodesPerElement;

                const label rLiL_i = rowL + il;
                const label rLiR_i = rowL + ir;
                const label rRiR_i = rowR + ir;
                const label rRiL_i = rowR + il;

                // right hand side; L and R
                p_rhs[il] -= aflux;
                p_rhs[ir] += aflux;

                // upwind advection left node
                const scalar alhsfacL = 0.5 * (tmDot + std::abs(tmDot)) * cpIp;
                p_lhs[rLiL_i] += alhsfacL;
                p_lhs[rRiL_i] -= alhsfacL;

                // upwind advection right node
                const scalar alhsfacR = 0.5 * (tmDot - std::abs(tmDot)) * cpIp;
                p_lhs[rRiR_i] -= alhsfacR;
                p_lhs[rLiR_i] += alhsfacR;

                //================================
                // Diffusion
                //================================

                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label rowL = il * nodesPerElement;
                    const label rowR = ir * nodesPerElement;

                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    scalar lhs_riC_i = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = p_scs_areav[ip * SPATIAL_DIM + j];

                        // -Gamma*dphi/dxj*A_j; fixed i over j loop; see
                        // below..
                        const scalar lhsfacDiff_i =
                            -lambdaEffIp * p_dndx[offSetDnDx + j] * axj;

                        // lhs; il then ir
                        lhs_riC_i += lhsfacDiff_i;
                    }

                    // deal with accumulated lhs and flux for
                    // -Gamma*dphi/dxj*Aj
                    p_lhs[rowL + ic] += lhs_riC_i;
                    p_lhs[rowR + ic] -= lhs_riC_i;

                    const scalar T = p_T[ic];
                    p_rhs[il] -= lhs_riC_i * T;
                    p_rhs[ir] += lhs_riC_i * T;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif
