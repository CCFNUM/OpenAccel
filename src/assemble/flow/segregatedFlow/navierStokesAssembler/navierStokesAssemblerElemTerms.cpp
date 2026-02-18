// File : navierStokesAssemblerElemTerms.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "flowModel.h"
#include "navierStokesAssembler.h"

namespace accel
{

void navierStokesAssembler::assembleElemTermsInterior_(const domain* domain,
                                                       Context* ctx)
{
    auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const bool compressible = domain->isMaterialCompressible();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS uu; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_beta;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_dudx;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_muEff;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_velocity_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;

    // ip values
    std::vector<scalar> uIp(SPATIAL_DIM);
    std::vector<scalar> coordIp(SPATIAL_DIM);

    // extrapolated value from the L/R direction
    std::vector<scalar> uIpL(SPATIAL_DIM);
    std::vector<scalar> uIpR(SPATIAL_DIM);

    // pointers for fast access
    scalar* p_uIp = &uIp[0];
    scalar* p_uIpL = &uIpL[0];
    scalar* p_uIpR = &uIpR[0];
    scalar* p_coordIp = &coordIp[0];

    // deal with state
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& rhoSTKFieldRef = model_->rhoRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& gradUSTKFieldRef = model_->URef().gradRef().stkFieldRef();
    const auto& blendingFactorSTKFieldRef =
        model_->URef().blendingFactorRef().stkFieldRef();

    const auto& mDotSTKFieldRef = model_->mDotRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

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
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_beta.resize(nodesPerElement * SPATIAL_DIM);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_dudx.resize(nodesPerElement * SPATIAL_DIM * SPATIAL_DIM);
        ws_rho.resize(nodesPerElement);
        ws_muEff.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_velocity_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_beta = &ws_beta[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_dudx = &ws_dudx[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_velocity_shape_function = &ws_velocity_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];

        // extract shape function
        if (isUShifted)
        {
            meSCS->shifted_shape_fcn(&p_velocity_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_velocity_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meSCS->shape_fcn(&p_coordinate_shape_function[0]);

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // ip data for this element; scs and scv
            const scalar* mDot =
                stk::mesh::field_data(mDotSTKFieldRef, elementBucket, iElement);

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            stk::mesh::Entity const* nodeRels =
                elementBucket.begin_nodes(iElement);
            label numNodes = elementBucket.num_nodes(iElement);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // pointers to real data
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* beta =
                    stk::mesh::field_data(blendingFactorSTKFieldRef, node);

                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const scalar* dudx =
                    stk::mesh::field_data(gradUSTKFieldRef, node);
                const scalar rho = *stk::mesh::field_data(rhoSTKFieldRef, node);
                const scalar muEff =
                    *stk::mesh::field_data(muEffSTKFieldRef, node);

                // gather scalars
                p_rho[ni] = rho;
                p_muEff[ni] = muEff;

                const label niNdim = ni * SPATIAL_DIM;

                // gather vectors/tensors
                const label row_p_dudx = niNdim * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_U[niNdim + i] = U[i];
                    p_beta[niNdim + i] = beta[i];
                    p_coordinates[niNdim + i] = coords[i];

                    const label row_dudx = i * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_dudx[row_p_dudx + row_dudx + j] = dudx[row_dudx + j];
                    }
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

            // compute dndx
            if (isUGradientShifted)
            {
                meSCS->shifted_grad_op(1,
                                       &p_coordinates[0],
                                       &p_dndx[0],
                                       &ws_deriv[0],
                                       &ws_det_j[0],
                                       &scs_error);
            }
            else
            {
                meSCS->grad_op(1,
                               &p_coordinates[0],
                               &p_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scs_error);
            }

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                const label ipNdim = ip * SPATIAL_DIM;

                const label offSetSF = ip * nodesPerElement;

                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // save off mDot
                const scalar tmDot = mDot[ip];

                // save off some offsets
                const label ilNdim = il * SPATIAL_DIM;
                const label irNdim = ir * SPATIAL_DIM;

                // zero out values of interest for this ip
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uIp[j] = 0.0;
                    p_coordIp[j] = 0.0;
                }

                // compute scs point values; offset to Shape Function; sneak in
                // divU
                scalar muEffIp = 0.0;
                scalar divU = 0.0;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_velocity_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];
                    muEffIp += r * p_muEff[ic];

                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const scalar ui = p_U[ic * SPATIAL_DIM + i];
                        p_uIp[i] += r * ui;
                        p_coordIp[i] +=
                            r_coord * p_coordinates[ic * SPATIAL_DIM + i];
                        divU += ui * p_dndx[offSetDnDx + i];
                    }
                }

                // final upwind extrapolation;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_uIpL[i] = p_U[ilNdim + i];
                    p_uIpR[i] = p_U[irNdim + i];
                }

                // assemble advection; rhs and upwind contributions; add divU
                // stress (explicit) if compressible fluid
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    scalar uiUpwind;
                    scalar dcorr = 0;
                    if (tmDot > 0)
                    {
                        uiUpwind = p_uIpL[i];

                        // deferred correction
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar dxj =
                                p_coordIp[j] -
                                p_coordinates[il * SPATIAL_DIM + j];
                            dcorr += p_beta[SPATIAL_DIM * il + i] * dxj *
                                     p_dudx[il * (SPATIAL_DIM * SPATIAL_DIM) +
                                            i * SPATIAL_DIM + j];
                        }
                    }
                    else
                    {
                        uiUpwind = p_uIpR[i];

                        // deferred correction
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar dxj =
                                p_coordIp[j] -
                                p_coordinates[ir * SPATIAL_DIM + j];
                            dcorr += p_beta[SPATIAL_DIM * ir + i] * dxj *
                                     p_dudx[ir * (SPATIAL_DIM * SPATIAL_DIM) +
                                            i * SPATIAL_DIM + j];
                        }
                    }

                    // total upwind advection
                    const scalar aflux = tmDot * (uiUpwind + dcorr);

                    // divU stress term
                    const scalar divUstress = compressible
                                                  ? 2.0 / 3.0 * muEffIp * divU *
                                                        p_scs_areav[ipNdim + i]
                                                  : 0.0;

                    const label indexL = ilNdim + i;
                    const label indexR = irNdim + i;

                    const label rowL = indexL * nodesPerElement * SPATIAL_DIM;
                    const label rowR = indexR * nodesPerElement * SPATIAL_DIM;

                    const label rLiL_i = rowL + ilNdim + i;
                    const label rLiR_i = rowL + irNdim + i;
                    const label rRiR_i = rowR + irNdim + i;
                    const label rRiL_i = rowR + ilNdim + i;

                    // right hand side; L and R
                    p_rhs[indexL] -= aflux + divUstress;
                    p_rhs[indexR] += aflux + divUstress;

                    // upwind advection left node
                    const scalar alhsfacL = 0.5 * (tmDot + std::abs(tmDot));
                    p_lhs[rLiL_i] += alhsfacL;
                    p_lhs[rRiL_i] -= alhsfacL;

                    // upwind advection right node
                    const scalar alhsfacR = 0.5 * (tmDot - std::abs(tmDot));
                    p_lhs[rRiR_i] -= alhsfacR;
                    p_lhs[rLiR_i] += alhsfacR;
                }

                // stress
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label icNdim = ic * SPATIAL_DIM;

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label indexL = ilNdim + i;
                        const label indexR = irNdim + i;

                        const label rowL =
                            indexL * nodesPerElement * SPATIAL_DIM;
                        const label rowR =
                            indexR * nodesPerElement * SPATIAL_DIM;

                        // viscous stress
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip + icNdim;
                        scalar lhs_riC_i = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = p_scs_areav[ipNdim + j];
                            const scalar uj = p_U[icNdim + j];

                            // -mu*dui/dxj*A_j; fixed i over j loop; see below..
                            const scalar lhsfacDiff_i =
                                -muEffIp * p_dndx[offSetDnDx + j] * axj;

                            // lhs; il then ir
                            lhs_riC_i += lhsfacDiff_i;

                            // -mu*duj/dxi*A_j (transpose vel grad term)
                            const scalar lhsfacDiff_j =
                                -muEffIp * p_dndx[offSetDnDx + i] * axj;

                            // lhs; il then ir
                            p_lhs[rowL + icNdim + j] += lhsfacDiff_j;
                            p_lhs[rowR + icNdim + j] -= lhsfacDiff_j;

                            // rhs; il then ir
                            p_rhs[indexL] -= lhsfacDiff_j * uj;
                            p_rhs[indexR] += lhsfacDiff_j * uj;
                        }

                        // deal with accumulated lhs and flux for -mu*dui/dxj*Aj
                        p_lhs[rowL + icNdim + i] += lhs_riC_i;
                        p_lhs[rowR + icNdim + i] -= lhs_riC_i;

                        const scalar ui = p_U[icNdim + i];
                        p_rhs[indexL] -= lhs_riC_i * ui;
                        p_rhs[indexR] += lhs_riC_i * ui;
                    }
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
