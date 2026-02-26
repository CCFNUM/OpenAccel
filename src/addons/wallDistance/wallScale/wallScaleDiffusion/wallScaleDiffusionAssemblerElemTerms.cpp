// File       : wallScaleDiffusionAssemblerElemTerms.cpp
// Created    : Thu Aug 21 2025 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "wallScaleDiffusionAssembler.h"

namespace accel
{

void wallScaleDiffusionAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // TODO: Account for BLOCKSIZE in
    // space for LHS/RHS; nodesPerElem*nodesPerElem* and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_yScale;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;

    // Get transport fields/side fields
    const auto& yScaleSTKFieldRef = phi_->stkFieldRef();

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
        ws_yScale.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);

        // pointer to lhs/rhs
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_yScale = &ws_yScale[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_shape_function = &ws_shape_function[0];

        // extract shape function
        if (isShifted)
        {
            meSCS->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_shape_function[0]);
        }

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
                p_yScale[ni] = *stk::mesh::field_data(yScaleSTKFieldRef, node);

                // gather vectors
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
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

                const label rowL = il * nodesPerElement;
                const label rowR = ir * nodesPerElement;

                const label rLiL_i = rowL + il;
                const label rLiR_i = rowL + ir;
                const label rRiR_i = rowR + ir;
                const label rRiL_i = rowR + il;

                //================================
                // Diffusion: Positive-Definite
                //================================

                scalar lhs_riC = 0.0;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    scalar lhs_riC_i = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = p_scs_areav[ip * SPATIAL_DIM + j];

                        // -dphi/dxj*A_j; fixed i over j loop; see
                        // below..
                        const scalar lhsfacDiff_i =
                            -p_dndx[offSetDnDx + j] * axj;

                        // accumulate
                        lhs_riC_i += lhsfacDiff_i;
                    }

                    // rhs
                    const scalar yScale = p_yScale[ic];
                    p_rhs[il] -= lhs_riC_i * yScale;
                    p_rhs[ir] += lhs_riC_i * yScale;

                    // accumulate absolute
                    lhs_riC += 0.5 * std::abs(lhs_riC_i);
                }

                // accmulate
                p_lhs[rLiL_i] += lhs_riC;
                p_lhs[rLiR_i] -= lhs_riC;
                p_lhs[rRiR_i] += lhs_riC;
                p_lhs[rRiL_i] -= lhs_riC;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
