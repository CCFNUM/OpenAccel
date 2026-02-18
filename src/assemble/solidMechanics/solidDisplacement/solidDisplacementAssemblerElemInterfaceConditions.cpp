// File       : solidDisplacementAssemblerElemInterfaceConditions.cpp
// Created    : Thu Dec 04 2025 12:30:00 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Nonconformal interface treatment for solid displacement
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

#include "dgInfo.h"
#include "interface.h"
#include "solidDisplacementAssembler.h"

namespace accel
{

void solidDisplacementAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
    {
        const auto& mesh = field_broker_->meshRef();
        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // space for LHS/RHS; nodesPerElement*nodesPerElement and
        // nodesPerElement
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // master element
        std::vector<scalar> ws_shape_function;

        // Get transport fields/side fields
        const auto& sidePhiFluxSTKFieldRef =
            phi_->sideFluxFieldRef().stkFieldRef();

        // Get geometric fields
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

        // shifted ip's for field
        const bool isShifted = phi_->isShifted();

        // shifted ip's for gradients?
        const bool isGradientShifted = phi_->isGradientShifted();

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);

        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // face master element
            MasterElement* meFC = MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // mapping from ip to nodes for this ordinal; face perspective (use
            // with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            // resize some things; matrix related
            const label lhsSize =
                nodesPerSide * SPATIAL_DIM * nodesPerSide * SPATIAL_DIM;
            const label rhsSize = nodesPerSide * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerSide);

            // algorithm related; element
            ws_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_shape_function = &ws_shape_function[0];

            // shape functions
            if (isShifted)
            {
                meFC->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nBoundarySides =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundarySides;
                 ++iSide)
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

                // get side
                stk::mesh::Entity side = sideBucket[iSide];

                // face data
                const scalar* areaVec = stk::mesh::field_data(
                    exposedAreaVecSTKFieldRef, sideBucket, iSide);
                const scalar* phiFluxVec = stk::mesh::field_data(
                    sidePhiFluxSTKFieldRef, sideBucket, iSide);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);

                // fill connected nodes
                for (label iNode = 0; iNode < numSideNodes; ++iNode)
                {
                    stk::mesh::Entity node = sideNodeRels[iNode];

                    // set connected nodes
                    connectedNodes[iNode] = node;
                }

                // loop over side ip's
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = faceIpNodeMap[ip];

                    scalar asq = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        asq += axj * axj;
                    }
                    const scalar amag = std::sqrt(asq);

                    //================================
                    // Diffusion: only diffusion
                    //================================

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        // matrix entries
                        label indexR = nearestNode * SPATIAL_DIM + i;

                        // flux value is stored into domain: multiply by -1
                        p_rhs[indexR] -=
                            (-phiFluxVec[SPATIAL_DIM * ip + i]) * amag;
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
        scalar penaltyFactor =
            interfaceSideInfoPtr->interfPtr()->penaltyFactor();

        const auto& mesh = field_broker_->meshRef();

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        // space for LHS/RHS
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values; both boundary and opposing surface
        std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
        std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);
        std::vector<scalar> cNx(SPATIAL_DIM);
        std::vector<scalar> oNx(SPATIAL_DIM);

        // mapping for -1:1 -> -0.5:0.5 volume element
        std::vector<scalar> currentElementIsoParCoords(SPATIAL_DIM);
        std::vector<scalar> opposingElementIsoParCoordsTrans(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_cNx = &cNx[0];
        scalar* p_oNx = &oNx[0];

        // nodal fields to gather
        std::vector<scalar> ws_c_face_disp;
        std::vector<scalar> ws_o_face_disp;
        std::vector<scalar> ws_c_elem_disp;
        std::vector<scalar> ws_o_elem_disp;
        std::vector<scalar> ws_c_elem_coordinates;
        std::vector<scalar> ws_o_elem_coordinates;
        std::vector<scalar> ws_c_mu;
        std::vector<scalar> ws_o_mu;
        std::vector<scalar> ws_c_lambda;
        std::vector<scalar> ws_o_lambda;

        // master element data
        std::vector<scalar> ws_c_dndx;
        std::vector<scalar> ws_o_dndx;
        std::vector<scalar> ws_c_det_j;
        std::vector<scalar> ws_o_det_j;
        std::vector<scalar> ws_c_general_shape_function;
        std::vector<scalar> ws_o_general_shape_function;

        // interpolated values at integration points
        scalar currentDispBip[SPATIAL_DIM];
        scalar opposingDispBip[SPATIAL_DIM];
        scalar currentMuBip;
        scalar opposingMuBip;
        scalar currentLambdaBip;
        scalar opposingLambdaBip;
        scalar currentInverseLength;
        scalar opposingInverseLength;
        scalar currentStiffness;
        scalar opposingStiffness;
        scalar penaltyIp;
        scalar scs_error;
        scalar c_amag;

        // Get transport fields
        const auto& dispSTKFieldRef = phi_->stkFieldRef();
        const auto& ESTKFieldRef = model_->ERef().stkFieldRef();
        const auto& nuSTKFieldRef = model_->nuRef().stkFieldRef();

        // Get geometric fields
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // Check if plane stress or plane strain
        const bool planeStress = domain->solidMechanics_.planeStress_;

        // rotation matrix (in case of rotational periodicity)
        const auto& rotMat = interfaceSideInfoPtr->rotationMatrix_;

        // extract vector of dgInfo
        const std::vector<std::vector<dgInfo*>>& dgInfoVec =
            interfaceSideInfoPtr->dgInfoVec_;

        for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size());
             iSide++)
        {
            const std::vector<dgInfo*>& faceDgInfoVec = dgInfoVec[iSide];

            // now loop over all the DgInfo objects on this particular exposed
            // face
            for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
            {
                dgInfo* dgInfo = faceDgInfoVec[k];

                if (dgInfo->gaussPointExposed_)
                    continue;

                // extract current/opposing face/element
                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;
                stk::mesh::Entity currentElement = dgInfo->currentElement_;
                stk::mesh::Entity opposingElement = dgInfo->opposingElement_;
                const label currentFaceOrdinal = dgInfo->currentFaceOrdinal_;
                const label opposingFaceOrdinal = dgInfo->opposingFaceOrdinal_;

                // master element; face and volume
                MasterElement* meFCCurrent = dgInfo->meFCCurrent_;
                MasterElement* meFCOpposing = dgInfo->meFCOpposing_;
                MasterElement* meSCSCurrent = dgInfo->meSCSCurrent_;
                MasterElement* meSCSOpposing = dgInfo->meSCSOpposing_;

                // local ip, ordinals, etc
                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                currentIsoParCoords = dgInfo->currentIsoParCoords_;
                opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

                // mapping from ip to nodes for this ordinal
                const label* ipNodeMap =
                    meSCSCurrent->ipNodeMap(currentFaceOrdinal);

                // extract some master element info
                const label currentNodesPerFace = meFCCurrent->nodesPerElement_;
                const label opposingNodesPerFace =
                    meFCOpposing->nodesPerElement_;
                const label currentNodesPerElement =
                    meSCSCurrent->nodesPerElement_;
                const label opposingNodesPerElement =
                    meSCSOpposing->nodesPerElement_;

                // resize some things; matrix related
                const label totalNodes =
                    currentNodesPerElement + opposingNodesPerElement;
                const label lhsSize =
                    totalNodes * SPATIAL_DIM * totalNodes * SPATIAL_DIM;
                const label rhsSize = totalNodes * SPATIAL_DIM;
                lhs.resize(lhsSize);
                rhs.resize(rhsSize);
                scratchIds.resize(rhsSize);
                scratchVals.resize(rhsSize);
                connectedNodes.resize(totalNodes);

                // algorithm related; element; dndx will be at a single gauss
                // point
                ws_c_elem_disp.resize(currentNodesPerElement * SPATIAL_DIM);
                ws_o_elem_disp.resize(opposingNodesPerElement * SPATIAL_DIM);
                ws_c_elem_coordinates.resize(currentNodesPerElement *
                                             SPATIAL_DIM);
                ws_o_elem_coordinates.resize(opposingNodesPerElement *
                                             SPATIAL_DIM);
                ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
                ws_c_det_j.resize(1);
                ws_o_det_j.resize(1);

                // algorithm related; face
                ws_c_face_disp.resize(currentNodesPerFace * SPATIAL_DIM);
                ws_o_face_disp.resize(opposingNodesPerFace * SPATIAL_DIM);
                ws_c_mu.resize(currentNodesPerFace);
                ws_o_mu.resize(opposingNodesPerFace);
                ws_c_lambda.resize(currentNodesPerFace);
                ws_o_lambda.resize(opposingNodesPerFace);
                ws_c_general_shape_function.resize(currentNodesPerFace);
                ws_o_general_shape_function.resize(opposingNodesPerFace);

                // pointers
                scalar* p_lhs = &lhs[0];
                scalar* p_rhs = &rhs[0];

                scalar* p_c_face_disp = &ws_c_face_disp[0];
                scalar* p_o_face_disp = &ws_o_face_disp[0];
                scalar* p_c_elem_disp = &ws_c_elem_disp[0];
                scalar* p_o_elem_disp = &ws_o_elem_disp[0];
                scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];
                scalar* p_c_mu = &ws_c_mu[0];
                scalar* p_o_mu = &ws_o_mu[0];
                scalar* p_c_lambda = &ws_c_lambda[0];
                scalar* p_o_lambda = &ws_o_lambda[0];

                // me pointers
                scalar* p_c_general_shape_function =
                    &ws_c_general_shape_function[0];
                scalar* p_o_general_shape_function =
                    &ws_o_general_shape_function[0];
                scalar* p_c_dndx = &ws_c_dndx[0];
                scalar* p_o_dndx = &ws_o_dndx[0];

                // populate current face_node_ordinals
                const label* c_face_node_ordinals =
                    meSCSCurrent->side_node_ordinals(currentFaceOrdinal);

                // gather current face data
                stk::mesh::Entity const* current_face_node_rels =
                    bulkData.begin_nodes(currentFace);
                const label current_num_face_nodes =
                    bulkData.num_nodes(currentFace);
                for (label ni = 0; ni < current_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = current_face_node_rels[ni];

                    // gather material properties
                    const scalar E = *stk::mesh::field_data(ESTKFieldRef, node);
                    const scalar nu =
                        *stk::mesh::field_data(nuSTKFieldRef, node);

                    // Compute Lame parameters
                    p_c_mu[ni] = E / (2.0 * (1.0 + nu));
                    if (planeStress)
                    {
                        p_c_lambda[ni] = nu * E / ((1.0 + nu) * (1.0 - nu));
                    }
                    else
                    {
                        p_c_lambda[ni] =
                            nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    }

                    // gather displacement vector
                    const scalar* disp =
                        stk::mesh::field_data(dispSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_face_disp[ni * SPATIAL_DIM + i] = disp[i];
                    }
                }

                // populate opposing face_node_ordinals
                const label* o_face_node_ordinals =
                    meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

                // gather opposing face data
                stk::mesh::Entity const* opposing_face_node_rels =
                    bulkData.begin_nodes(opposingFace);
                const label opposing_num_face_nodes =
                    bulkData.num_nodes(opposingFace);
                for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = opposing_face_node_rels[ni];

                    // gather material properties
                    const scalar E = *stk::mesh::field_data(ESTKFieldRef, node);
                    const scalar nu =
                        *stk::mesh::field_data(nuSTKFieldRef, node);

                    // Compute Lame parameters
                    p_o_mu[ni] = E / (2.0 * (1.0 + nu));
                    if (planeStress)
                    {
                        p_o_lambda[ni] = nu * E / ((1.0 + nu) * (1.0 - nu));
                    }
                    else
                    {
                        p_o_lambda[ni] =
                            nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    }

                    // gather displacement vector
                    const scalar* disp =
                        stk::mesh::field_data(dispSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_face_disp[ni * SPATIAL_DIM + i] = disp[i];
                    }
                }

                // gather current element data
                stk::mesh::Entity const* current_elem_node_rels =
                    bulkData.begin_nodes(currentElement);
                const label current_num_elem_nodes =
                    bulkData.num_nodes(currentElement);
                for (label ni = 0; ni < current_num_elem_nodes; ++ni)
                {
                    stk::mesh::Entity node = current_elem_node_rels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather displacement
                    const scalar* disp =
                        stk::mesh::field_data(dispSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_elem_disp[ni * SPATIAL_DIM + i] = disp[i];
                    }

                    // gather coordinates
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }
                }

                // gather opposing element data
                stk::mesh::Entity const* opposing_elem_node_rels =
                    bulkData.begin_nodes(opposingElement);
                const label opposing_num_elem_nodes =
                    bulkData.num_nodes(opposingElement);
                for (label ni = 0; ni < opposing_num_elem_nodes; ++ni)
                {
                    stk::mesh::Entity node = opposing_elem_node_rels[ni];

                    // set connected nodes
                    connectedNodes[ni + current_num_elem_nodes] = node;

                    // gather displacement
                    const scalar* disp =
                        stk::mesh::field_data(dispSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_disp[ni * SPATIAL_DIM + i] = disp[i];
                    }

                    // gather coordinates
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }
                }

                // apply transformations for rotational periodicity
                interfaceSideInfoPtr->rotateVectorListCompact<SPATIAL_DIM>(
                    ws_o_face_disp, opposing_num_face_nodes);
                interfaceSideInfoPtr->rotateVectorList<SPATIAL_DIM>(
                    ws_o_elem_disp, opposing_num_elem_nodes);
                interfaceSideInfoPtr->transformCoordinateList(
                    ws_o_elem_coordinates, opposing_num_elem_nodes);

                // pointer to face data
                const scalar* c_areaVec = stk::mesh::field_data(
                    exposedAreaVecSTKFieldRef, currentFace);

                c_amag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar c_axj =
                        c_areaVec[currentGaussPointId * SPATIAL_DIM + j];
                    c_amag += c_axj * c_axj;
                }
                c_amag = std::sqrt(c_amag);

                // compute normal
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_cNx[i] =
                        c_areaVec[currentGaussPointId * SPATIAL_DIM + i] /
                        c_amag;
                }

                // compute opposing normal
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_oNx[i] = -p_cNx[i];
                }

                // project from side to element
                meSCSCurrent->sidePcoords_to_elemPcoords(
                    currentFaceOrdinal,
                    1,
                    &currentIsoParCoords[0],
                    &currentElementIsoParCoords[0]);
                meSCSOpposing->sidePcoords_to_elemPcoords(
                    opposingFaceOrdinal,
                    1,
                    &opposingIsoParCoords[0],
                    &opposingElementIsoParCoordsTrans[0]);

                // compute dndx
                scs_error = 0.0;
                meSCSCurrent->general_face_grad_op(
                    currentFaceOrdinal,
                    &currentElementIsoParCoords[0],
                    &p_c_elem_coordinates[0],
                    &p_c_dndx[0],
                    &ws_c_det_j[0],
                    &scs_error);
                meSCSOpposing->general_face_grad_op(
                    opposingFaceOrdinal,
                    &opposingElementIsoParCoordsTrans[0],
                    &p_o_elem_coordinates[0],
                    &p_o_dndx[0],
                    &ws_o_det_j[0],
                    &scs_error);

                // current inverse length scale
                currentInverseLength = 0.0;
                for (label ic = 0; ic < current_num_face_nodes; ++ic)
                {
                    const label faceNodeNumber = c_face_node_ordinals[ic];
                    const label offSetDnDx = faceNodeNumber * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];
                        currentInverseLength += dndxj * nxj;
                    }
                }

                // opposing inverse length scale
                opposingInverseLength = 0.0;
                for (label ic = 0; ic < opposing_num_face_nodes; ++ic)
                {
                    const label faceNodeNumber = o_face_node_ordinals[ic];
                    const label offSetDnDx = faceNodeNumber * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_oNx[j];
                        const scalar dndxj = p_o_dndx[offSetDnDx + j];
                        opposingInverseLength += dndxj * nxj;
                    }
                }

                // interpolate face data at integration point
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    currentDispBip[i] = 0.0;
                    opposingDispBip[i] = 0.0;
                }

                meFCCurrent->interpolatePoint(SPATIAL_DIM,
                                              &currentIsoParCoords[0],
                                              &p_c_face_disp[0],
                                              &currentDispBip[0]);
                meFCOpposing->interpolatePoint(SPATIAL_DIM,
                                               &opposingIsoParCoords[0],
                                               &p_o_face_disp[0],
                                               &opposingDispBip[0]);

                currentMuBip = 0.0;
                meFCCurrent->interpolatePoint(
                    1, &currentIsoParCoords[0], &ws_c_mu[0], &currentMuBip);

                opposingMuBip = 0.0;
                meFCOpposing->interpolatePoint(
                    1, &opposingIsoParCoords[0], &ws_o_mu[0], &opposingMuBip);

                currentLambdaBip = 0.0;
                meFCCurrent->interpolatePoint(1,
                                              &currentIsoParCoords[0],
                                              &ws_c_lambda[0],
                                              &currentLambdaBip);

                opposingLambdaBip = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_o_lambda[0],
                                               &opposingLambdaBip);

                // compute average stiffness for penalty (using bulk modulus)
                currentStiffness = currentLambdaBip + 2.0 * currentMuBip / 3.0;
                opposingStiffness =
                    opposingLambdaBip + 2.0 * opposingMuBip / 3.0;

                // penalty computation (similar to diffusion penalty)
                penaltyIp = penaltyFactor *
                            (currentStiffness * currentInverseLength +
                             opposingStiffness * opposingInverseLength) /
                            2.0;

                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                    p_lhs[p] = 0.0;
                for (label p = 0; p < rhsSize; ++p)
                    p_rhs[p] = 0.0;

                // extract nearest node
                const label nn = ipNodeMap[currentGaussPointId];

                // compute general shape function at current and opposing
                // integration points
                meFCCurrent->general_shape_fcn(1,
                                               &currentIsoParCoords[0],
                                               &ws_c_general_shape_function[0]);
                meFCOpposing->general_shape_fcn(
                    1,
                    &opposingIsoParCoords[0],
                    &ws_o_general_shape_function[0]);

                // assemble for each displacement component
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    // penalty term: enforces displacement continuity
                    const scalar ncPenalty =
                        penaltyIp * (currentDispBip[i] - opposingDispBip[i]);

                    // assemble residual
                    const label indexR = nn * SPATIAL_DIM + i;
                    p_rhs[indexR] -= ncPenalty * c_amag;

                    // set-up row for matrix
                    const label rowR = indexR * totalNodes * SPATIAL_DIM;

                    // sensitivities; current face (penalty)
                    const scalar lhsFacC = penaltyIp * c_amag;
                    for (label ic = 0; ic < currentNodesPerFace; ++ic)
                    {
                        const label icNdim =
                            c_face_node_ordinals[ic] * SPATIAL_DIM + i;
                        const scalar r = p_c_general_shape_function[ic];
                        p_lhs[rowR + icNdim] += r * lhsFacC;
                    }

                    // sensitivities; opposing face (penalty)
                    const scalar lhsFacO = penaltyIp * c_amag;
                    for (label ic = 0; ic < opposingNodesPerFace; ++ic)
                    {
                        const label icNdim = (o_face_node_ordinals[ic] +
                                              currentNodesPerElement) *
                                                 SPATIAL_DIM +
                                             i;
                        const scalar r = p_o_general_shape_function[ic];
                        p_lhs[rowR + icNdim] -= r * lhsFacO;
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} // namespace accel

#endif /* HAS_INTERFACE */
