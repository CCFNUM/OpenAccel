// File       : wallScaleDiffusionAssemblerElemInterfaceConditions.cpp
// Created    : Thu Jun 19 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

#include "wallScaleDiffusionAssembler.h"

namespace accel
{

void wallScaleDiffusionAssembler::assembleElemTermsInterfaces_(
    const domain* domain,
    Context* ctx)
{
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isConformalTreatment())
            continue;

        if (interf->isInternal())
        {
            assembleElemTermsInterfaceSide_(
                domain, interf->masterInfoPtr(), ctx);
            assembleElemTermsInterfaceSide_(
                domain, interf->slaveInfoPtr(), ctx);
        }
        else
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(domain->index());

            if (interf->isFluidSolidType())
            {
                // treat as a no-slip wall
                assembleElemTermsInterfaceSideNoSlipWall_(
                    domain, interfaceSideInfoPtr, ctx);
            }
            else
            {
                assembleElemTermsInterfaceSide_(
                    domain, interfaceSideInfoPtr, ctx);
            }
        }
    }
}

void wallScaleDiffusionAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    scalar penaltyFactor = interfaceSideInfoPtr->interfPtr()->penaltyFactor();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
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
    std::vector<scalar> opposingElementIsoParCoords(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_cNx = &cNx[0];
    scalar* p_oNx = &oNx[0];

    // nodal fields to gather
    std::vector<scalar> ws_c_face_yScale;
    std::vector<scalar> ws_o_face_yScale;
    std::vector<scalar> ws_c_elem_yScale;
    std::vector<scalar> ws_o_elem_yScale;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // Get transport fields/side fields
    const auto& yScaleSTKFieldRef = phi_->stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // rotation matrix (in case of rotational periodicity)
    const utils::matrix& rotMat = interfaceSideInfoPtr->rotationMatrix_;

    // extract vector of dgInfo
    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
        interfaceSideInfoPtr->dgInfoVec_;

    for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size()); iSide++)
    {
        const std::vector<dgInfo*>& faceDgInfoVec = dgInfoVec[iSide];

        // now loop over all the DgInfo objects on this
        // particular exposed face
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
            const label opposingNodesPerFace = meFCOpposing->nodesPerElement_;
            const label currentNodesPerElement = meSCSCurrent->nodesPerElement_;
            const label opposingNodesPerElement =
                meSCSOpposing->nodesPerElement_;

            // resize some things; matrix related
            const label totalNodes =
                currentNodesPerElement + opposingNodesPerElement;
            const label lhsSize = totalNodes * totalNodes;
            const label rhsSize = totalNodes;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(totalNodes);

            // algorithm related; element; dndx will be at a
            // single gauss point...
            ws_c_elem_yScale.resize(currentNodesPerElement);
            ws_o_elem_yScale.resize(opposingNodesPerElement);
            ws_c_elem_coordinates.resize(currentNodesPerElement * SPATIAL_DIM);
            ws_o_elem_coordinates.resize(opposingNodesPerElement * SPATIAL_DIM);
            ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
            ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
            ws_c_det_j.resize(1);
            ws_o_det_j.resize(1);

            // algorithm related; face
            ws_c_face_yScale.resize(currentNodesPerFace);
            ws_o_face_yScale.resize(opposingNodesPerFace);
            ws_c_general_shape_function.resize(currentNodesPerFace);
            ws_o_general_shape_function.resize(opposingNodesPerFace);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];

            scalar* p_c_face_yScale = &ws_c_face_yScale[0];
            scalar* p_o_face_yScale = &ws_o_face_yScale[0];
            scalar* p_c_elem_yScale = &ws_c_elem_yScale[0];
            scalar* p_o_elem_yScale = &ws_o_elem_yScale[0];
            scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
            scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];

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

                // gather; scalar
                p_c_face_yScale[ni] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);
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

                // gather; scalar
                p_o_face_yScale[ni] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);
            }

            // gather current element data; sneak in first of
            // connected nodes
            stk::mesh::Entity const* current_elem_node_rels =
                bulkData.begin_nodes(currentElement);
            const label current_num_elem_nodes =
                bulkData.num_nodes(currentElement);
            for (label ni = 0; ni < current_num_elem_nodes; ++ni)
            {
                stk::mesh::Entity node = current_elem_node_rels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather; scalar
                p_c_elem_yScale[ni] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_c_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                }
            }

            // gather opposing element data; sneak in second
            // connected nodes
            stk::mesh::Entity const* opposing_elem_node_rels =
                bulkData.begin_nodes(opposingElement);
            const label opposing_num_elem_nodes =
                bulkData.num_nodes(opposingElement);
            for (label ni = 0; ni < opposing_num_elem_nodes; ++ni)
            {
                stk::mesh::Entity node = opposing_elem_node_rels[ni];

                // set connected nodes
                connectedNodes[ni + current_num_elem_nodes] = node;

                // gather; scalar
                p_o_elem_yScale[ni] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_o_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                }
            }

            // apply transformations
            interfaceSideInfoPtr->rotateVectorListCompact<1>(
                ws_o_face_yScale, opposing_num_face_nodes);
            interfaceSideInfoPtr->rotateVectorList<1>(ws_o_elem_yScale,
                                                      opposing_num_elem_nodes);
            interfaceSideInfoPtr->transformCoordinateList(
                ws_o_elem_coordinates, opposing_num_elem_nodes);

            // pointer to face data
            const scalar* c_areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, currentFace);

            scalar c_amag = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar c_axj =
                    c_areaVec[currentGaussPointId * SPATIAL_DIM + j];
                c_amag += c_axj * c_axj;
            }
            c_amag = std::sqrt(c_amag);

            // now compute normal
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                p_cNx[i] =
                    c_areaVec[currentGaussPointId * SPATIAL_DIM + i] / c_amag;
            }

            // compute opposing normal: in theory it is assumed
            // that the current and opposing sub-control surfaces
            // are sufficiently planar
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                p_oNx[i] = -p_cNx[i];
            }

            // project from side to element; method deals with
            // the -1:1 isInElement range to the proper
            // underlying CVFEM range
            meSCSCurrent->sidePcoords_to_elemPcoords(
                currentFaceOrdinal,
                1,
                &currentIsoParCoords[0],
                &currentElementIsoParCoords[0]);
            meSCSOpposing->sidePcoords_to_elemPcoords(
                opposingFaceOrdinal,
                1,
                &opposingIsoParCoords[0],
                &opposingElementIsoParCoords[0]);

            // compute dndx
            scalar scs_error = 0.0;
            meSCSCurrent->general_face_grad_op(currentFaceOrdinal,
                                               &currentElementIsoParCoords[0],
                                               &p_c_elem_coordinates[0],
                                               &p_c_dndx[0],
                                               &ws_c_det_j[0],
                                               &scs_error);
            meSCSOpposing->general_face_grad_op(opposingFaceOrdinal,
                                                &opposingElementIsoParCoords[0],
                                                &p_o_elem_coordinates[0],
                                                &p_o_dndx[0],
                                                &ws_o_det_j[0],
                                                &scs_error);

            // current inverse length scale; can loop over face
            // nodes to avoid "nodesOnFace" array
            scalar currentInverseLength = 0.0;
            for (label ic = 0; ic < current_num_face_nodes; ++ic)
            {
                const label faceNodeNumber = c_face_node_ordinals[ic];
                const label offSetDnDx =
                    faceNodeNumber * SPATIAL_DIM; // single intg. point
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_cNx[j];
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];
                    currentInverseLength += dndxj * nxj;
                }
            }

            // opposing inverse length scale; can loop over face
            // nodes to avoid "nodesOnFace" array
            scalar opposingInverseLength = 0.0;
            for (label ic = 0; ic < opposing_num_face_nodes; ++ic)
            {
                const label faceNodeNumber = o_face_node_ordinals[ic];
                const label offSetDnDx =
                    faceNodeNumber * SPATIAL_DIM; // single intg. point
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_oNx[j];
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];
                    opposingInverseLength += dndxj * nxj;
                }
            }

            // interpolate face data; current and opposing...
            scalar currentYScaleBip = 0.0;
            meFCCurrent->interpolatePoint(1,
                                          &currentIsoParCoords[0],
                                          &ws_c_face_yScale[0],
                                          &currentYScaleBip);

            scalar opposingYScaleBip = 0.0;
            meFCOpposing->interpolatePoint(1,
                                           &opposingIsoParCoords[0],
                                           &ws_o_face_yScale[0],
                                           &opposingYScaleBip);

            // compute diffusion vector; current
            scalar currentDiffFluxBip = 0;
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_cNx[j];
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];

                    const scalar yScale = p_c_elem_yScale[ic];

                    // -dphi/dxj*Aj
                    currentDiffFluxBip += -dndxj * nxj * yScale;
                }
            }

            // compute diffusion vector; opposing
            scalar opposingDiffFluxBip = 0;
            for (label ic = 0; ic < opposingNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_oNx[j];
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];

                    const scalar yScale = p_o_elem_yScale[ic];

                    // -dphi/dxj*Aj
                    opposingDiffFluxBip += -dndxj * nxj * yScale;
                }
            }

            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // extract nearset node
            const label nn = ipNodeMap[currentGaussPointId];

            // compute general shape function at current and
            // opposing integration points
            meFCCurrent->general_shape_fcn(
                1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
            meFCOpposing->general_shape_fcn(
                1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);

            // compute penalty
            const scalar penaltyIp =
                penaltyFactor * (currentInverseLength + opposingInverseLength) /
                2.0;

            // non conformal diffusive flux
            const scalar ncDiffFlux =
                (currentDiffFluxBip - opposingDiffFluxBip) / 2.0;

            // assemble residual; form proper rhs index for
            // current face assembly
            const label indexR = nn;
            p_rhs[indexR] -= ((ncDiffFlux + penaltyIp * (currentYScaleBip -
                                                         opposingYScaleBip)) *
                              c_amag);

            // set-up row for matrix
            const label rowR = indexR * totalNodes;

            // sensitivities; current face (penalty); use general shape function
            // for this single ip
            const scalar lhsFacC = penaltyIp * c_amag;
            for (label ic = 0; ic < currentNodesPerFace; ++ic)
            {
                const label icNdim = c_face_node_ordinals[ic];
                const scalar r = p_c_general_shape_function[ic];
                p_lhs[rowR + icNdim] += r * lhsFacC;
            }

            // sensitivities; current element (diffusion)
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const label icNdim = ic;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_cNx[j];
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];

                    // -dphi/dxj*nj*dS
                    p_lhs[rowR + icNdim] += -dndxj * nxj * c_amag / 2.0;
                }
            }

            // sensitivities; opposing face (penalty); use general shape
            // function for this single ip
            const scalar lhsFacO = penaltyIp * c_amag;
            for (label ic = 0; ic < opposingNodesPerFace; ++ic)
            {
                const label icNdim =
                    (o_face_node_ordinals[ic] + currentNodesPerElement);
                const scalar r = p_o_general_shape_function[ic];
                p_lhs[rowR + icNdim] -= r * lhsFacO;
            }

            // sensitivities; opposing element (diffusion)
            for (label ic = 0; ic < opposingNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const label icNdim = (ic + currentNodesPerElement);
                for (label l = 0; l < SPATIAL_DIM; ++l)
                {
                    const scalar nxl = p_oNx[l];
                    const scalar dndxl = p_o_dndx[offSetDnDx + l];

                    // -dphi/dxj*nj*dS
                    p_lhs[rowR + icNdim] -= -dndxl * nxl * c_amag / 2.0;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void wallScaleDiffusionAssembler::assembleElemTermsInterfaceSideNoSlipWall_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_yScale;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get transport fields/side fields
    const auto& yScaleSTKFieldRef = phi_->stkFieldRef();
    const auto& nodalSideYScaleSTKFieldRef =
        phi_->nodeSideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

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

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_yScale.resize(nodesPerElement);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_yScale = &ws_yScale[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_dndx = &ws_dndx[0];

        const stk::mesh::Bucket::size_type nBoundarySides = sideBucket.size();

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

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element; n/a
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label iNode = 0; iNode < numNodes; ++iNode)
            {
                stk::mesh::Entity node = elemNodeRels[iNode];

                // set connected nodes
                connectedNodes[iNode] = node;

                // gather scalars
                p_bcMultiplier[iNode] = 1.0;

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = iNode * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }

                // gather 1-dim fields
                p_yScale[iNode] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);

            // Fill gamma + overwrite phi values at boundary with dirichlet
            // values
            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                const label ic = faceNodeOrdinals[iNode];

                stk::mesh::Entity node = sideNodeRels[iNode];

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather; scalar
                p_yScale[ic] =
                    *stk::mesh::field_data(nodalSideYScaleSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coordinates[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coordinates[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over side ip's
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                const label faceOffSet = ip * SPATIAL_DIM;

                //================================
                // Diffusion: only diffusion
                //================================

                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[faceOffSet + j];
                        const scalar dndxj = p_dndx[offSetDnDx + j];

                        // matrix entries
                        label indexR = nearestNode;
                        label rowR = indexR * nodesPerElement;

                        const scalar yScale = p_yScale[ic];

                        // -dphi/dxj*Aj
                        scalar lhsfac = -dndxj * axj;
                        p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic];
                        p_rhs[indexR] -= lhsfac * yScale;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */

#endif /* HAS_INTERFACE */
