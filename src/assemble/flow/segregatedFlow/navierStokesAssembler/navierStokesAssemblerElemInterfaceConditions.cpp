// File       : navierStokesAssemblerElemInterfaceConditions.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

// code
#include "dgInfo.h"
#include "flowModel.h"
#include "interface.h"
#include "mesh.h"
#include "navierStokesAssembler.h"

namespace accel
{

void navierStokesAssembler::assembleElemTermsInterfaces_(const domain* domain,
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
                // If this is a fluid/solid interface, deal with as a no-slip
                // wall
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

void navierStokesAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const bool slipNonOverlap =
        interfaceSideInfoPtr->interfPtr()->isSlipNonOverlap();
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

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
    std::vector<scalar> opposingElementIsoParCoordsTrans(SPATIAL_DIM);

    // c/o velocity and normal flux
    std::vector<scalar> currentUBip(SPATIAL_DIM);
    std::vector<scalar> opposingUBip(SPATIAL_DIM);
    std::vector<scalar> currentDiffFluxBip(SPATIAL_DIM);
    std::vector<scalar> opposingDiffFluxBip(SPATIAL_DIM);

    // interpolate nodal values to point-in-elem
    const label sizeOfScalarField = 1;
    const label sizeOfVectorField = SPATIAL_DIM;

    // pointers to fixed values
    scalar* p_cNx = &cNx[0];
    scalar* p_oNx = &oNx[0];

    // nodal fields to gather
    std::vector<scalar> ws_c_face_velocity;
    std::vector<scalar> ws_o_face_velocity;
    std::vector<scalar> ws_c_elem_velocity;
    std::vector<scalar> ws_o_elem_velocity;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;
    std::vector<scalar> ws_c_muEff;
    std::vector<scalar> ws_o_muEff;
    std::vector<scalar> ws_bcMultiplier;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // Get transport fields/side fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto* nodalSideUSTKFieldPtr =
        interfaceSideInfoPtr->hasNonoverlap_ && !slipNonOverlap
            ? model_->URef().nodeSideFieldRef().stkFieldPtr()
            : nullptr;

    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();

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

        // now loop over all the DgInfo objects on this particular
        // exposed face
        for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
        {
            dgInfo* dgInfo = faceDgInfoVec[k];

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

            // if gauss point is exposed (non-overlapping), then
            // treat as a wall. Currently the interface is assumed not moving.
            // So the non-overlap portition will not account for any movement
            if (dgInfo->gaussPointExposed_)
            {
                if (slipNonOverlap)
                {
                    // mapping from ip to nodes for this ordinal
                    const label* ipNodeMap =
                        meSCSCurrent->ipNodeMap(currentFaceOrdinal);

                    // extract nearset node
                    const label nearestNode = ipNodeMap[currentGaussPointId];

                    // extract some master element info
                    const label currentNodesPerSide =
                        meFCCurrent->nodesPerElement_;
                    const label currentNodesPerElement =
                        meSCSCurrent->nodesPerElement_;

                    // resize some things; matrix related
                    const label lhsSize = currentNodesPerElement * SPATIAL_DIM *
                                          currentNodesPerElement * SPATIAL_DIM;
                    const label rhsSize = currentNodesPerElement * SPATIAL_DIM;
                    lhs.resize(lhsSize);
                    rhs.resize(rhsSize);
                    scratchIds.resize(rhsSize);
                    scratchVals.resize(rhsSize);
                    connectedNodes.resize(currentNodesPerElement);

                    // algorithm related; element; dndx will be at a
                    // single gauss point...
                    ws_c_elem_velocity.resize(currentNodesPerElement *
                                              SPATIAL_DIM);
                    ws_c_elem_coordinates.resize(currentNodesPerElement *
                                                 SPATIAL_DIM);

                    ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                    ws_c_det_j.resize(1);

                    // algorithm related; face
                    ws_c_muEff.resize(currentNodesPerSide);
                    ws_c_general_shape_function.resize(currentNodesPerSide);

                    // pointers
                    scalar* p_lhs = &lhs[0];
                    scalar* p_rhs = &rhs[0];

                    scalar* p_c_elem_velocity = &ws_c_elem_velocity[0];
                    scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                    scalar* p_c_muEff = &ws_c_muEff[0];

                    // me pointers
                    scalar* p_c_general_shape_function =
                        &ws_c_general_shape_function[0];
                    scalar* p_c_dndx = &ws_c_dndx[0];

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
                        p_c_muEff[ni] =
                            *stk::mesh::field_data(muEffSTKFieldRef, node);
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

                        // gather; vector
                        const scalar* U =
                            stk::mesh::field_data(USTKFieldRef, node);
                        const scalar* coords =
                            stk::mesh::field_data(coordsSTKFieldRef, node);
                        const label niNdim = ni * SPATIAL_DIM;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_c_elem_velocity[niNdim + i] = U[i];
                            p_c_elem_coordinates[niNdim + i] = coords[i];
                        }
                    }

                    // pointer to face data
                    const scalar* c_areaVec = stk::mesh::field_data(
                        exposedAreaVecSTKFieldRef, currentFace);

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
                            c_areaVec[currentGaussPointId * SPATIAL_DIM + i] /
                            c_amag;
                    }

                    // project from side to element; method deals with
                    // the -1:1 isInElement range to the proper
                    // underlying CVFEM range
                    meSCSCurrent->sidePcoords_to_elemPcoords(
                        currentFaceOrdinal,
                        1,
                        &currentIsoParCoords[0],
                        &currentElementIsoParCoords[0]);

                    // compute dndx
                    scalar scs_error = 0.0;
                    meSCSCurrent->general_face_grad_op(
                        currentFaceOrdinal,
                        &currentElementIsoParCoords[0],
                        &p_c_elem_coordinates[0],
                        &p_c_dndx[0],
                        &ws_c_det_j[0],
                        &scs_error);

                    // interpolate face data; current only ...
                    scalar currentMuEffBip = 0.0;
                    meFCCurrent->interpolatePoint(sizeOfScalarField,
                                                  &currentIsoParCoords[0],
                                                  &ws_c_muEff[0],
                                                  &currentMuEffBip);

                    // zero lhs/rhs
                    for (label p = 0; p < lhsSize; ++p)
                    {
                        p_lhs[p] = 0.0;
                    }
                    for (label p = 0; p < rhsSize; ++p)
                    {
                        p_rhs[p] = 0.0;
                    }

                    // zero tangential stress
                    for (label ic = 0; ic < currentNodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            ic * SPATIAL_DIM; // single intg. point
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj =
                                c_areaVec[currentGaussPointId * SPATIAL_DIM +
                                          j];
                            const scalar dndxj = p_c_dndx[offSetDnDx + j];
                            const scalar uxj =
                                p_c_elem_velocity[ic * SPATIAL_DIM + j];
                            const scalar divUstress = 2.0 / 3.0 *
                                                      currentMuEffBip * dndxj *
                                                      uxj * axj * comp;
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                /* matrix entries */
                                label indexR = nearestNode * SPATIAL_DIM + i;
                                label rowR = indexR * currentNodesPerElement *
                                             SPATIAL_DIM;
                                const scalar dndxi = p_c_dndx[offSetDnDx + i];
                                const scalar uxi =
                                    p_c_elem_velocity[ic * SPATIAL_DIM + i];
                                const scalar nxi = p_cNx[i];
                                const scalar nxinxi = nxi * nxi;

                                /* -mu*dui/dxj*Aj*ni*ni; sneak in divU
                                 * (explicit): */
                                scalar lhsfac =
                                    -currentMuEffBip * dndxj * axj * nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;
                                p_rhs[indexR] -=
                                    lhsfac * uxi + divUstress * nxinxi;

                                /* -mu*duj/dxi*Aj*ni*ni */
                                lhsfac =
                                    -currentMuEffBip * dndxi * axj * nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;
                                p_rhs[indexR] -= lhsfac * uxj;

                                /* now we need the +nx*ny*Fy + nx*nz*Fz part */
                                for (label l = 0; l < SPATIAL_DIM; ++l)
                                {
                                    if (i != l)
                                    {
                                        const scalar nxinxl = nxi * p_cNx[l];
                                        const scalar uxl =
                                            p_c_elem_velocity[ic * SPATIAL_DIM +
                                                              l];
                                        const scalar dndxl =
                                            p_c_dndx[offSetDnDx + l];

                                        /* -ni*nl*mu*dul/dxj*Aj; sneak in divU
                                         * (explicit):
                                         */
                                        lhsfac = -currentMuEffBip * dndxj *
                                                 axj * nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + l] +=
                                            lhsfac;
                                        p_rhs[indexR] -=
                                            lhsfac * uxl + divUstress * nxinxl;

                                        /* -ni*nl*mu*duj/dxl*Aj */
                                        lhsfac = -currentMuEffBip * dndxl *
                                                 axj * nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + j] +=
                                            lhsfac;
                                        p_rhs[indexR] -= lhsfac * uxj;
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    // mapping from ip to nodes for this ordinal
                    const label* ipNodeMap =
                        meSCSCurrent->ipNodeMap(currentFaceOrdinal);

                    // extract nearset node
                    const label nearestNode = ipNodeMap[currentGaussPointId];

                    // extract some master element info
                    const label currentNodesPerSide =
                        meFCCurrent->nodesPerElement_;
                    const label currentNodesPerElement =
                        meSCSCurrent->nodesPerElement_;

                    // resize some things; matrix related
                    const label lhsSize = currentNodesPerElement * SPATIAL_DIM *
                                          currentNodesPerElement * SPATIAL_DIM;
                    const label rhsSize = currentNodesPerElement * SPATIAL_DIM;
                    lhs.resize(lhsSize);
                    rhs.resize(rhsSize);
                    scratchIds.resize(rhsSize);
                    scratchVals.resize(rhsSize);
                    connectedNodes.resize(currentNodesPerElement);

                    // algorithm related; element; dndx will be at a
                    // single gauss point...
                    ws_c_elem_velocity.resize(currentNodesPerElement *
                                              SPATIAL_DIM);
                    ws_c_elem_coordinates.resize(currentNodesPerElement *
                                                 SPATIAL_DIM);
                    ws_bcMultiplier.resize(currentNodesPerElement);
                    ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                    ws_c_det_j.resize(1);

                    // algorithm related; face
                    ws_c_muEff.resize(currentNodesPerSide);
                    ws_c_general_shape_function.resize(currentNodesPerSide);

                    // pointers
                    scalar* p_lhs = &lhs[0];
                    scalar* p_rhs = &rhs[0];

                    scalar* p_c_elem_velocity = &ws_c_elem_velocity[0];
                    scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                    scalar* p_c_muEff = &ws_c_muEff[0];
                    scalar* p_bcMultiplier = &ws_bcMultiplier[0];

                    // me pointers
                    scalar* p_c_general_shape_function =
                        &ws_c_general_shape_function[0];
                    scalar* p_c_dndx = &ws_c_dndx[0];

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

                        p_bcMultiplier[ni] = 1.0;

                        // gather; vector
                        const scalar* U =
                            stk::mesh::field_data(USTKFieldRef, node);
                        const scalar* coords =
                            stk::mesh::field_data(coordsSTKFieldRef, node);
                        const label niNdim = ni * SPATIAL_DIM;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_c_elem_velocity[niNdim + i] = U[i];
                            p_c_elem_coordinates[niNdim + i] = coords[i];
                        }
                    }

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
                        const label icnn = c_face_node_ordinals[ni];

                        stk::mesh::Entity node = current_face_node_rels[ni];

                        // gather; scalar
                        p_c_muEff[ni] =
                            *stk::mesh::field_data(muEffSTKFieldRef, node);

                        p_bcMultiplier[icnn] = 0.0;

                        // gather; vector
                        const scalar* U =
                            stk::mesh::field_data(*nodalSideUSTKFieldPtr, node);
                        const label niNdim = icnn * SPATIAL_DIM;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_c_elem_velocity[niNdim + i] = U[i];
                        }
                    }

                    // pointer to face data
                    const scalar* c_areaVec = stk::mesh::field_data(
                        exposedAreaVecSTKFieldRef, currentFace);

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
                            c_areaVec[currentGaussPointId * SPATIAL_DIM + i] /
                            c_amag;
                    }

                    // project from side to element; method deals with
                    // the -1:1 isInElement range to the proper
                    // underlying CVFEM range
                    meSCSCurrent->sidePcoords_to_elemPcoords(
                        currentFaceOrdinal,
                        1,
                        &currentIsoParCoords[0],
                        &currentElementIsoParCoords[0]);

                    // compute dndx
                    scalar scs_error = 0.0;
                    meSCSCurrent->general_face_grad_op(
                        currentFaceOrdinal,
                        &currentElementIsoParCoords[0],
                        &p_c_elem_coordinates[0],
                        &p_c_dndx[0],
                        &ws_c_det_j[0],
                        &scs_error);

                    // interpolate face data; current only ...
                    scalar currentMuEffBip = 0.0;
                    meFCCurrent->interpolatePoint(sizeOfScalarField,
                                                  &currentIsoParCoords[0],
                                                  &ws_c_muEff[0],
                                                  &currentMuEffBip);

                    // zero lhs/rhs
                    for (label p = 0; p < lhsSize; ++p)
                    {
                        p_lhs[p] = 0.0;
                    }
                    for (label p = 0; p < rhsSize; ++p)
                    {
                        p_rhs[p] = 0.0;
                    }

                    // stress: zero normal stress
                    for (label ic = 0; ic < currentNodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            ic * SPATIAL_DIM; // single intg. point
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj =
                                c_areaVec[currentGaussPointId * SPATIAL_DIM +
                                          j];
                            const scalar dndxj = p_c_dndx[offSetDnDx + j];
                            const scalar uxj =
                                p_c_elem_velocity[ic * SPATIAL_DIM + j];
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                /* matrix entries */
                                label indexR = nearestNode * SPATIAL_DIM + i;
                                label rowR = indexR * currentNodesPerElement *
                                             SPATIAL_DIM;
                                const scalar dndxi = p_c_dndx[offSetDnDx + i];
                                const scalar uxi =
                                    p_c_elem_velocity[ic * SPATIAL_DIM + i];
                                const scalar nxi = p_cNx[i];
                                const scalar om_nxinxi = 1.0 - nxi * nxi;

                                /* -mu*dui/dxj*Aj(1.0-nini); sneak in divU
                                 * (explicit) */
                                scalar lhsfac =
                                    -currentMuEffBip * dndxj * axj * om_nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + i] +=
                                    lhsfac * p_bcMultiplier[ic];
                                p_rhs[indexR] -= lhsfac * uxi;

                                /* -mu*duj/dxi*Aj(1.0-nini) */
                                lhsfac =
                                    -currentMuEffBip * dndxi * axj * om_nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + j] +=
                                    lhsfac * p_bcMultiplier[ic];
                                p_rhs[indexR] -= lhsfac * uxj;

                                /* now we need the -nx*ny*Fy - nx*nz*Fz part */
                                for (label l = 0; l < SPATIAL_DIM; ++l)
                                {
                                    if (i != l)
                                    {
                                        const scalar nxinxl = nxi * p_cNx[l];
                                        const scalar uxl =
                                            p_c_elem_velocity[ic * SPATIAL_DIM +
                                                              l];
                                        const scalar dndxl =
                                            p_c_dndx[offSetDnDx + l];

                                        /* +ni*nl*mu*dul/dxj*Aj; sneak in divU
                                         * (explicit):
                                         */
                                        lhsfac = currentMuEffBip * dndxj * axj *
                                                 nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + l] +=
                                            lhsfac * p_bcMultiplier[ic];
                                        p_rhs[indexR] -= lhsfac * uxl;

                                        /* +ni*nl*mu*duj/dxl*Aj */
                                        lhsfac = currentMuEffBip * dndxl * axj *
                                                 nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + j] +=
                                            lhsfac * p_bcMultiplier[ic];
                                        p_rhs[indexR] -= lhsfac * uxj;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // mapping from ip to nodes for this ordinal
                const label* ipNodeMap =
                    meSCSCurrent->ipNodeMap(currentFaceOrdinal);

                // extract some master element info
                const label currentNodesPerSide = meFCCurrent->nodesPerElement_;
                const label opposingNodesPerSide =
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

                // algorithm related; element; dndx will be at a
                // single gauss point...
                ws_c_elem_velocity.resize(currentNodesPerElement * SPATIAL_DIM);
                ws_o_elem_velocity.resize(opposingNodesPerElement *
                                          SPATIAL_DIM);
                ws_c_elem_coordinates.resize(currentNodesPerElement *
                                             SPATIAL_DIM);
                ws_o_elem_coordinates.resize(opposingNodesPerElement *
                                             SPATIAL_DIM);
                ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
                ws_c_det_j.resize(1);
                ws_o_det_j.resize(1);

                // algorithm related; face
                ws_c_face_velocity.resize(currentNodesPerSide * SPATIAL_DIM);
                ws_o_face_velocity.resize(opposingNodesPerSide * SPATIAL_DIM);
                ws_c_muEff.resize(currentNodesPerSide);
                ws_o_muEff.resize(opposingNodesPerSide);
                ws_c_general_shape_function.resize(currentNodesPerSide);
                ws_o_general_shape_function.resize(opposingNodesPerSide);

                // pointers
                scalar* p_lhs = &lhs[0];
                scalar* p_rhs = &rhs[0];

                scalar* p_c_face_velocity = &ws_c_face_velocity[0];
                scalar* p_o_face_velocity = &ws_o_face_velocity[0];
                scalar* p_c_elem_velocity = &ws_c_elem_velocity[0];
                scalar* p_o_elem_velocity = &ws_o_elem_velocity[0];
                scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];
                scalar* p_c_muEff = &ws_c_muEff[0];
                scalar* p_o_muEff = &ws_o_muEff[0];

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
                    p_c_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);

                    // gather; vector
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label offSet = i * current_num_face_nodes + ni;
                        p_c_face_velocity[offSet] = U[i];
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

                    // gather; scalar
                    p_o_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);

                    // gather; vector
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label offSet = i * opposing_num_face_nodes + ni;
                        p_o_face_velocity[offSet] = U[i];
                    }
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

                    // gather; vector
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label niNdim = ni * SPATIAL_DIM;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_elem_velocity[niNdim + i] = U[i];
                        p_c_elem_coordinates[niNdim + i] = coords[i];
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

                    // gather; vector
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label niNdim = ni * SPATIAL_DIM;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_velocity[niNdim + i] = U[i];
                        p_o_elem_coordinates[niNdim + i] = coords[i];
                    }
                }

                // apply transformations
                interfaceSideInfoPtr->rotateVectorListCompact<SPATIAL_DIM>(
                    ws_o_face_velocity, opposing_num_face_nodes);
                interfaceSideInfoPtr->rotateVectorList<SPATIAL_DIM>(
                    ws_o_elem_velocity, opposing_num_elem_nodes);
                interfaceSideInfoPtr->transformCoordinateList(
                    ws_o_elem_coordinates, opposing_num_elem_nodes);

                // pointer to face data
                const scalar* c_areaVec = stk::mesh::field_data(
                    exposedAreaVecSTKFieldRef, currentFace);
                const scalar* ncmDot =
                    stk::mesh::field_data(mDotSideSTKFieldRef, currentFace);

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
                        c_areaVec[currentGaussPointId * SPATIAL_DIM + i] /
                        c_amag;
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
                    &opposingElementIsoParCoordsTrans[0]);

                // compute dndx
                scalar scs_error = 0.0;
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
                meFCCurrent->interpolatePoint(sizeOfVectorField,
                                              &currentIsoParCoords[0],
                                              &ws_c_face_velocity[0],
                                              &currentUBip[0]);

                meFCOpposing->interpolatePoint(sizeOfVectorField,
                                               &opposingIsoParCoords[0],
                                               &ws_o_face_velocity[0],
                                               &opposingUBip[0]);

                scalar currentMuEffBip = 0.0;
                meFCCurrent->interpolatePoint(sizeOfScalarField,
                                              &currentIsoParCoords[0],
                                              &ws_c_muEff[0],
                                              &currentMuEffBip);

                scalar opposingMuEffBip = 0.0;
                meFCOpposing->interpolatePoint(sizeOfScalarField,
                                               &opposingIsoParCoords[0],
                                               &ws_o_muEff[0],
                                               &opposingMuEffBip);

                // zero-out diffusion fluxes
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    currentDiffFluxBip[j] = 0;
                    opposingDiffFluxBip[j] = 0;
                }

                // compute viscous stress tensor; current
                for (label ic = 0; ic < currentNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];
                        const scalar uxj =
                            p_c_elem_velocity[ic * SPATIAL_DIM + j];

                        const scalar divUstress = 2.0 / 3.0 * currentMuEffBip *
                                                  dndxj * uxj * nxj * comp;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dndxi = p_c_dndx[offSetDnDx + i];
                            const scalar uxi =
                                p_c_elem_velocity[ic * SPATIAL_DIM + i];

                            // -mu*dui/dxj*Aj with divU
                            currentDiffFluxBip[i] +=
                                -currentMuEffBip * dndxj * nxj * uxi +
                                divUstress;

                            // -mu*duj/dxi*Aj
                            currentDiffFluxBip[i] +=
                                -currentMuEffBip * dndxi * nxj * uxj;
                        }
                    }
                }

                // compute viscous stress tensor; opposing
                for (label ic = 0; ic < opposingNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_oNx[j];
                        const scalar dndxj = p_o_dndx[offSetDnDx + j];
                        const scalar uxj =
                            p_o_elem_velocity[ic * SPATIAL_DIM + j];

                        const scalar divUstress = 2.0 / 3.0 * opposingMuEffBip *
                                                  dndxj * uxj * nxj * comp;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dndxi = p_o_dndx[offSetDnDx + i];
                            const scalar uxi =
                                p_o_elem_velocity[ic * SPATIAL_DIM + i];

                            // -mu*dui/dxj*Aj with divU
                            opposingDiffFluxBip[i] +=
                                -opposingMuEffBip * dndxj * nxj * uxi +
                                divUstress;

                            // -mu*duj/dxi*Aj
                            opposingDiffFluxBip[i] +=
                                -opposingMuEffBip * dndxi * nxj * uxj;
                        }
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
                meFCCurrent->general_shape_fcn(1,
                                               &currentIsoParCoords[0],
                                               &ws_c_general_shape_function[0]);
                meFCOpposing->general_shape_fcn(
                    1,
                    &opposingIsoParCoords[0],
                    &ws_o_general_shape_function[0]);

                // save mDot
                const scalar tmDot = ncmDot[currentGaussPointId];
                const scalar abs_tmDot = std::abs(tmDot);

                // compute penalty
                const scalar penaltyIp =
                    penaltyFactor * 0.5 *
                    (currentMuEffBip * currentInverseLength +
                     opposingMuEffBip * opposingInverseLength);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    // non conformal diffusive flux
                    const scalar ncDiffFlux =
                        (currentDiffFluxBip[i] - opposingDiffFluxBip[i]) / 2.0;

                    // non conformal advection
                    const scalar ncAdv =
                        tmDot * (currentUBip[i] + opposingUBip[i]) / 2.0 +
                        abs_tmDot * (currentUBip[i] - opposingUBip[i]) / 2.0;

                    // assemble residual; form proper rhs index for
                    // current face assembly
                    const label indexR = nn * SPATIAL_DIM + i;
                    p_rhs[indexR] -=
                        ((ncDiffFlux +
                          penaltyIp * (currentUBip[i] - opposingUBip[i])) *
                             c_amag +
                         ncAdv);

                    // set-up row for matrix
                    const label rowR = indexR * totalNodes * SPATIAL_DIM;

                    // sensitivities; current face (penalty and
                    // advection); use general shape function for
                    // this single ip
                    const scalar lhsFacC =
                        penaltyIp * c_amag + (abs_tmDot + tmDot) / 2.0;
                    for (label ic = 0; ic < currentNodesPerSide; ++ic)
                    {
                        const label icNdim =
                            c_face_node_ordinals[ic] * SPATIAL_DIM;
                        const scalar r = p_c_general_shape_function[ic];
                        p_lhs[rowR + icNdim + i] += r * lhsFacC;
                    }

                    // sensitivities; current element (diffusion)
                    for (label ic = 0; ic < currentNodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            ic * SPATIAL_DIM; // single intg. point
                        const label icNdim = ic * SPATIAL_DIM;
                        const scalar dndxi = p_c_dndx[offSetDnDx + i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar nxj = p_cNx[j];
                            const scalar dndxj = p_c_dndx[offSetDnDx + j];

                            // -mu*dui/dxj*nj*dS (divU neglected)
                            p_lhs[rowR + icNdim + i] +=
                                -currentMuEffBip * dndxj * nxj * c_amag / 2.0;

                            // -mu*duj/dxi*nj*dS
                            p_lhs[rowR + icNdim + j] +=
                                -currentMuEffBip * dndxi * nxj * c_amag / 2.0;
                        }
                    }

                    // sensitivities; opposing face (penalty and
                    // advection); use general shape function for
                    // this single ip
                    const scalar lhsFacO =
                        penaltyIp * c_amag + (abs_tmDot - tmDot) / 2.0;
                    for (label ic = 0; ic < opposingNodesPerSide; ++ic)
                    {
                        const label icNdim = (o_face_node_ordinals[ic] +
                                              currentNodesPerElement) *
                                             SPATIAL_DIM;
                        const scalar r = p_o_general_shape_function[ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_lhs[rowR + icNdim + j] -=
                                r * lhsFacO * rotMat(i, j);
                        }
                    }

                    // sensitivities; opposing element (diffusion)
                    for (label ic = 0; ic < opposingNodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            ic * SPATIAL_DIM; // single intg. point
                        const label icNdim =
                            (ic + currentNodesPerElement) * SPATIAL_DIM;
                        const scalar dndxi = p_o_dndx[offSetDnDx + i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            for (label l = 0; l < SPATIAL_DIM; ++l)
                            {
                                const scalar nxl = p_oNx[l];
                                const scalar dndxl = p_o_dndx[offSetDnDx + l];

                                // -mu*dui/dxj*nj*dS (divU neglected)
                                p_lhs[rowR + icNdim + j] -=
                                    -opposingMuEffBip * dndxl * nxl * c_amag /
                                    2.0 * rotMat(i, j);

                                // -mu*duj/dxi*nj*dS
                                p_lhs[rowR + icNdim + j] -=
                                    -opposingMuEffBip * dndxi * nxl * c_amag /
                                    2.0 * rotMat(l, j);
                            }
                        }
                    }
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsInterfaceSideNoSlipWall_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    if (domain->turbulence_.option_ == turbulenceOption::laminar)
    {
        // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM
        // and nodesPerElem*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values; both boundary and opposing surface
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_nx = &nx[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_coords;
        std::vector<scalar> ws_muEff;
        std::vector<scalar> ws_bcMultiplier;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
        const auto& nodalSideUSTKFieldRef =
            model_->URef().nodeSideFieldRef().stkFieldRef();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // define vector of parent topos; should always be UNITY in size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();

        // shifted ip's for gradients?
        const bool isUGradientShifted = model_->URef().isGradientShifted();

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
            const label lhsSize =
                nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
            const label rhsSize = nodesPerElement * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerElement);

            // algorithm related; element/face
            ws_U.resize(nodesPerElement * SPATIAL_DIM);
            ws_coords.resize(nodesPerElement * SPATIAL_DIM);
            ws_muEff.resize(nodesPerSide);
            ws_bcMultiplier.resize(nodesPerElement);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
            ws_det_j.resize(numScsBip);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_bcMultiplier = &ws_bcMultiplier[0];
            scalar* p_coords = &ws_coords[0];
            scalar* p_muEff = &ws_muEff[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];
            scalar* p_dndx = &ws_dndx[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
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

                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                // extract the connected element to this exposed face; should be
                // single in size!
                stk::mesh::Entity const* faceElemRels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);

                // get element; its face ordinal number and populate
                // face_node_ordinals
                stk::mesh::Entity element = faceElemRels[0];
                const label faceOrdinal =
                    bulkData.begin_element_ordinals(side)[0];

                // mapping from ip to nodes for this ordinal;
                const label* ipNodeMap =
                    meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

                // populate faceNodeOrdinals
                const label* faceNodeOrdinals =
                    meSCS->side_node_ordinals(faceOrdinal);

                //==========================================
                // gather nodal data off of element
                //==========================================
                stk::mesh::Entity const* elemNodeRels =
                    bulkData.begin_nodes(element);
                label numNodes = bulkData.num_nodes(element);

                // sanity check on num nodes
                STK_ThrowAssert(numNodes == nodesPerElement);
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = elemNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_bcMultiplier[ni] = 1.0;

                    // gather vectors
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                        p_coords[offSet + j] = coords[j];
                    }
                }

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    const label ic = faceNodeOrdinals[ni];

                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);

                    // set 0 the boundary nodes
                    p_bcMultiplier[ic] = 0.0;

                    // gather vectors
                    scalar* U =
                        stk::mesh::field_data(nodalSideUSTKFieldRef, node);

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        // overwrite velocity value at boundary with the
                        // dirichlet value
                        p_U[ic * SPATIAL_DIM + i] = U[i];
                    }
                }

                // compute dndx
                scalar scs_error = 0.0;
                if (isUGradientShifted)
                {
                    meSCS->shifted_face_grad_op(1,
                                                faceOrdinal,
                                                &p_coords[0],
                                                &p_dndx[0],
                                                &ws_det_j[0],
                                                &scs_error);
                }
                else
                {
                    meSCS->face_grad_op(1,
                                        faceOrdinal,
                                        &p_coords[0],
                                        &p_dndx[0],
                                        &ws_det_j[0],
                                        &scs_error);
                }

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // offset for bip area vector and types of shape function
                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    // zero out vector quantities
                    scalar asq = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[faceOffSet + j];
                        asq += axj * axj;
                    }
                    const scalar amag = std::sqrt(asq);

                    // interpolate to bip
                    scalar muEffBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        muEffBip += r * p_muEff[ic];
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_nx[i] = areaVec[faceOffSet + i] / amag;
                    }

                    //================================
                    // stress: zero normal stress
                    //================================
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[faceOffSet + j];
                            const scalar dndxj = p_dndx[offSetDnDx + j];
                            const scalar uxj = p_U[ic * SPATIAL_DIM + j];
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                label indexR = nearestNode * SPATIAL_DIM + i;
                                label rowR =
                                    indexR * nodesPerElement * SPATIAL_DIM;
                                const scalar dndxi = p_dndx[offSetDnDx + i];
                                const scalar uxi = p_U[ic * SPATIAL_DIM + i];
                                const scalar nxi = p_nx[i];
                                const scalar om_nxinxi = 1.0 - nxi * nxi;

                                // -mu*dui/dxj*Aj(1.0-nini)
                                scalar lhsfac =
                                    -muEffBip * dndxj * axj * om_nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + i] +=
                                    lhsfac * p_bcMultiplier[ic];
                                p_rhs[indexR] -= lhsfac * uxi;

                                // -mu*duj/dxi*Aj(1.0-nini)
                                lhsfac = -muEffBip * dndxi * axj * om_nxinxi;
                                p_lhs[rowR + ic * SPATIAL_DIM + j] +=
                                    lhsfac * p_bcMultiplier[ic];
                                p_rhs[indexR] -= lhsfac * uxj;

                                // now we need the -nx*ny*Fy - nx*nz*Fz part
                                for (label l = 0; l < SPATIAL_DIM; ++l)
                                {
                                    if (i != l)
                                    {
                                        const scalar nxinxl = nxi * p_nx[l];
                                        const scalar uxl =
                                            p_U[ic * SPATIAL_DIM + l];
                                        const scalar dndxl =
                                            p_dndx[offSetDnDx + l];

                                        // +ni*nl*mu*dul/dxj*Aj
                                        lhsfac =
                                            muEffBip * dndxj * axj * nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + l] +=
                                            lhsfac * p_bcMultiplier[ic];
                                        p_rhs[indexR] -= lhsfac * uxl;

                                        // +ni*nl*mu*duj/dxl*Aj
                                        lhsfac =
                                            muEffBip * dndxl * axj * nxinxl;
                                        p_lhs[rowR + ic * SPATIAL_DIM + j] +=
                                            lhsfac * p_bcMultiplier[ic];
                                        p_rhs[indexR] -= lhsfac * uxj;
                                    }
                                }
                            }
                        }
                    }
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
        // space for LHS/RHS; nodesPerSide*SPATIAL_DIM*nodesPerSide*SPATIAL_DIM
        // and nodesPerSide*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // bip values
        std::vector<scalar> uBip(SPATIAL_DIM);
        std::vector<scalar> unitNormal(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_uBip = &uBip[0];
        scalar* p_unitNormal = &unitNormal[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& sideUSTKFieldRef =
            model_->URef().sideFieldRef().stkFieldRef();

        const auto& uWallCoeffsSTKFieldRef =
            model_->uWallCoeffsRef().stkFieldRef();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();

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
            ws_U.resize(nodesPerSide * SPATIAL_DIM);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nBoundaryFaces =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iBFace = 0;
                 iBFace < nBoundaryFaces;
                 ++iBFace)
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

                // get face
                stk::mesh::Entity side = sideBucket[iBFace];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* face_node_rels =
                    bulkData.begin_nodes(side);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = face_node_rels[ni];
                    connectedNodes[ni] = node;

                    // gather vectors
                    scalar* U = stk::mesh::field_data(USTKFieldRef, node);

                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                    }
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const scalar* uWallCoeffsBip =
                    stk::mesh::field_data(uWallCoeffsSTKFieldRef, side);
                const scalar* UbcVec =
                    stk::mesh::field_data(sideUSTKFieldRef, side);

                // loop over face nodes
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    const label nearestNode = faceIpNodeMap[ip];

                    // zero out vector quantities; squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] = 0.0;

                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    // form unit normal
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_unitNormal[j] = areaVec[offSetAveraVec + j] / aMag;
                    }

                    // interpolate to bip
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        const label offSetFN = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_uBip[j] += r * p_U[offSetFN + j];
                        }
                    }

                    //================================
                    // stress: zero normal stress
                    //================================
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        label indexR = nearestNode * SPATIAL_DIM + i;
                        label rowR = indexR * nodesPerSide * SPATIAL_DIM;

                        scalar uiTan = 0.0;
                        scalar uiBcTan = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar ninj =
                                p_unitNormal[i] * p_unitNormal[j];
                            if (i == j)
                            {
                                const scalar om_nini = 1.0 - ninj;
                                uiTan += om_nini * p_uBip[j];
                                uiBcTan +=
                                    om_nini * UbcVec[ip * SPATIAL_DIM + j];
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    p_lhs[rowR + ic * SPATIAL_DIM + i] +=
                                        uWallCoeffsBip[ip] * om_nini *
                                        p_velocity_face_shape_function
                                            [offSetSF_face + ic];
                                }
                            }
                            else
                            {
                                uiTan -= ninj * p_uBip[j];
                                uiBcTan -= ninj * UbcVec[ip * SPATIAL_DIM + j];
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    p_lhs[rowR + ic * SPATIAL_DIM + j] -=
                                        uWallCoeffsBip[ip] * ninj *
                                        p_velocity_face_shape_function
                                            [offSetSF_face + ic];
                                }
                            }
                        }
                        p_rhs[indexR] -= uWallCoeffsBip[ip] * (uiTan - uiBcTan);
                    }
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

} // namespace accel

#endif /* HAS_INTERFACE */
