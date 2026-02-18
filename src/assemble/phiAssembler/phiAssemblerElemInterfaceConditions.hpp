// File       : phiAssemblerElemInterfaceConditions.hpp
// Created    : Thu Feb 29 2024 14:10:40 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Element based (stencil) assembly kernel implementation for
//              interfaces
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

#include "dgInfo.h"
#include "interface.h"

namespace accel
{

template <size_t N>
void phiAssembler<N>::assembleElemTermsInterfaces_(const domain* domain,
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

            assembleElemTermsInterfaceSide_(domain, interfaceSideInfoPtr, ctx);
        }
    }
}

template <size_t N>
void phiAssembler<N>::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    const bool includeAdv =
        ((transportMode_ != diffusion) &&
         (!interfaceSideInfoPtr->interfPtr()->isFluidSolidType()) &&
         (domain->type() == domainType::fluid));

    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const bool mediumIndependent = phi_->mediumIndependent();

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

    // c/o phi and normal flux
    std::vector<scalar> currentPhiBip(N);
    std::vector<scalar> opposingPhiBip(N);
    std::vector<scalar> currentDiffFluxBip(N);
    std::vector<scalar> opposingDiffFluxBip(N);

    // pointers to fixed values
    scalar* p_cNx = &cNx[0];
    scalar* p_oNx = &oNx[0];

    // nodal fields to gather
    std::vector<scalar> ws_c_face_phi;
    std::vector<scalar> ws_o_face_phi;
    std::vector<scalar> ws_c_elem_phi;
    std::vector<scalar> ws_o_elem_phi;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;
    std::vector<scalar> ws_c_Gamma;
    std::vector<scalar> ws_o_Gamma;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // Get transport fields/side fields
    const auto& phiSTKFieldRef = phi_->stkFieldRef();

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

    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType() &&
        !mediumIndependent)
    {
        // this is like zero-gradient
    }
    else
    {
        for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size());
             iSide++)
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
                const label opposingNodesPerFace =
                    meFCOpposing->nodesPerElement_;
                const label currentNodesPerElement =
                    meSCSCurrent->nodesPerElement_;
                const label opposingNodesPerElement =
                    meSCSOpposing->nodesPerElement_;

                // resize some things; matrix related
                const label totalNodes =
                    currentNodesPerElement + opposingNodesPerElement;
                const label lhsSize = totalNodes * N * totalNodes * N;
                const label rhsSize = totalNodes * N;
                lhs.resize(lhsSize);
                rhs.resize(rhsSize);
                scratchIds.resize(rhsSize);
                scratchVals.resize(rhsSize);
                connectedNodes.resize(totalNodes);

                // algorithm related; element; dndx will be at a
                // single gauss point...
                ws_c_elem_phi.resize(currentNodesPerElement * N);
                ws_o_elem_phi.resize(opposingNodesPerElement * N);
                ws_c_elem_coordinates.resize(currentNodesPerElement *
                                             SPATIAL_DIM);
                ws_o_elem_coordinates.resize(opposingNodesPerElement *
                                             SPATIAL_DIM);
                ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
                ws_c_det_j.resize(1);
                ws_o_det_j.resize(1);

                // algorithm related; face
                ws_c_face_phi.resize(currentNodesPerFace * N);
                ws_o_face_phi.resize(opposingNodesPerFace * N);
                ws_c_Gamma.resize(currentNodesPerFace);
                ws_o_Gamma.resize(opposingNodesPerFace);
                ws_c_general_shape_function.resize(currentNodesPerFace);
                ws_o_general_shape_function.resize(opposingNodesPerFace);

                // pointers
                scalar* p_lhs = &lhs[0];
                scalar* p_rhs = &rhs[0];

                scalar* p_c_face_phi = &ws_c_face_phi[0];
                scalar* p_o_face_phi = &ws_o_face_phi[0];
                scalar* p_c_elem_phi = &ws_c_elem_phi[0];
                scalar* p_o_elem_phi = &ws_o_elem_phi[0];
                scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];
                scalar* p_c_Gamma = &ws_c_Gamma[0];
                scalar* p_o_Gamma = &ws_o_Gamma[0];

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
                    p_c_Gamma[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                    // gather N-dim fields
                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        const label offSet = i * current_num_face_nodes + ni;
                        p_c_face_phi[offSet] = phi[i];
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
                    p_o_Gamma[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                    // gather N-dim fields
                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        const label offSet = i * opposing_num_face_nodes + ni;
                        p_o_face_phi[offSet] = phi[i];
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
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }

                    // gather N-dim fields
                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        p_c_elem_phi[ni * N + i] = phi[i];
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
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }

                    // gather N-dim fields
                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        p_o_elem_phi[ni * N + i] = phi[i];
                    }
                }

                // apply transformations
                interfaceSideInfoPtr->rotateVectorListCompact<N>(
                    ws_o_face_phi, opposing_num_face_nodes);
                interfaceSideInfoPtr->rotateVectorList<N>(
                    ws_o_elem_phi, opposing_num_elem_nodes);
                interfaceSideInfoPtr->transformCoordinateList(
                    ws_o_elem_coordinates, opposing_num_elem_nodes);

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
                meFCCurrent->interpolatePoint(N,
                                              &currentIsoParCoords[0],
                                              &ws_c_face_phi[0],
                                              &currentPhiBip[0]);

                meFCOpposing->interpolatePoint(N,
                                               &opposingIsoParCoords[0],
                                               &ws_o_face_phi[0],
                                               &opposingPhiBip[0]);

                scalar currentGammaBip = 0.0;
                meFCCurrent->interpolatePoint(1,
                                              &currentIsoParCoords[0],
                                              &ws_c_Gamma[0],
                                              &currentGammaBip);

                scalar opposingGammaBip = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_o_Gamma[0],
                                               &opposingGammaBip);

                // zero-out diffusion fluxes
                for (label i = 0; i < N; ++i)
                {
                    currentDiffFluxBip[i] = 0.0;
                    opposingDiffFluxBip[i] = 0.0;
                }

                // compute diffusion vector; current
                for (label ic = 0; ic < currentNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];

                        for (label i = 0; i < N; ++i)
                        {
                            const scalar phi = p_c_elem_phi[ic * N + i];

                            // -Gamma*dphi/dxj*Aj
                            currentDiffFluxBip[i] +=
                                -currentGammaBip * dndxj * nxj * phi;
                        }
                    }
                }

                // compute diffusion vector; opposing
                for (label ic = 0; ic < opposingNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_oNx[j];
                        const scalar dndxj = p_o_dndx[offSetDnDx + j];

                        for (label i = 0; i < N; ++i)
                        {
                            const scalar phi = p_o_elem_phi[ic * N + i];

                            // -Gamma*dphi/dxj*Aj
                            opposingDiffFluxBip[i] +=
                                -opposingGammaBip * dndxj * nxj * phi;
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
                const scalar tmDot =
                    includeAdv ? (stk::mesh::field_data(
                                     *mDotSideSTKFieldPtr_,
                                     currentFace))[currentGaussPointId]
                               : 0.0;
                const scalar abs_tmDot = std::abs(tmDot);

                // compute penalty
                const scalar penaltyIp =
                    penaltyFactor *
                    (currentGammaBip * currentInverseLength +
                     opposingGammaBip * opposingInverseLength) /
                    2.0;

                for (label i = 0; i < N; ++i)
                {
                    // non conformal diffusive flux
                    const scalar ncDiffFlux =
                        (currentDiffFluxBip[i] - opposingDiffFluxBip[i]) / 2.0;

                    // non conformal advection
                    const scalar ncAdv =
                        tmDot * (currentPhiBip[i] + opposingPhiBip[i]) / 2.0 +
                        abs_tmDot * (currentPhiBip[i] - opposingPhiBip[i]) /
                            2.0;

                    // assemble residual; form proper rhs index for
                    // current face assembly
                    const label indexR = nn * N + i;
                    p_rhs[indexR] -=
                        ((ncDiffFlux +
                          penaltyIp * (currentPhiBip[i] - opposingPhiBip[i])) *
                             c_amag +
                         ncAdv);

                    // set-up row for matrix
                    const label rowR = indexR * totalNodes * N;

                    // sensitivities; current face (penalty and
                    // advection); use general shape function for
                    // this single ip
                    const scalar lhsFacC =
                        penaltyIp * c_amag + (abs_tmDot + tmDot) / 2.0;
                    for (label ic = 0; ic < currentNodesPerFace; ++ic)
                    {
                        const label icNdim = c_face_node_ordinals[ic] * N;
                        const scalar r = p_c_general_shape_function[ic];
                        p_lhs[rowR + icNdim + i] += r * lhsFacC;
                    }

                    // sensitivities; current element (diffusion)
                    for (label ic = 0; ic < currentNodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            ic * SPATIAL_DIM; // single intg. point
                        const label icNdim = ic * N;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar nxj = p_cNx[j];
                            const scalar dndxj = p_c_dndx[offSetDnDx + j];

                            // -Gamma*dphi/dxj*nj*dS
                            p_lhs[rowR + icNdim + i] +=
                                -currentGammaBip * dndxj * nxj * c_amag / 2.0;
                        }
                    }

                    // sensitivities; opposing face (penalty and
                    // advection); use general shape function for
                    // this single ip
                    const scalar lhsFacO =
                        penaltyIp * c_amag + (abs_tmDot - tmDot) / 2.0;
                    for (label ic = 0; ic < opposingNodesPerFace; ++ic)
                    {
                        const label icNdim = (o_face_node_ordinals[ic] +
                                              currentNodesPerElement) *
                                             N;
                        const scalar r = p_o_general_shape_function[ic];
                        for (label j = 0; j < N; ++j)
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
                        const label icNdim = (ic + currentNodesPerElement) * N;
                        for (label j = 0; j < N; ++j)
                        {
                            for (label l = 0; l < SPATIAL_DIM; ++l)
                            {
                                const scalar nxl = p_oNx[l];
                                const scalar dndxl = p_o_dndx[offSetDnDx + l];

                                // -Gamma*dphi/dxj*nj*dS
                                p_lhs[rowR + icNdim + j] -=
                                    -opposingGammaBip * dndxl * nxl * c_amag /
                                    2.0 * rotMat(i, j);
                            }
                        }
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

} // namespace accel

#endif /* HAS_INTERFACE */
