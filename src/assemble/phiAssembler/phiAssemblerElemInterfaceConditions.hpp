// File       : phiAssemblerElemInterfaceConditions.hpp
// Created    : Thu Feb 29 2024 14:10:40 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Element based (stencil) assembly kernel implementation for
//              interfaces
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

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
            // Internal interface (same domain): scatter both sides.
            assembleElemTermsInterfaceSide_(
                domain, interf->masterInfoPtr(), ctx);
            assembleElemTermsInterfaceSide_(
                domain, interf->slaveInfoPtr(), ctx);
        }
        else
        {
            // Cross-domain penalty / DG path: each domain scatters its
            // own side (existing behaviour).
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

    // two-sided interface treatment only if both zones allow it
    const label opposingZone =
        interfaceSideInfoPtr->interfPtr()->masterZoneIndex() == domain->index()
            ? interfaceSideInfoPtr->interfPtr()->slaveZoneIndex()
            : interfaceSideInfoPtr->interfPtr()->masterZoneIndex();
    const bool mediumIndependent = phi_->interfaceTwoSided(domain->index()) &&
                                   phi_->interfaceTwoSided(opposingZone);

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

    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType() &&
        !mediumIndependent)
    {
        // this is like zero-gradient
    }
    else
    {
        const scalar multiplier = 0.5;

        for (auto& faceIpInfoVec : interfaceSideInfoPtr->ipInfoVec())
        {
            for (size_t k = 0; k < faceIpInfoVec.size(); ++k)
            {
                ipInfo* ip = faceIpInfoVec[k];

                if (ip->isExposed_)
                    continue;

                // extract current/opposing face/element
                stk::mesh::Entity currentFace = ip->currentFace_;
                stk::mesh::Entity opposingFace = ip->opposingFace_;
                stk::mesh::Entity currentElement = ip->currentElement_;
                stk::mesh::Entity opposingElement = ip->opposingElement_;
                const label currentFaceOrdinal = ip->currentFaceOrdinal_;
                const label opposingFaceOrdinal = ip->opposingFaceOrdinal_;

                // master element; face and volume
                MasterElement* meFCCurrent = ip->meFCCurrent_;
                MasterElement* meFCOpposing = ip->meFCOpposing_;
                MasterElement* meSCSCurrent = ip->meSCSCurrent_;
                MasterElement* meSCSOpposing = ip->meSCSOpposing_;

                // local ip, ordinals, etc
                const label currentGaussPointId = ip->currentGaussPointId_;
                const label opposingGaussPointId = ip->opposingGaussPointId_;

                currentIsoParCoords = ip->currentIsoParCoords_;
                opposingIsoParCoords = ip->opposingIsoParCoords_;

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

                // algorithm related; element
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

                    p_c_Gamma[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        p_c_face_phi[i * current_num_face_nodes + ni] = phi[i];
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

                    p_o_Gamma[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        p_o_face_phi[i * opposing_num_face_nodes + ni] = phi[i];
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
                    connectedNodes[ni] = node;

                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_c_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }

                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);
                    for (label i = 0; i < N; ++i)
                    {
                        p_c_elem_phi[ni * N + i] = phi[i];
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
                    connectedNodes[ni + current_num_elem_nodes] = node;

                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }

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

                // contact surface area fraction
                const scalar fcs = ip->areaFraction_;

                scalar c_amag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar c_axj =
                        c_areaVec[currentGaussPointId * SPATIAL_DIM + j];
                    c_amag += c_axj * c_axj;
                }
                c_amag = std::sqrt(c_amag);

                // compute normal from parent scs area vector
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_cNx[i] =
                        c_areaVec[currentGaussPointId * SPATIAL_DIM + i] /
                        c_amag;
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_oNx[i] = -p_cNx[i];
                }

                // use standard isopar coords for current IP
                // (face-centroid IP location for this gauss point)
                meFCCurrent->general_shape_fcn(1,
                                               &currentIsoParCoords[0],
                                               &ws_c_general_shape_function[0]);

                // compute opposing general shape function
                meFCOpposing->general_shape_fcn(
                    1,
                    &opposingIsoParCoords[0],
                    &ws_o_general_shape_function[0]);

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

                // current inverse length scale
                scalar currentInverseLength = 0.0;
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
                scalar opposingInverseLength = 0.0;
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

                // interpolate face data
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
                    const label offSetDnDx = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];
                        for (label i = 0; i < N; ++i)
                        {
                            currentDiffFluxBip[i] += -currentGammaBip * dndxj *
                                                     nxj *
                                                     p_c_elem_phi[ic * N + i];
                        }
                    }
                }

                // compute diffusion vector; opposing
                for (label ic = 0; ic < opposingNodesPerElement; ++ic)
                {
                    const label offSetDnDx = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_oNx[j];
                        const scalar dndxj = p_o_dndx[offSetDnDx + j];
                        for (label i = 0; i < N; ++i)
                        {
                            opposingDiffFluxBip[i] += -opposingGammaBip *
                                                      dndxj * nxj *
                                                      p_o_elem_phi[ic * N + i];
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

                // extract nearest node
                const label nn = ipNodeMap[currentGaussPointId];

                // mass flow rate at parent IP; advection
                const scalar tmDot =
                    includeAdv ? (stk::mesh::field_data(
                                     *mDotSideSTKFieldPtr_,
                                     currentFace))[currentGaussPointId]
                               : 0.0;
                const scalar abs_tmDot = std::abs(tmDot);

                // compute penalty (DG only)
                const scalar penaltyIp =
                    penaltyFactor *
                    (currentGammaBip * currentInverseLength +
                     opposingGammaBip * opposingInverseLength) /
                    2.0;
                const scalar penaltyMultiplier = 2.0 * (1.0 - multiplier);

                for (label i = 0; i < N; ++i)
                {
                    // diffusive flux
                    const scalar ncDiffFlux =
                        currentDiffFluxBip[i] * multiplier -
                        opposingDiffFluxBip[i] * (1.0 - multiplier);

                    // advective flux: upwind
                    const scalar ncAdv =
                        tmDot * (currentPhiBip[i] + opposingPhiBip[i]) / 2.0 +
                        abs_tmDot * (currentPhiBip[i] - opposingPhiBip[i]) /
                            2.0;

                    // assemble residual weighted by f_cs
                    const label indexR = nn * N + i;
                    const scalar fluxResid =
                        ((ncDiffFlux +
                          penaltyMultiplier * penaltyIp *
                              (currentPhiBip[i] - opposingPhiBip[i])) *
                             fcs * c_amag +
                         fcs * ncAdv);
                    p_rhs[indexR] -= fluxResid;

                    // diagnostics hook. Only meaningful for scalar phi
                    // (N==1).
                    if constexpr (N == 1)
                    {
                        this->recordInterfaceFlux_(
                            currentFace, currentGaussPointId, fluxResid);
                    }

                    // set-up row for matrix
                    const label rowR = indexR * totalNodes * N;

                    // sensitivities; current face (penalty and
                    // advection)
                    const scalar lhsFacC =
                        penaltyMultiplier * penaltyIp * fcs * c_amag +
                        fcs * (abs_tmDot + tmDot) / 2.0;
                    for (label ic = 0; ic < currentNodesPerFace; ++ic)
                    {
                        const label icNdim = c_face_node_ordinals[ic] * N;
                        const scalar r = p_c_general_shape_function[ic];
                        p_lhs[rowR + icNdim + i] += r * lhsFacC;
                    }

                    // sensitivities; current element (diffusion)
                    for (label ic = 0; ic < currentNodesPerElement; ++ic)
                    {
                        const label offSetDnDx = ic * SPATIAL_DIM;
                        const label icNdim = ic * N;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar nxj = p_cNx[j];
                            const scalar dndxj = p_c_dndx[offSetDnDx + j];
                            p_lhs[rowR + icNdim + i] += -currentGammaBip *
                                                        dndxj * nxj * fcs *
                                                        c_amag * multiplier;
                        }
                    }

                    // sensitivities; opposing face (penalty and
                    // advection)
                    const scalar lhsFacO =
                        penaltyMultiplier * penaltyIp * fcs * c_amag +
                        fcs * (abs_tmDot - tmDot) / 2.0;
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
                        const label offSetDnDx = ic * SPATIAL_DIM;
                        const label icNdim = (ic + currentNodesPerElement) * N;
                        for (label j = 0; j < N; ++j)
                        {
                            for (label l = 0; l < SPATIAL_DIM; ++l)
                            {
                                const scalar nxl = p_oNx[l];
                                const scalar dndxl = p_o_dndx[offSetDnDx + l];
                                p_lhs[rowR + icNdim + j] -=
                                    -opposingGammaBip * dndxl * nxl * fcs *
                                    c_amag * (1.0 - multiplier) * rotMat(i, j);
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
