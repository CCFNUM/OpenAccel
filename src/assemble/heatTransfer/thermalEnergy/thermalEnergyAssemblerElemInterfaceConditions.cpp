// File       : thermalEnergyAssemblerElemInterfaceConditions.cpp
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef WITH_THERMAL_TEMPERATURE

#include "interface.h"
#include "ipInfo.h"
#include "simulation.h"
#include "thermalEnergyAssembler.h"

namespace accel
{

void thermalEnergyAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    bool htcTreatment = false;
    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
    {
        if (domain->type() == domainType::fluid)
        {
            if (domain->turbulence_.option_ != turbulenceOption::laminar)
            {
                htcTreatment = true;
            }
        }
        else
        {
            if (interfaceSideInfoPtr->interfPtr()->masterZoneIndex() ==
                domain->index())
            {
                label slaveZoneIndex =
                    interfaceSideInfoPtr->interfPtr()->slaveZoneIndex();
                if (model_->realmRef()
                        .simulationRef()
                        .domainVector()[slaveZoneIndex]
                        ->turbulence_.option_ != turbulenceOption::laminar)
                {
                    htcTreatment = true;
                }
            }
            else
            {
                label masterZoneIndex =
                    interfaceSideInfoPtr->interfPtr()->masterZoneIndex();
                if (model_->realmRef()
                        .simulationRef()
                        .domainVector()[masterZoneIndex]
                        ->turbulence_.option_ != turbulenceOption::laminar)
                {
                    htcTreatment = true;
                }
            }
        }
    }

    if (htcTreatment)
    {
        assembleElemTermsInterfaceSideHTC_(domain, interfaceSideInfoPtr, ctx);
    }
    else
    {
        const bool includeAdv =
            domain->type() == domainType::fluid &&
            !interfaceSideInfoPtr->interfPtr()->isFluidSolidType();

        const auto& mesh = field_broker_->meshRef();

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        scalar penaltyFactor =
            interfaceSideInfoPtr->interfPtr()->penaltyFactor();

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

        // pointers to fixed values
        scalar* p_cNx = &cNx[0];
        scalar* p_oNx = &oNx[0];

        // nodal fields to gather
        std::vector<scalar> ws_c_face_h;
        std::vector<scalar> ws_o_face_h;
        std::vector<scalar> ws_c_elem_h;
        std::vector<scalar> ws_o_elem_h;
        std::vector<scalar> ws_c_face_T;
        std::vector<scalar> ws_o_face_T;
        std::vector<scalar> ws_c_elem_T;
        std::vector<scalar> ws_o_elem_T;
        std::vector<scalar> ws_c_elem_coordinates;
        std::vector<scalar> ws_o_elem_coordinates;
        std::vector<scalar> ws_c_lambdaEff;
        std::vector<scalar> ws_o_lambdaEff;
        std::vector<scalar> ws_c_cp;
        std::vector<scalar> ws_o_cp;

        // master element data
        std::vector<scalar> ws_c_dndx;
        std::vector<scalar> ws_o_dndx;
        std::vector<scalar> ws_c_det_j;
        std::vector<scalar> ws_o_det_j;
        std::vector<scalar> ws_c_general_shape_function;
        std::vector<scalar> ws_o_general_shape_function;

        // Get transport fields/side fields
        const auto& hSTKFieldRef = phi_->stkFieldRef();
        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
        auto& qDotSideSTKFieldRef =
            model_->qDotRef().sideFieldRef().stkFieldRef();

        // Get geometric fields
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // extract vector of interface IP info
        const scalar multiplier = 0.5;

        const auto& ipInfoVec = interfaceSideInfoPtr->ipInfoVec();

        for (label iSide = 0; iSide < static_cast<label>(ipInfoVec.size());
             iSide++)
        {
            const auto& faceIpInfoVec = ipInfoVec[iSide];

            // Local cs-fragment index within this face

            for (size_t k = 0; k < faceIpInfoVec.size(); ++k)
            {
                const ipInfo* ip = faceIpInfoVec[k];

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
                const label lhsSize = totalNodes * totalNodes;
                const label rhsSize = totalNodes;
                lhs.resize(lhsSize);
                rhs.resize(rhsSize);
                scratchIds.resize(rhsSize);
                scratchVals.resize(rhsSize);
                connectedNodes.resize(totalNodes);

                // algorithm related; element; dndx will be at a
                // single gauss point...
                ws_c_elem_h.resize(currentNodesPerElement);
                ws_o_elem_h.resize(opposingNodesPerElement);
                ws_c_elem_T.resize(currentNodesPerElement);
                ws_o_elem_T.resize(opposingNodesPerElement);
                ws_c_elem_coordinates.resize(currentNodesPerElement *
                                             SPATIAL_DIM);
                ws_o_elem_coordinates.resize(opposingNodesPerElement *
                                             SPATIAL_DIM);
                ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
                ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
                ws_c_det_j.resize(1);
                ws_o_det_j.resize(1);

                // algorithm related; face
                ws_c_face_h.resize(currentNodesPerFace);
                ws_o_face_h.resize(opposingNodesPerFace);
                ws_c_face_T.resize(currentNodesPerFace);
                ws_o_face_T.resize(opposingNodesPerFace);
                ws_c_lambdaEff.resize(currentNodesPerFace);
                ws_o_lambdaEff.resize(opposingNodesPerFace);
                ws_c_cp.resize(currentNodesPerFace);
                ws_o_cp.resize(opposingNodesPerFace);
                ws_c_general_shape_function.resize(currentNodesPerFace);
                ws_o_general_shape_function.resize(opposingNodesPerFace);

                // pointers
                scalar* p_lhs = &lhs[0];
                scalar* p_rhs = &rhs[0];

                scalar* p_c_face_h = &ws_c_face_h[0];
                scalar* p_o_face_h = &ws_o_face_h[0];
                scalar* p_c_elem_h = &ws_c_elem_h[0];
                scalar* p_o_elem_h = &ws_o_elem_h[0];
                scalar* p_c_face_T = &ws_c_face_T[0];
                scalar* p_o_face_T = &ws_o_face_T[0];
                scalar* p_c_elem_T = &ws_c_elem_T[0];
                scalar* p_o_elem_T = &ws_o_elem_T[0];
                scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
                scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];
                scalar* p_c_lambdaEff = &ws_c_lambdaEff[0];
                scalar* p_o_lambdaEff = &ws_o_lambdaEff[0];
                scalar* p_c_cp = &ws_c_cp[0];
                scalar* p_o_cp = &ws_o_cp[0];

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
                    p_c_face_h[ni] = *stk::mesh::field_data(hSTKFieldRef, node);
                    p_c_face_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                    p_c_lambdaEff[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                    p_c_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                    p_c_lambdaEff[ni] *= p_c_cp[ni];
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
                    p_o_face_h[ni] = *stk::mesh::field_data(hSTKFieldRef, node);
                    p_o_face_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                    p_o_lambdaEff[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                    p_o_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                    p_o_lambdaEff[ni] *= p_o_cp[ni];
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

                    // gather scalar
                    p_c_elem_h[ni] = *stk::mesh::field_data(hSTKFieldRef, node);
                    p_c_elem_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

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

                    // gather scalar
                    p_o_elem_h[ni] = *stk::mesh::field_data(hSTKFieldRef, node);
                    p_o_elem_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                    // gather; vector
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_o_elem_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }
                }

                // apply transformations
                interfaceSideInfoPtr->transformCoordinateList(
                    ws_o_elem_coordinates, opposing_num_elem_nodes);

                // pointer to face data
                scalar* qDot =
                    stk::mesh::field_data(qDotSideSTKFieldRef, currentFace);
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

                // compute opposing normal: in theory it is assumed
                // that the current and opposing sub-control surfaces are
                // sufficiently planar
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
                scalar currentHBip = 0.0;
                meFCCurrent->interpolatePoint(
                    1, &currentIsoParCoords[0], &ws_c_face_h[0], &currentHBip);

                scalar opposingHBip = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_o_face_h[0],
                                               &opposingHBip);

                scalar currentTBip = 0.0;
                meFCCurrent->interpolatePoint(
                    1, &currentIsoParCoords[0], &ws_c_face_T[0], &currentTBip);

                scalar opposingTBip = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_o_face_T[0],
                                               &opposingTBip);

                scalar currentLambdaEffBip = 0.0;
                meFCCurrent->interpolatePoint(1,
                                              &currentIsoParCoords[0],
                                              &ws_c_lambdaEff[0],
                                              &currentLambdaEffBip);

                scalar opposingLambdaEffBip = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_o_lambdaEff[0],
                                               &opposingLambdaEffBip);

                scalar currentCpBip = 0.0;
                meFCCurrent->interpolatePoint(
                    1, &currentIsoParCoords[0], &ws_c_cp[0], &currentCpBip);

                scalar opposingCpBip = 0.0;
                meFCOpposing->interpolatePoint(
                    1, &opposingIsoParCoords[0], &ws_o_cp[0], &opposingCpBip);

                // GGI scatter (single-sided)

                // compute diffusion vector; current
                scalar currentDiffFluxBip = 0;
                for (label ic = 0; ic < currentNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];

                        const scalar T = p_c_elem_T[ic];

                        // -Gamma*dphi/dxj*Aj
                        currentDiffFluxBip +=
                            -currentLambdaEffBip * dndxj * nxj * T;
                    }
                }

                // compute diffusion vector; opposing
                scalar opposingDiffFluxBip = 0;
                for (label ic = 0; ic < opposingNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_oNx[j];
                        const scalar dndxj = p_o_dndx[offSetDnDx + j];

                        const scalar T = p_o_elem_T[ic];

                        // -Gamma*dphi/dxj*Aj
                        opposingDiffFluxBip +=
                            -opposingLambdaEffBip * dndxj * nxj * T;
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

                // compute penalty (DG only)
                const scalar penaltyIp =
                    penaltyFactor *
                    (currentLambdaEffBip * currentInverseLength +
                     opposingLambdaEffBip * opposingInverseLength) /
                    2.0;
                const scalar penaltyMultiplier = 2.0 * (1.0 - multiplier);

                // non conformal diffusive flux
                const scalar ncDiffFlux =
                    currentDiffFluxBip * multiplier -
                    opposingDiffFluxBip * (1.0 - multiplier);

                // non conformal advection
                const scalar ncAdv =
                    tmDot * (currentHBip + opposingHBip) / 2.0 +
                    abs_tmDot * (currentHBip - opposingHBip) / 2.0;

                // assemble residual; form proper rhs index for
                // current face assembly
                const label indexR = nn;
                p_rhs[indexR] -=
                    ((ncDiffFlux + penaltyIp * penaltyMultiplier *
                                       (currentTBip - opposingTBip)) *
                         fcs * c_amag +
                     fcs * ncAdv);

                // fill the nc-heat flow rate; accumulate so that
                // multiple ggiInfos sharing one master gauss point
                // sum into a single IP slot (non-conformal mesh
                // can have multiple slave SCSes overlapping one
                // master SCS). Zeroing happens up-front in
                // assembleElemTermsInterfaces_.
                qDot[currentGaussPointId] +=
                    ((ncDiffFlux + penaltyIp * penaltyMultiplier *
                                       (currentTBip - opposingTBip)) *
                         fcs * c_amag +
                     fcs * ncAdv);

                // set-up row for matrix
                const label rowR = indexR * totalNodes;

                // sensitivities; current face (penalty and
                // advection); use general shape function for
                // this single ip
                const scalar lhsFacC = penaltyMultiplier * penaltyIp /
                                           currentCpBip * fcs * c_amag +
                                       fcs * (abs_tmDot + tmDot) / 2.0;
                for (label ic = 0; ic < currentNodesPerFace; ++ic)
                {
                    const label icNdim = c_face_node_ordinals[ic];
                    const scalar r = p_c_general_shape_function[ic];
                    p_lhs[rowR + icNdim] += r * lhsFacC;
                }

                // sensitivities; current element (diffusion)
                for (label ic = 0; ic < currentNodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point
                    const label icNdim = ic;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nxj = p_cNx[j];
                        const scalar dndxj = p_c_dndx[offSetDnDx + j];

                        // -Gamma*dphi/dxj*nj*dS
                        p_lhs[rowR + icNdim] += -currentLambdaEffBip * dndxj *
                                                nxj * fcs * c_amag *
                                                multiplier / currentCpBip;
                    }
                }

                // sensitivities; opposing face (penalty and
                // advection)
                const scalar lhsFacO = penaltyMultiplier * penaltyIp /
                                           opposingCpBip * fcs * c_amag +
                                       fcs * (abs_tmDot - tmDot) / 2.0;
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
                    const label offSetDnDx =
                        ic * SPATIAL_DIM; // single intg. point
                    const label icNdim = (ic + currentNodesPerElement);
                    for (label l = 0; l < SPATIAL_DIM; ++l)
                    {
                        const scalar nxl = p_oNx[l];
                        const scalar dndxl = p_o_dndx[offSetDnDx + l];

                        // -Gamma*dphi/dxj*nj*dS
                        p_lhs[rowR + icNdim] -=
                            -opposingLambdaEffBip * dndxl * nxl * fcs * c_amag *
                            (1.0 - multiplier) / opposingCpBip;
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void thermalEnergyAssembler::assembleElemTermsInterfaceSideHTC_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

    // mapping for -1:1 -> -0.5:0.5 volume element
    std::vector<scalar> currentElementIsoParCoords(SPATIAL_DIM);

    // nodal fields to gather
    std::vector<scalar> ws_c_face_T;
    std::vector<scalar> ws_o_face_T;
    std::vector<scalar> ws_c_cp;
    std::vector<scalar> ws_o_cp;

    // master element data
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // Get transport fields/side fields
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* TWallCoeffsSTKFieldPtr = model_->TWallCoeffsRef().stkFieldPtr();
    auto& qDotSideSTKFieldRef = model_->qDotRef().sideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // GGI path is not yet validated for the total-energy HTC problem;
    // preserve the original errorMsg behavior at the top of the function and
    // run the unified loop for DG below.

    // Unified loop over per-IP info.  Storage is owned by the base
    // interfaceSideInfo; concrete side classes store derived records upcast to
    // ipInfo*, and ip->areaFraction_ is the default 1.0 so the math is
    // identical to the original DG implementation.
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
            MasterElement* meSCSCurrent = ip->meSCSCurrent_;
            MasterElement* meFCOpposing = ip->meFCOpposing_;
            MasterElement* meSCSOpposing = ip->meSCSOpposing_;

            // local ip, ordinals, etc
            const label currentGaussPointId = ip->currentGaussPointId_;
            currentIsoParCoords = ip->currentIsoParCoords_;
            opposingIsoParCoords = ip->opposingIsoParCoords_;

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap =
                meSCSCurrent->ipNodeMap(currentFaceOrdinal);

            // extract some master element info
            const label currentNodesPerFace = meFCCurrent->nodesPerElement_;
            const label currentNodesPerElement = meSCSCurrent->nodesPerElement_;
            const label opposingNodesPerFace = meFCOpposing->nodesPerElement_;
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

            // algorithm related; face
            ws_c_face_T.resize(currentNodesPerFace);
            ws_o_face_T.resize(opposingNodesPerFace);
            ws_c_cp.resize(currentNodesPerFace);
            ws_o_cp.resize(opposingNodesPerFace);
            ws_c_general_shape_function.resize(currentNodesPerFace);
            ws_o_general_shape_function.resize(opposingNodesPerFace);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];

            scalar* p_c_face_T = &ws_c_face_T[0];
            scalar* p_o_face_T = &ws_o_face_T[0];
            scalar* p_c_cp = &ws_c_cp[0];
            scalar* p_o_cp = &ws_o_cp[0];

            // me pointers
            scalar* p_c_general_shape_function =
                &ws_c_general_shape_function[0];
            scalar* p_o_general_shape_function =
                &ws_o_general_shape_function[0];

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
                p_c_face_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                p_c_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
            }

            // gather opposing face data
            stk::mesh::Entity const* opposing_face_node_rels =
                bulkData.begin_nodes(opposingFace);
            const label opposing_num_face_nodes =
                bulkData.num_nodes(opposingFace);
            for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
            {
                stk::mesh::Entity node = opposing_face_node_rels[ni];

                // gather; scalar
                p_o_face_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                p_o_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
            }

            // populate opposing face_node_ordinals
            const label* o_face_node_ordinals =
                meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

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
            }

            // pointer to face data
            scalar* qDot =
                stk::mesh::field_data(qDotSideSTKFieldRef, currentFace);
            const scalar* c_TWallCoeff =
                stk::mesh::field_data(*TWallCoeffsSTKFieldPtr, currentFace);
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

            // project from side to element; method deals with
            // the -1:1 isInElement range to the proper
            // underlying CVFEM range
            meSCSCurrent->sidePcoords_to_elemPcoords(
                currentFaceOrdinal,
                1,
                &currentIsoParCoords[0],
                &currentElementIsoParCoords[0]);

            // interpolate face data; current and opposing...
            scalar currentTBip = 0.0;
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_face_T[0], &currentTBip);
            scalar opposingTBip = 0.0;
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_face_T[0], &opposingTBip);

            scalar currentCpBip = 0.0;
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_cp[0], &currentCpBip);
            scalar opposingCpBip = 0.0;
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_cp[0], &opposingCpBip);

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

            // Flux = HTC * A * Delta T
            scalar diffFlux = -c_TWallCoeff[currentGaussPointId] * c_amag *
                              (opposingTBip - currentTBip);

            // assemble residual; form proper rhs index for
            // current face assembly
            const label indexR = nn;
            p_rhs[indexR] -= diffFlux;

            // fill the nc-heat flow rate (accumulate to handle
            // multiple ggiInfos per gauss point)
            qDot[currentGaussPointId] += diffFlux;

            // set-up row for matrix
            const label rowR = indexR * totalNodes;

            // sensitivities; current face
            const scalar lhsFacC =
                c_TWallCoeff[currentGaussPointId] * c_amag / currentCpBip;
            for (label ic = 0; ic < currentNodesPerFace; ++ic)
            {
                const label icNdim = c_face_node_ordinals[ic];
                const scalar r = p_c_general_shape_function[ic];
                p_lhs[rowR + icNdim] += r * lhsFacC;
            }

            // sensitivities; opposing face
            const scalar lhsFacO =
                c_TWallCoeff[currentGaussPointId] * c_amag / opposingCpBip;
            for (label ic = 0; ic < opposingNodesPerFace; ++ic)
            {
                const label icNdim =
                    (o_face_node_ordinals[ic] + currentNodesPerElement);
                const scalar r = p_o_general_shape_function[ic];
                p_lhs[rowR + icNdim] -= r * lhsFacO;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif /* WITH_THERMAL_TEMPERATURE */
