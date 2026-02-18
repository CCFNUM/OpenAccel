// File       : thermalEnergyAssemblerElemInterfaceConditions.cpp
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WITH_THERMAL_TEMPERATURE
#ifdef HAS_INTERFACE

#include "dgInfo.h"
#include "interface.h"
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
    bool hybridTreatment = false;
    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
    {
        if (domain->type() == domainType::fluid)
        {
            if (domain->turbulence_.transitional_)
            {
                hybridTreatment = true;
            }
            else if (domain->turbulence_.option_ != turbulenceOption::laminar)
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
                        ->turbulence_.transitional_)
                {
                    hybridTreatment = true;
                }
                else if (model_->realmRef()
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
                        ->turbulence_.transitional_)
                {
                    hybridTreatment = true;
                }
                else if (model_->realmRef()
                             .simulationRef()
                             .domainVector()[masterZoneIndex]
                             ->turbulence_.option_ != turbulenceOption::laminar)
                {
                    htcTreatment = true;
                }
            }
        }
    }

    if (hybridTreatment)
    {
#if 0
        assembleElemTermsInterfaceSideHybrid_(
            domain, interfaceSideInfoPtr, ctx);
#else
        assembleElemTermsInterfaceSideHTC_(domain, interfaceSideInfoPtr, ctx);
#endif
    }
    else if (htcTreatment)
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

        // extract vector of dgInfo
        const std::vector<std::vector<dgInfo*>>& dgInfoVec =
            interfaceSideInfoPtr->dgInfoVec_;

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
                interfaceSideInfoPtr->rotateVectorListCompact<1>(
                    ws_o_face_h, opposing_num_face_nodes);
                interfaceSideInfoPtr->rotateVectorList<1>(
                    ws_o_elem_h, opposing_num_elem_nodes);
                interfaceSideInfoPtr->transformCoordinateList(
                    ws_o_elem_coordinates, opposing_num_elem_nodes);

                // pointer to face data
                scalar* qDot =
                    stk::mesh::field_data(qDotSideSTKFieldRef, currentFace);
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
                    (currentLambdaEffBip * currentInverseLength +
                     opposingLambdaEffBip * opposingInverseLength) /
                    2.0;

                // non conformal diffusive flux
                const scalar ncDiffFlux =
                    (currentDiffFluxBip - opposingDiffFluxBip) / 2.0;

                // non conformal advection
                const scalar ncAdv =
                    tmDot * (currentHBip + opposingHBip) / 2.0 +
                    abs_tmDot * (currentHBip - opposingHBip) / 2.0;

                // assemble residual; form proper rhs index for
                // current face assembly
                const label indexR = nn;
                p_rhs[indexR] -=
                    ((ncDiffFlux + penaltyIp * (currentTBip - opposingTBip)) *
                         c_amag +
                     ncAdv);

                // fill the nc-heat flow rate
                qDot[currentGaussPointId] =
                    ((ncDiffFlux + penaltyIp * (currentTBip - opposingTBip)) *
                         c_amag +
                     ncAdv);

                // set-up row for matrix
                const label rowR = indexR * totalNodes;

                // sensitivities; current face (penalty and
                // advection); use general shape function for
                // this single ip
                const scalar lhsFacC = penaltyIp / currentCpBip * c_amag +
                                       (abs_tmDot + tmDot) / 2.0;
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
                                                nxj * c_amag / 2.0 /
                                                currentCpBip;
                    }
                }

                // sensitivities; opposing face (penalty and
                // advection); use general shape function for
                // this single ip
                const scalar lhsFacO = penaltyIp / opposingCpBip * c_amag +
                                       (abs_tmDot - tmDot) / 2.0;
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
                        p_lhs[rowR + icNdim] -= -opposingLambdaEffBip * dndxl *
                                                nxl * c_amag / 2.0 /
                                                opposingCpBip;
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void thermalEnergyAssembler::assembleElemTermsInterfaceSideHybrid_(
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
    std::vector<scalar> ws_c_face_gamma;
    std::vector<scalar> ws_o_face_gamma;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;
    std::vector<scalar> ws_c_lambdaEff;
    std::vector<scalar> ws_o_lambdaEff;
    std::vector<scalar> ws_c_cp;
    std::vector<scalar> ws_o_cp;
    std::vector<scalar> ws_c_face_U;
    std::vector<scalar> ws_o_face_U;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // workspace for velocity interpolation
    std::vector<scalar> UBip(SPATIAL_DIM);

    // Get transport fields/side fields
    const auto& hSTKFieldRef = phi_->stkFieldRef();
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto* TWallCoeffsSTKFieldPtr = model_->TWallCoeffsRef().stkFieldPtr();
    const auto& gammaSTKFieldRef = model_->gammaRef().stkFieldRef();
    auto& qDotSideSTKFieldRef = model_->qDotRef().sideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // extract vector of dgInfo
    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
        interfaceSideInfoPtr->dgInfoVec_;

    // Compressibility flag: 1.0 for compressible, 0.0 otherwise
    // Used to avoid branching in hot loop
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;
    const bool currentIsFluid = (domain->type() == domainType::fluid);
    const scalar currentFluidFlag = currentIsFluid ? 1.0 : 0.0;

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
            MasterElement* meSCSCurrent = dgInfo->meSCSCurrent_;
            MasterElement* meFCOpposing = dgInfo->meFCOpposing_;
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
            ws_c_elem_h.resize(currentNodesPerElement);
            ws_o_elem_h.resize(opposingNodesPerElement);
            ws_c_elem_T.resize(currentNodesPerElement);
            ws_o_elem_T.resize(opposingNodesPerElement);
            ws_c_elem_coordinates.resize(currentNodesPerElement * SPATIAL_DIM);
            ws_o_elem_coordinates.resize(opposingNodesPerElement * SPATIAL_DIM);
            ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
            ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
            ws_c_det_j.resize(1);
            ws_o_det_j.resize(1);

            // algorithm related; face
            ws_c_face_h.resize(currentNodesPerFace);
            ws_o_face_h.resize(opposingNodesPerFace);
            ws_c_face_T.resize(currentNodesPerFace);
            ws_o_face_T.resize(opposingNodesPerFace);
            ws_c_face_gamma.resize(currentNodesPerFace);
            ws_o_face_gamma.resize(opposingNodesPerFace);
            ws_c_lambdaEff.resize(currentNodesPerFace);
            ws_o_lambdaEff.resize(opposingNodesPerFace);
            ws_c_cp.resize(currentNodesPerFace);
            ws_o_cp.resize(opposingNodesPerFace);
            ws_c_face_U.resize(currentNodesPerFace * SPATIAL_DIM);
            ws_o_face_U.resize(opposingNodesPerFace * SPATIAL_DIM);
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
            scalar* p_c_face_gamma = &ws_c_face_gamma[0];
            scalar* p_o_face_gamma = &ws_o_face_gamma[0];
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
                p_c_face_gamma[ni] =
                    *stk::mesh::field_data(gammaSTKFieldRef, node);
                p_c_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_c_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_c_lambdaEff[ni] *= p_c_cp[ni];

                // gather; vector (velocity)
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_c_face_U[ni * SPATIAL_DIM + j] = U[j];
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
                p_o_face_h[ni] = *stk::mesh::field_data(hSTKFieldRef, node);
                p_o_face_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                p_o_face_gamma[ni] =
                    *stk::mesh::field_data(gammaSTKFieldRef, node);
                p_o_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_o_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_o_lambdaEff[ni] *= p_o_cp[ni];

                // gather; vector (velocity)
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_o_face_U[ni * SPATIAL_DIM + j] = U[j];
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
            interfaceSideInfoPtr->rotateVectorListCompact<1>(
                ws_o_face_h, opposing_num_face_nodes);
            interfaceSideInfoPtr->rotateVectorList<1>(ws_o_elem_h,
                                                      opposing_num_elem_nodes);
            interfaceSideInfoPtr->transformCoordinateList(
                ws_o_elem_coordinates, opposing_num_elem_nodes);

            // pointer to face data
            scalar* qDot =
                stk::mesh::field_data(qDotSideSTKFieldRef, currentFace);
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
                &opposingElementIsoParCoordsTrans[0]);

            // compute dndx
            scalar scs_error = 0.0;
            meSCSCurrent->general_face_grad_op(currentFaceOrdinal,
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
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_face_h[0], &opposingHBip);

            scalar currentTBip = 0.0;
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_face_T[0], &currentTBip);

            scalar opposingTBip = 0.0;
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_face_T[0], &opposingTBip);

            scalar gammaBip = 0.0;
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_face_gamma[0], &gammaBip);

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

            // compute diffusion vector; current
            scalar currentDiffFluxBip = 0;
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point

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
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point

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

            // extract nearset node
            const label nn = ipNodeMap[currentGaussPointId];

            // compute general shape function at current and
            // opposing integration points
            meFCCurrent->general_shape_fcn(
                1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
            meFCOpposing->general_shape_fcn(
                1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);

            // LAM PART

            // compute penalty
            const scalar penaltyIp =
                penaltyFactor *
                (currentLambdaEffBip * currentInverseLength +
                 opposingLambdaEffBip * opposingInverseLength) /
                2.0;

            // non conformal diffusive flux
            const scalar ncDiffFlux =
                (currentDiffFluxBip - opposingDiffFluxBip) / 2.0;

            const scalar lamFlux =
                ((ncDiffFlux + penaltyIp * (currentTBip - opposingTBip)) *
                 c_amag);

            // HTC PART
            // Since we solve for enthalpy in BOTH fluid and solid domains,
            // both sides need Cp division: q = HTC * (h_solid/cp_solid -
            // h_fluid/cp_fluid)
            const scalar currentCpFactor = currentCpBip + SMALL;
            const scalar opposingCpFactor = opposingCpBip + SMALL;

            // Compute dynamic temperature T_dyn = u²/(2*Cp) for compressible
            // flow Velocity is zero in solid, so we interpolate from both sides
            // and the solid contribution vanishes naturally (branch-free)

            // Interpolate velocity from current side (zero if solid)
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_face_U[0], &UBip[0]);
            scalar uMagSq = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                uMagSq += UBip[j] * UBip[j];
            }

            // Interpolate velocity from opposing side (zero if solid) and add
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_face_U[0], &UBip[0]);
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                uMagSq += UBip[j] * UBip[j];
            }

            // Select fluid Cp using multiplicative flag (branch-free)
            const scalar fluidCp = currentFluidFlag * currentCpBip +
                                   (1.0 - currentFluidFlag) * opposingCpBip;

            // Compute T_dyn (comp=0 for incompressible, making T_dyn vanish.
            // Will be 0 also for a solid domain)
            const scalar T_dyn = comp * uMagSq / (2.0 * (fluidCp + SMALL));

            const scalar* c_TWallCoeff =
                stk::mesh::field_data(*TWallCoeffsSTKFieldPtr, currentFace);

            // Flux = HTC * A * (Delta T - T_dyn)
            const scalar htcFlux = -c_TWallCoeff[currentGaussPointId] * c_amag *
                                   (opposingTBip - currentTBip - T_dyn);

            // BLEND
            const scalar finalFlux =
                gammaBip * htcFlux + (1.0 - gammaBip) * lamFlux;

            // assemble residual; form proper rhs index for
            // current face assembly
            const label indexR = nn;
            p_rhs[indexR] -= finalFlux;

            // fill the nc-heat flow rate
            qDot[currentGaussPointId] = finalFlux;

            // set-up row for matrix
            const label rowR = indexR * totalNodes;

            // sensitivities; current face (penalty and
            // advection); use general shape function for
            // this single ip
            const scalar lhsFacC = penaltyIp / currentCpBip * c_amag;

            const scalar htcLhsFacC =
                c_TWallCoeff[currentGaussPointId] * c_amag / currentCpFactor;

            for (label ic = 0; ic < currentNodesPerFace; ++ic)
            {
                const label icNdim = c_face_node_ordinals[ic];
                const scalar r = p_c_general_shape_function[ic];
                p_lhs[rowR + icNdim] +=
                    gammaBip * r * htcLhsFacC + (1.0 - gammaBip) * r * lhsFacC;
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

                    // -Gamma*dphi/dxj*nj*dS
                    p_lhs[rowR + icNdim] +=
                        (1.0 - gammaBip) * (-currentLambdaEffBip * dndxj * nxj *
                                            c_amag / 2.0 / currentCpBip);
                }
            }

            // sensitivities; opposing face (penalty and
            // advection); use general shape function for
            // this single ip
            const scalar lhsFacO = penaltyIp / opposingCpBip * c_amag;

            const scalar htcLhsFacO =
                c_TWallCoeff[currentGaussPointId] * c_amag / opposingCpFactor;

            for (label ic = 0; ic < opposingNodesPerFace; ++ic)
            {
                const label icNdim =
                    (o_face_node_ordinals[ic] + currentNodesPerElement);
                const scalar r = p_o_general_shape_function[ic];
                p_lhs[rowR + icNdim] -=
                    gammaBip * r * htcLhsFacO + (1.0 - gammaBip) * r * lhsFacO;
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

                    // -Gamma*dphi/dxj*nj*dS
                    p_lhs[rowR + icNdim] -=
                        (1.0 - gammaBip) * (-opposingLambdaEffBip * dndxl *
                                            nxl * c_amag / 2.0 / opposingCpBip);
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
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
    std::vector<scalar> ws_c_face_U;
    std::vector<scalar> ws_o_face_U;

    // master element data
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // workspace for velocity interpolation
    std::vector<scalar> UBip(SPATIAL_DIM);

    // Get transport fields/side fields
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto* TWallCoeffsSTKFieldPtr = model_->TWallCoeffsRef().stkFieldPtr();
    auto& qDotSideSTKFieldRef = model_->qDotRef().sideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // extract vector of dgInfo
    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
        interfaceSideInfoPtr->dgInfoVec_;

    // Compressibility flag: 1.0 for compressible, 0.0 otherwise
    // Used to avoid branching in hot loop
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;
    const bool currentIsFluid = (domain->type() == domainType::fluid);
    const scalar currentFluidFlag = currentIsFluid ? 1.0 : 0.0;

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
            MasterElement* meSCSCurrent = dgInfo->meSCSCurrent_;
            MasterElement* meFCOpposing = dgInfo->meFCOpposing_;
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
            ws_c_face_U.resize(currentNodesPerFace * SPATIAL_DIM);
            ws_o_face_U.resize(opposingNodesPerFace * SPATIAL_DIM);
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

                // gather; vector (velocity)
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_c_face_U[ni * SPATIAL_DIM + j] = U[j];
                }
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

                // gather; vector (velocity)
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_o_face_U[ni * SPATIAL_DIM + j] = U[j];
                }
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

            // interpolate Cp at integration points
            scalar currentCpBip = 0.0;
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_cp[0], &currentCpBip);
            scalar opposingCpBip = 0.0;
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_cp[0], &opposingCpBip);

            // Since we solve for enthalpy in BOTH fluid and solid domains,
            // both sides need Cp division: q = HTC * (h_solid/cp_solid -
            // h_fluid/cp_fluid)
            const scalar currentCpFactor = currentCpBip + SMALL;
            const scalar opposingCpFactor = opposingCpBip + SMALL;

            // Compute dynamic temperature T_dyn = u²/(2*Cp) for compressible
            // flow Velocity is zero in solid, so we interpolate from both sides
            // and the solid contribution vanishes naturally (branch-free)

            // Interpolate velocity from current side (zero if solid)
            meFCCurrent->interpolatePoint(
                1, &currentIsoParCoords[0], &ws_c_face_U[0], &UBip[0]);
            scalar uMagSq = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                uMagSq += UBip[j] * UBip[j];
            }

            // Interpolate velocity from opposing side (zero if solid) and add
            meFCOpposing->interpolatePoint(
                1, &opposingIsoParCoords[0], &ws_o_face_U[0], &UBip[0]);
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                uMagSq += UBip[j] * UBip[j];
            }

            // Select fluid Cp using multiplicative flag (branch-free)
            const scalar fluidCp = currentFluidFlag * currentCpBip +
                                   (1.0 - currentFluidFlag) * opposingCpBip;

            // Compute T_dyn (comp=0 for incompressible, making T_dyn vanish.
            // Will be 0 also for a solid domain)
            const scalar T_dyn = comp * uMagSq / (2.0 * (fluidCp + SMALL));

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

            // Flux = HTC * A * (Delta T - T_dyn)
            scalar diffFlux = -c_TWallCoeff[currentGaussPointId] * c_amag *
                              (opposingTBip - currentTBip - T_dyn);

            // assemble residual; form proper rhs index for
            // current face assembly
            const label indexR = nn;
            p_rhs[indexR] -= diffFlux;

            // fill the nc-heat flow rate
            qDot[currentGaussPointId] = diffFlux;

            // set-up row for matrix
            const label rowR = indexR * totalNodes;

            // sensitivities; current face
            // For fluid side (enthalpy): divide by Cp
            // For solid side (temperature): no Cp division
            const scalar lhsFacC =
                c_TWallCoeff[currentGaussPointId] * c_amag / currentCpFactor;
            for (label ic = 0; ic < currentNodesPerFace; ++ic)
            {
                const label icNdim = c_face_node_ordinals[ic];
                const scalar r = p_c_general_shape_function[ic];
                p_lhs[rowR + icNdim] += r * lhsFacC;
            }

            // sensitivities; opposing face
            // For fluid side (enthalpy): divide by Cp
            // For solid side (temperature): no Cp division
            const scalar lhsFacO =
                c_TWallCoeff[currentGaussPointId] * c_amag / opposingCpFactor;
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

#endif /* HAS_INTERFACE */
#endif /* WITH_THERMAL_TEMPERATURE */
