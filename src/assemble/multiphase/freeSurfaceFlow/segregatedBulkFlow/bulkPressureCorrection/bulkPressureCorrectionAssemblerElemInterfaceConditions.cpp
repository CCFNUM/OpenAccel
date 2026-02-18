// File       : bulkPressureCorrectionAssemblerElemInterfaceConditions.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

// code
#include "bulkPressureCorrectionAssembler.h"
#include "dgInfo.h"
#include "freeSurfaceFlowModel.h"
#include "interface.h"
#include "mesh.h"
#include "navierStokesAssembler.h"
#include "zoneTransformation.h"

namespace accel
{

void bulkPressureCorrectionAssembler::assembleElemTermsInterfaces_(
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

void bulkPressureCorrectionAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    // required for frame motion (MFR)
    const auto coriolisMatrix =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* p_mat = coriolisMatrix.data();

    const auto origin =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().origin_
            : utils::vector::Zero();
    const scalar* p_ori = origin.data();

    // mesh motion
    const bool meshMoving = domain->zonePtr()->meshMoving();

    const bool compressible = domain->isMaterialCompressible();
    const scalar comp = compressible ? 1.0 : 0.0;

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

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
    std::vector<scalar> cRhoUBip(SPATIAL_DIM);
    std::vector<scalar> oRhoUBip(SPATIAL_DIM);
    std::vector<scalar> cUmBip(SPATIAL_DIM);
    std::vector<scalar> cDuBip(SPATIAL_DIM);
    std::vector<scalar> oDuBip(SPATIAL_DIM);
    std::vector<scalar> currentCoordsBip(SPATIAL_DIM);

    // pressure stabilization
    std::vector<scalar> cGjpBip(SPATIAL_DIM);
    std::vector<scalar> oGjpBip(SPATIAL_DIM);
    std::vector<scalar> cDpdxBip(SPATIAL_DIM);
    std::vector<scalar> oDpdxBip(SPATIAL_DIM);

    // mapping for -1:1 -> -0.5:0.5 volume element
    std::vector<scalar> currentElementIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingElementIsoParCoords(SPATIAL_DIM);

    // interpolate nodal values to point-in-elem
    const label sizeOfScalarField = 1;
    const label sizeOfVectorField = SPATIAL_DIM;

    // pointers to fixed values
    scalar* p_cNx = &cNx[0];
    scalar* p_oNx = &oNx[0];

    // nodal fields to gather; face
    std::vector<scalar> ws_c_p;
    std::vector<scalar> ws_o_p;
    std::vector<scalar> ws_c_Gjp;
    std::vector<scalar> ws_o_Gjp;
    std::vector<scalar> ws_c_rhoU;
    std::vector<scalar> ws_o_rhoU;
    std::vector<scalar> ws_c_Um;
    std::vector<scalar> ws_c_rho;
    std::vector<scalar> ws_o_rho;
    std::vector<scalar> ws_c_alpha;
    std::vector<scalar> ws_o_alpha;
    std::vector<scalar> ws_c_psi;
    std::vector<scalar> ws_o_psi;
    std::vector<scalar> ws_c_face_coordinates;

    // element
    std::vector<scalar> ws_c_elem_p;
    std::vector<scalar> ws_o_elem_p;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;
    std::vector<scalar> ws_c_du;
    std::vector<scalar> ws_o_du;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;
    std::vector<scalar> ws_c_general_shape_function;
    std::vector<scalar> ws_o_general_shape_function;

    // Get transport fields/side fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& gradPSTKFieldRef = model_->pRef().gradRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef(phaseIndex_).sideFieldRef().stkFieldRef();

    const auto* psiSTKFieldPtr =
        compressible ? model_->psiRef(phaseIndex_).stkFieldPtr() : nullptr;

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get pressure diffusivity coefficient field
    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    scalar densityScale = model_->rhoRef(phaseIndex_).scale();

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

            // if gauss point is exposed (non-overlapping), then
            // treat as a wall
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
            const label currentNodesPerSide = meFCCurrent->nodesPerElement_;
            const label opposingNodesPerSide = meFCOpposing->nodesPerElement_;
            const label currentNodesPerElement = meSCSCurrent->nodesPerElement_;
            const label opposingNodesPerElement =
                meSCSOpposing->nodesPerElement_;

            // arithmetic weights
            const scalar f_c = 1.0 / static_cast<scalar>(currentNodesPerSide);
            const scalar f_o = 1.0 / static_cast<scalar>(opposingNodesPerSide);

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
            ws_c_p.resize(currentNodesPerSide);
            ws_o_p.resize(opposingNodesPerSide);
            ws_c_du.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_du.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_Gjp.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_Gjp.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_rhoU.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_rhoU.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_Um.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_c_rho.resize(currentNodesPerSide);
            ws_o_rho.resize(opposingNodesPerSide);
            ws_c_alpha.resize(currentNodesPerSide);
            ws_o_alpha.resize(opposingNodesPerSide);
            ws_c_psi.resize(currentNodesPerSide);
            ws_o_psi.resize(opposingNodesPerSide);
            ws_c_face_coordinates.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_c_general_shape_function.resize(currentNodesPerSide);
            ws_o_general_shape_function.resize(opposingNodesPerSide);

            // algorithm related; element; dndx will be at a single
            // gauss point
            ws_c_elem_p.resize(currentNodesPerElement);
            ws_o_elem_p.resize(opposingNodesPerElement);
            ws_c_elem_coordinates.resize(currentNodesPerElement * SPATIAL_DIM);
            ws_o_elem_coordinates.resize(opposingNodesPerElement * SPATIAL_DIM);
            ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
            ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
            ws_c_det_j.resize(1);
            ws_o_det_j.resize(1);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];

            // face
            scalar* p_c_p = &ws_c_p[0];
            scalar* p_o_p = &ws_o_p[0];
            scalar* p_c_du = &ws_c_du[0];
            scalar* p_o_du = &ws_o_du[0];
            scalar* p_c_Gjp = &ws_c_Gjp[0];
            scalar* p_o_Gjp = &ws_o_Gjp[0];
            scalar* p_c_rhoU = &ws_c_rhoU[0];
            scalar* p_o_rhoU = &ws_o_rhoU[0];
            scalar* p_c_Um = &ws_c_Um[0];
            scalar* p_c_rho = &ws_c_rho[0];
            scalar* p_o_rho = &ws_o_rho[0];
            scalar* p_c_alpha = &ws_c_alpha[0];
            scalar* p_o_alpha = &ws_o_alpha[0];
            scalar* p_c_psi = &ws_c_psi[0];
            scalar* p_o_psi = &ws_o_psi[0];
            scalar* p_c_face_coordinates = &ws_c_face_coordinates[0];

            // element
            scalar* p_c_elem_p = &ws_c_elem_p[0];
            scalar* p_o_elem_p = &ws_o_elem_p[0];
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
                p_c_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_c_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_c_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
                p_c_psi[ni] =
                    compressible ? *stk::mesh::field_data(*psiSTKFieldPtr, node)
                                 : 0.0;

                // gather; vector
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * current_num_face_nodes + ni;
                    p_c_rhoU[offSet] = U[i];
                    p_c_Gjp[offSet] = Gjp[i];
                    p_c_face_coordinates[offSet] = coords[i];
                    p_c_du[offSet] = du[i];
                }

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label offSet = i * current_num_face_nodes + ni;
                        p_c_Um[offSet] = Um[i];
                    }
                }
                else
                {
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label offSet = i * current_num_face_nodes + ni;
                        p_c_Um[offSet] = 0.0;
                    }
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
                p_o_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_o_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_o_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
                p_o_psi[ni] =
                    compressible ? *stk::mesh::field_data(*psiSTKFieldPtr, node)
                                 : 0.0;

                // gather; vector
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * opposing_num_face_nodes + ni;
                    p_o_rhoU[offSet] = U[i];
                    p_o_Gjp[offSet] = Gjp[i];
                    p_o_du[offSet] = du[i];
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

                // gather; scalar
                p_c_elem_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label niNdim = ni * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
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

                // gather; scalar
                p_o_elem_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label niNdim = ni * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_o_elem_coordinates[niNdim + i] = coords[i];
                }
            }

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

            // transform opposing normal back to opposing side
            interfaceSideInfoPtr->reverseRotateVector<SPATIAL_DIM>(oNx);

            // project from side to element; method deals with the
            // -1:1 isInElement range to the proper underlying CVFEM
            // range
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

            // bip gradients; zero out
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                cDpdxBip[j] = 0.0;
                oDpdxBip[j] = 0.0;
            }

            // current pressure gradient
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const scalar p = p_c_elem_p[ic];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];
                    cDpdxBip[j] += dndxj * p;
                }
            }

            // opposing pressure gradient
            for (label ic = 0; ic < opposingNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const scalar p = p_o_elem_p[ic];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];
                    oDpdxBip[j] += dndxj * p;
                }
            }

            // product of density and velocity; current (take
            // over previous nodal value for velocity)
            for (label ni = 0; ni < current_num_face_nodes; ++ni)
            {
                const scalar rho = p_c_rho[ni];
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * current_num_face_nodes + ni;
                    p_c_rhoU[offSet] *= rho;
                }
            }

            // opposite
            for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
            {
                const scalar rho = p_o_rho[ni];
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * opposing_num_face_nodes + ni;
                    p_o_rhoU[offSet] *= rho;
                }
            }

            // interpolate to boundary ips
            scalar cPBip = 0.0;
            meFCCurrent->interpolatePoint(
                sizeOfScalarField, &currentIsoParCoords[0], &ws_c_p[0], &cPBip);

            scalar oPBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_p[0],
                                           &oPBip);

            // interpolate density
            scalar cRhoBip = 0.0;
            meFCCurrent->interpolatePoint(sizeOfScalarField,
                                          &currentIsoParCoords[0],
                                          &ws_c_rho[0],
                                          &cRhoBip);

            scalar oRhoBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_rho[0],
                                           &oRhoBip);

            // interpolate volume fraction
            scalar cAlphaBip = 0.0;
            meFCCurrent->interpolatePoint(sizeOfScalarField,
                                          &currentIsoParCoords[0],
                                          &ws_c_alpha[0],
                                          &cAlphaBip);

            scalar oAlphaBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_alpha[0],
                                           &oAlphaBip);

            // interpolate compressibility
            scalar cPsiBip = 0.0;
            meFCCurrent->interpolatePoint(sizeOfScalarField,
                                          &currentIsoParCoords[0],
                                          &ws_c_psi[0],
                                          &cPsiBip);

            scalar oPsiBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_psi[0],
                                           &oPsiBip);

            // interpolate coords
            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_face_coordinates[0],
                                          &currentCoordsBip[0]);

            // interpolate velocity
            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_rhoU[0],
                                          &cRhoUBip[0]);

            meFCOpposing->interpolatePoint(sizeOfVectorField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_rhoU[0],
                                           &oRhoUBip[0]);

            // interpolate mesh velocity
            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_Um[0],
                                          &cUmBip[0]);

            // interpolate pressure diffusivity
            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_du[0],
                                          &cDuBip[0]);

            meFCOpposing->interpolatePoint(sizeOfVectorField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_du[0],
                                           &oDuBip[0]);

            // projected nodal gradient: use arithmetic interpolations: zero-out
            // first
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                cGjpBip[i] = 0;
                oGjpBip[i] = 0;
            }

            for (label ic = 0; ic < currentNodesPerSide; ++ic)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    const label offSet = i * currentNodesPerSide + ic;
                    cGjpBip[i] += f_c * p_c_Gjp[offSet];
                }
            }

            for (label ic = 0; ic < opposingNodesPerSide; ++ic)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    const label offSet = i * opposingNodesPerSide + ic;
                    oGjpBip[i] += f_o * p_o_Gjp[offSet];
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

            // save mDot
            const scalar tmDot = (stk::mesh::field_data(
                mDotSideSTKFieldRef, currentFace))[currentGaussPointId];

            const scalar abs_tmDot = std::abs(tmDot);

            scalar currentDiffBipmag = 0;
            scalar opposingDiffBipmag = 0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                currentDiffBipmag +=
                    cRhoBip * cDuBip[i] / static_cast<scalar>(SPATIAL_DIM);
                opposingDiffBipmag +=
                    oRhoBip * oDuBip[i] / static_cast<scalar>(SPATIAL_DIM);
            }

            const scalar penaltyIp =
                penaltyFactor * 0.5 *
                (currentDiffBipmag * currentInverseLength +
                 opposingDiffBipmag * opposingInverseLength);

            scalar ncFlux = 0.0;
            scalar ncPstabFlux = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar cRhoU = cRhoUBip[j];
                const scalar oRhoU = oRhoUBip[j];
                ncFlux += 0.5 * (cRhoU * p_cNx[j] - oRhoU * p_oNx[j]);

                const scalar cPstab =
                    cRhoBip * cDuBip[j] * (cDpdxBip[j] - cGjpBip[j]);
                const scalar oPstab =
                    oRhoBip * oDuBip[j] * (oDpdxBip[j] - oGjpBip[j]);
                ncPstabFlux += 0.5 * (cPstab * p_cNx[j] - oPstab * p_oNx[j]);
            }

            scalar mDot =
                (ncFlux - ncPstabFlux + penaltyIp * (cPBip - oPBip)) * c_amag;

            // transform mDot to relative frame

            // 1.) frame motion
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    mDot -= cRhoBip * p_mat[i * SPATIAL_DIM + j] *
                            (currentCoordsBip[j] - p_ori[j]) *
                            c_areaVec[currentGaussPointId * SPATIAL_DIM + i];
                }
            }

            // 2.) mesh motion
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                mDot -= cRhoBip * cUmBip[i] * p_cNx[i] * c_amag;
            }

            // form residual
            const label nn = ipNodeMap[currentGaussPointId];

#ifndef NDEBUG
            // SCL check: accumulate grid flux to nearest node
            if (sclCheckSTKFieldPtr)
            {
                scalar gridFlux = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    gridFlux += cRhoBip * cUmBip[i] * p_cNx[i] * c_amag *
                                cAlphaBip / densityScale;
                }
                stk::mesh::Entity node = current_elem_node_rels[nn];
                scalar* scl = stk::mesh::field_data(*sclCheckSTKFieldPtr, node);
                *scl -= gridFlux;
            }
#endif /* NDEBUG */

            p_rhs[nn] -= mDot * cAlphaBip / densityScale;

            // set-up row for matrix
            const label rowR = nn * totalNodes;

            // sensitivities; current face (penalty and advection); use general
            // shape function for this single ip
            scalar lhsFacC = penaltyIp * c_amag + (abs_tmDot + tmDot) / 2.0 *
                                                      cPsiBip / cRhoBip * comp;
            meFCCurrent->general_shape_fcn(
                1, &currentIsoParCoords[0], &ws_c_general_shape_function[0]);
            for (label ic = 0; ic < currentNodesPerSide; ++ic)
            {
                const label icnn = c_face_node_ordinals[ic];
                const scalar r = p_c_general_shape_function[ic];
                p_lhs[rowR + icnn] += r * lhsFacC * cAlphaBip / densityScale;
            }

            // sensitivities; current element (diffusion)
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                scalar lhscd = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_cNx[j];
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];
                    lhscd -= cRhoBip * cDuBip[j] * dndxj * nxj;
                }
                p_lhs[rowR + ic] +=
                    0.5 * lhscd * c_amag * cAlphaBip / densityScale;
            }

            // sensitivities; opposing face (penalty); use general
            // shape function for this single ip
            scalar lhsFacO = penaltyIp * c_amag + (abs_tmDot - tmDot) / 2.0 *
                                                      oPsiBip / oRhoBip * comp;
            meFCOpposing->general_shape_fcn(
                1, &opposingIsoParCoords[0], &ws_o_general_shape_function[0]);
            for (label ic = 0; ic < opposingNodesPerSide; ++ic)
            {
                const label icnn = o_face_node_ordinals[ic];
                const scalar r = p_o_general_shape_function[ic];
                p_lhs[rowR + icnn + currentNodesPerElement] -=
                    r * lhsFacO * oAlphaBip / densityScale;
            }

            // sensitivities; opposing element (diffusion)
            for (label ic = 0; ic < opposingNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                scalar lhscd = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_oNx[j];
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];
                    lhscd -= oRhoBip * oDuBip[j] * dndxj * nxj;
                }
                p_lhs[rowR + ic + currentNodesPerElement] -=
                    0.5 * lhscd * c_amag * oAlphaBip / densityScale;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::assembleElemTermsInterfaceSideNoSlipWall_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // required for frame motion (MFR)
    const auto coriolisMatrix =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* p_mat = coriolisMatrix.data();

    const auto origin =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().origin_
            : utils::vector::Zero();
    const scalar* p_ori = origin.data();

    // mesh motion
    const bool meshMoving = domain->zonePtr()->meshMoving();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> umBip(SPATIAL_DIM);
    std::vector<scalar> coordBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_umBip = &umBip[0];
    scalar* p_coordBip = &coordBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_alpha;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& sideUSTKFieldRef = model_->URef().sideFieldRef().stkFieldRef();

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    scalar densityScale = model_->rhoRef(phaseIndex_).scale();

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

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
        const label nodesPerSide = sideBucket.topology().num_nodes();
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
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
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

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = Um[j];
                    }
                }
                else
                {
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = 0.0;
                    }
                }
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = face_elem_ords[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //======================================
            // gather nodal data off of element
            //======================================
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

                // gather vectors
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                const label nearestNode = ipNodeMap[ip];

                // zero out vector quantities; form aMag
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_umBip[j] = 0.0;
                    p_coordBip[j] = 0.0;
                }

                // interpolate to bip
                scalar rhoBip = 0;
                scalar alphaBip = 0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        p_coordinate_face_shape_function[offSetSF_face + ic];

                    // use velocity shape functions
                    rhoBip += r_vel * p_rho[ic];
                    alphaBip += r_vel * p_alpha[ic];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // use coordinates shape functions
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                        // use velocity shape functions
                        p_umBip[j] += r_vel * p_Um[icNdim + j];
                    }
                }

                // form mDot; rho*uj*Aj
                scalar mDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                    mDot += (rhoBip * UbcVec[ip * SPATIAL_DIM + j]) * axj;
                }

                // transform mDot to relative frame

                // 1.) frame motion
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        mDot -= rhoBip * p_mat[i * SPATIAL_DIM + j] *
                                (p_coordBip[j] - p_ori[j]) *
                                areaVec[ip * SPATIAL_DIM + i];
                    }
                }

                // 2.) mesh motion
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    mDot -= rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                }

#ifndef NDEBUG
                // SCL check: accumulate grid flux to nearest node
                if (sclCheckSTKFieldPtr)
                {
                    scalar gridFlux = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        gridFlux += rhoBip * p_umBip[i] *
                                    areaVec[ip * SPATIAL_DIM + i] * alphaBip /
                                    densityScale;
                    }
                    stk::mesh::Entity node = elemNodeRels[nearestNode];
                    scalar* scl =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, node);
                    *scl -= gridFlux;
                }
#endif /* NDEBUG */

                // residual
                p_rhs[nearestNode] -= mDot * alphaBip / densityScale;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif /* HAS_INTERFACE */
