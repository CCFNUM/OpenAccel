// File       : solidDisplacementAssemblerElemInterfaceConditions.cpp
// Created    : Thu Dec 04 2025 12:30:00 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Nonconformal interface treatment for solid displacement
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "interface.h"
#include "ipInfo.h"
#include "solidDisplacementAssembler.h"

namespace accel
{

namespace
{

bool assembleFEMInterfacePenalty(
    const label currentNodes,
    const label opposingNodes,
    const scalar* currentShape,
    const scalar* opposingShape,
    const scalar* opposingToCurrent,
    const scalar penaltyMeasure,
    const scalar* currentDisplacement,
    const scalar* opposingDisplacement,
    scalar* currentResidual,
    scalar* currentCurrentTangent,
    scalar* currentOpposingTangent)
{
    if (currentNodes <= 0 || opposingNodes <= 0 || !currentShape ||
        !opposingShape || !opposingToCurrent || !currentDisplacement ||
        !opposingDisplacement || !currentResidual ||
        !currentCurrentTangent || !currentOpposingTangent)
    {
        return false;
    }

    const label currentDofs = currentNodes * SPATIAL_DIM;
    const label opposingDofs = opposingNodes * SPATIAL_DIM;
    std::fill(currentResidual, currentResidual + currentDofs, 0.0);
    std::fill(currentCurrentTangent,
              currentCurrentTangent + currentDofs * currentDofs,
              0.0);
    std::fill(currentOpposingTangent,
              currentOpposingTangent + currentDofs * opposingDofs,
              0.0);

    scalar currentValue[SPATIAL_DIM] = {};
    scalar opposingValue[SPATIAL_DIM] = {};
    scalar jump[SPATIAL_DIM] = {};

    for (label node = 0; node < currentNodes; ++node)
    {
        for (label component = 0; component < SPATIAL_DIM; ++component)
        {
            currentValue[component] +=
                currentShape[node] *
                currentDisplacement[node * SPATIAL_DIM + component];
        }
    }

    for (label node = 0; node < opposingNodes; ++node)
    {
        for (label component = 0; component < SPATIAL_DIM; ++component)
        {
            opposingValue[component] +=
                opposingShape[node] *
                opposingDisplacement[node * SPATIAL_DIM + component];
        }
    }

    for (label row = 0; row < SPATIAL_DIM; ++row)
    {
        jump[row] = currentValue[row];
        for (label column = 0; column < SPATIAL_DIM; ++column)
        {
            jump[row] -=
                opposingToCurrent[row * SPATIAL_DIM + column] *
                opposingValue[column];
        }
    }

    for (label currentTest = 0; currentTest < currentNodes; ++currentTest)
    {
        for (label rowComponent = 0; rowComponent < SPATIAL_DIM;
             ++rowComponent)
        {
            const label row =
                currentTest * SPATIAL_DIM + rowComponent;
            const scalar testMeasure =
                penaltyMeasure * currentShape[currentTest];
            currentResidual[row] = testMeasure * jump[rowComponent];

            for (label currentTrial = 0; currentTrial < currentNodes;
                 ++currentTrial)
            {
                const label column =
                    currentTrial * SPATIAL_DIM + rowComponent;
                currentCurrentTangent[row * currentDofs + column] =
                    testMeasure * currentShape[currentTrial];
            }

            for (label opposingTrial = 0; opposingTrial < opposingNodes;
                 ++opposingTrial)
            {
                for (label columnComponent = 0;
                     columnComponent < SPATIAL_DIM;
                     ++columnComponent)
                {
                    const label column =
                        opposingTrial * SPATIAL_DIM + columnComponent;
                    currentOpposingTangent[row * opposingDofs + column] =
                        -testMeasure * opposingShape[opposingTrial] *
                        opposingToCurrent[rowComponent * SPATIAL_DIM +
                                          columnComponent];
                }
            }
        }
    }

    return true;
}

} // namespace

void solidDisplacementAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    if (field_broker_->controlsRef().isCvfemSolidMechanics())
    {
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

        const solidMechanicsOption solidMechOption =
            domain->solidMechanics_.option_;

        // rotation matrix (in case of rotational periodicity)
        const auto& rotMat = interfaceSideInfoPtr->rotationMatrix_;

        // GGI path is not yet validated for the solid-displacement problem;
        // preserve the original errorMsg behavior at the top of the function
        // and run the unified loop for DG below.

        // Unified loop over per-IP info.  Storage is owned by the base
        // interfaceSideInfo; concrete side classes store derived records upcast
        // to ipInfo*, and ip->areaFraction_ is the default 1.0 so the math is
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
                const label lhsSize =
                    totalNodes * SPATIAL_DIM * totalNodes * SPATIAL_DIM;
                const label rhsSize = totalNodes * SPATIAL_DIM;
                lhs.resize(lhsSize);
                rhs.resize(rhsSize);
                scratchIds.resize(rhsSize);
                scratchVals.resize(rhsSize);
                connectedNodes.resize(totalNodes);

                // algorithm related; element; dndx will be at a
                // single gauss point
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

                // compute average stiffness for penalty
                if (solidMechOption == solidMechanicsOption::linearElastic)
                {
                    currentStiffness =
                        currentLambdaBip + 2.0 * currentMuBip / 3.0;
                    opposingStiffness =
                        opposingLambdaBip + 2.0 * opposingMuBip / 3.0;
                }
                else if (solidMechOption ==
                         solidMechanicsOption::neoHookean)
                {
                    // Gradient tensor F = I + grad(u) at interface
                    // BIP, current side. Layout from
                    // general_face_grad_op (1 IP): p_c_dndx[ic *
                    // SPATIAL_DIM + j] = dN_ic/dx_j
                    scalar c_gradTens[SPATIAL_DIM * SPATIAL_DIM] = {};
                    for (label ic = 0; ic < current_num_elem_nodes; ++ic)
                    {
                        const label offSet = ic * SPATIAL_DIM;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dxi =
                                p_c_elem_disp[ic * SPATIAL_DIM + i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                                c_gradTens[i * SPATIAL_DIM + j] +=
                                    p_c_dndx[offSet + j] * dxi;
                        }
                    }
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                        c_gradTens[i * SPATIAL_DIM + i] += 1.0;

                    // Gradient tensor F = I + grad(u), opposing
                    // side
                    scalar o_gradTens[SPATIAL_DIM * SPATIAL_DIM] = {};
                    for (label ic = 0; ic < opposing_num_elem_nodes; ++ic)
                    {
                        const label offSet = ic * SPATIAL_DIM;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dxi =
                                p_o_elem_disp[ic * SPATIAL_DIM + i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                                o_gradTens[i * SPATIAL_DIM + j] +=
                                    p_o_dndx[offSet + j] * dxi;
                        }
                    }
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                        o_gradTens[i * SPATIAL_DIM + i] += 1.0;

                    // I1 = tr(F^T F), beta = 1 - 1/(1+I1), J =
                    // det(F)
                    scalar c_Iconst = 0.0;
                    scalar o_Iconst = 0.0;
                    if (SPATIAL_DIM == 3)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            c_Iconst += c_gradTens[i] * c_gradTens[i] +
                                        c_gradTens[i + 3] * c_gradTens[i + 3] +
                                        c_gradTens[i + 6] * c_gradTens[i + 6];
                            o_Iconst += o_gradTens[i] * o_gradTens[i] +
                                        o_gradTens[i + 3] * o_gradTens[i + 3] +
                                        o_gradTens[i + 6] * o_gradTens[i + 6];
                        }
                    }
                    else
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            c_Iconst += c_gradTens[i] * c_gradTens[i] +
                                        c_gradTens[i + 2] * c_gradTens[i + 2];
                            o_Iconst += o_gradTens[i] * o_gradTens[i] +
                                        o_gradTens[i + 2] * o_gradTens[i + 2];
                        }
                    }

                    const scalar c_beta = 1.0 - 1.0 / (1.0 + c_Iconst);
                    const scalar o_beta = 1.0 - 1.0 / (1.0 + o_Iconst);

                    scalar c_J;
                    scalar o_J;
                    if (SPATIAL_DIM == 3)
                    {
                        c_J = c_gradTens[0] * (c_gradTens[4] * c_gradTens[8] -
                                               c_gradTens[7] * c_gradTens[5]) -
                              c_gradTens[1] * (c_gradTens[3] * c_gradTens[8] -
                                               c_gradTens[6] * c_gradTens[5]) +
                              c_gradTens[2] * (c_gradTens[3] * c_gradTens[7] -
                                               c_gradTens[6] * c_gradTens[4]);
                        o_J = o_gradTens[0] * (o_gradTens[4] * o_gradTens[8] -
                                               o_gradTens[7] * o_gradTens[5]) -
                              o_gradTens[1] * (o_gradTens[3] * o_gradTens[8] -
                                               o_gradTens[6] * o_gradTens[5]) +
                              o_gradTens[2] * (o_gradTens[3] * o_gradTens[7] -
                                               o_gradTens[6] * o_gradTens[4]);
                    }
                    else
                    {
                        c_J = c_gradTens[0] * c_gradTens[3] -
                              c_gradTens[2] * c_gradTens[1];
                        o_J = o_gradTens[0] * o_gradTens[3] -
                              o_gradTens[2] * o_gradTens[1];
                    }

                    // Effective tangent stiffness: lambda +
                    // 2*(mu*beta/J)/3
                    currentStiffness =
                        currentLambdaBip +
                        2.0 * (currentMuBip * c_beta / c_J) / 3.0;
                    opposingStiffness =
                        opposingLambdaBip +
                        2.0 * (opposingMuBip * o_beta / o_J) / 3.0;
                }

                // penalty computation (similar to diffusion
                // penalty)
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

                // compute general shape function at current and
                // opposing integration points
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
                    // penalty term: enforces displacement
                    // continuity
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
    }
    else
    {
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const auto& displacementField = phi_->stkFieldRef();
    const auto& areaVectorField = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
    {
        const auto& tractionField = phi_->sideFluxFieldRef().stkFieldRef();
        const stk::mesh::Selector selectedSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);
        const auto& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selectedSides);

        std::vector<scalar> shapeFunctions;
        for (const stk::mesh::Bucket* bucket : sideBuckets)
        {
            MasterElement* faceMasterElement =
                MasterElementRepo::get_surface_master_element(
                    bucket->topology());
            const label nodesPerSide = faceMasterElement->nodesPerElement_;
            const label integrationPoints = faceMasterElement->numIntPoints_;
            const label dofsPerSide = nodesPerSide * SPATIAL_DIM;

            lhs.assign(dofsPerSide * dofsPerSide, 0.0);
            rhs.resize(dofsPerSide);
            scratchIds.resize(dofsPerSide);
            scratchVals.resize(dofsPerSide);
            connectedNodes.resize(nodesPerSide);
            shapeFunctions.resize(integrationPoints * nodesPerSide);

            if (phi_->isShifted())
                faceMasterElement->shifted_shape_fcn(shapeFunctions.data());
            else
                faceMasterElement->shape_fcn(shapeFunctions.data());

            for (stk::mesh::Entity side : *bucket)
            {
                std::fill(rhs.begin(), rhs.end(), 0.0);

                const stk::mesh::Entity* sideNodes =
                    bulkData.begin_nodes(side);
                STK_ThrowAssert(bulkData.num_nodes(side) == nodesPerSide);
                for (label node = 0; node < nodesPerSide; ++node)
                    connectedNodes[node] = sideNodes[node];

                const scalar* traction =
                    stk::mesh::field_data(tractionField, side);
                const scalar* areaVector =
                    stk::mesh::field_data(areaVectorField, side);

                for (label ip = 0; ip < integrationPoints; ++ip)
                {
                    scalar areaMagnitudeSquared = 0.0;
                    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                    {
                        const scalar area =
                            areaVector[ip * SPATIAL_DIM + dim];
                        areaMagnitudeSquared += area * area;
                    }
                    const scalar areaMagnitude =
                        std::sqrt(areaMagnitudeSquared);

                    for (label node = 0; node < nodesPerSide; ++node)
                    {
                        const scalar weight =
                            shapeFunctions[ip * nodesPerSide + node] *
                            areaMagnitude;
                        for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                        {
                            rhs[node * SPATIAL_DIM + dim] +=
                                weight *
                                traction[ip * SPATIAL_DIM + dim];
                        }
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
        return;
    }

    STK_ThrowRequireMsg(
        domain->solidMechanics_.option_ == solidMechanicsOption::linearElastic,
        "The FEM solid interface currently supports linear_elastic only");

    const auto& EField = model_->ERef().stkFieldRef();
    const auto& nuField = model_->nuRef().stkFieldRef();
    const auto& coordinatesField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const bool planeStress = domain->solidMechanics_.planeStress_;
    const scalar penaltyFactor =
        interfaceSideInfoPtr->interfPtr()->penaltyFactor();
    const scalar* opposingToCurrent =
        interfaceSideInfoPtr->rotationMatrix_.data();

    std::vector<scalar> currentFaceDisplacement;
    std::vector<scalar> opposingFaceDisplacement;
    std::vector<scalar> currentElementCoordinates;
    std::vector<scalar> opposingElementCoordinates;
    std::vector<scalar> currentDndx;
    std::vector<scalar> opposingDndx;
    std::vector<scalar> currentDetJ(1);
    std::vector<scalar> opposingDetJ(1);
    std::vector<scalar> currentShape;
    std::vector<scalar> opposingShape;
    std::vector<scalar> currentResidual;
    std::vector<scalar> currentCurrentTangent;
    std::vector<scalar> currentOpposingTangent;
    std::vector<scalar> currentElementIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingElementIsoParCoords(SPATIAL_DIM);

    for (const auto& faceIpInfoVec : interfaceSideInfoPtr->ipInfoVec())
    {
        for (const ipInfo* ip : faceIpInfoVec)
        {
            if (ip->isExposed_)
                continue;

            MasterElement* currentFaceME = ip->meFCCurrent_;
            MasterElement* opposingFaceME = ip->meFCOpposing_;
            MasterElement* currentElementME = ip->meSCSCurrent_;
            MasterElement* opposingElementME = ip->meSCSOpposing_;

            const label currentFaceNodes = currentFaceME->nodesPerElement_;
            const label opposingFaceNodes = opposingFaceME->nodesPerElement_;
            const label currentElementNodes =
                currentElementME->nodesPerElement_;
            const label opposingElementNodes =
                opposingElementME->nodesPerElement_;
            const label totalNodes =
                currentElementNodes + opposingElementNodes;
            const label totalDofs = totalNodes * SPATIAL_DIM;

            lhs.assign(totalDofs * totalDofs, 0.0);
            rhs.assign(totalDofs, 0.0);
            scratchIds.resize(totalDofs);
            scratchVals.resize(totalDofs);
            connectedNodes.resize(totalNodes);

            currentFaceDisplacement.resize(currentFaceNodes * SPATIAL_DIM);
            opposingFaceDisplacement.resize(opposingFaceNodes * SPATIAL_DIM);
            currentElementCoordinates.resize(currentElementNodes *
                                             SPATIAL_DIM);
            opposingElementCoordinates.resize(opposingElementNodes *
                                              SPATIAL_DIM);
            currentDndx.resize(currentElementNodes * SPATIAL_DIM);
            opposingDndx.resize(opposingElementNodes * SPATIAL_DIM);
            currentShape.resize(currentFaceNodes);
            opposingShape.resize(opposingFaceNodes);

            const label currentFaceDofs = currentFaceNodes * SPATIAL_DIM;
            const label opposingFaceDofs = opposingFaceNodes * SPATIAL_DIM;
            currentResidual.resize(currentFaceDofs);
            currentCurrentTangent.resize(currentFaceDofs * currentFaceDofs);
            currentOpposingTangent.resize(currentFaceDofs *
                                          opposingFaceDofs);

            const stk::mesh::Entity* currentFaceNodeRels =
                bulkData.begin_nodes(ip->currentFace_);
            const stk::mesh::Entity* opposingFaceNodeRels =
                bulkData.begin_nodes(ip->opposingFace_);
            STK_ThrowAssert(bulkData.num_nodes(ip->currentFace_) ==
                            currentFaceNodes);
            STK_ThrowAssert(bulkData.num_nodes(ip->opposingFace_) ==
                            opposingFaceNodes);

            scalar currentMu = 0.0;
            scalar currentLambda = 0.0;
            for (label node = 0; node < currentFaceNodes; ++node)
            {
                const scalar* displacement = stk::mesh::field_data(
                    displacementField, currentFaceNodeRels[node]);
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    currentFaceDisplacement[node * SPATIAL_DIM + dim] =
                        displacement[dim];
                }
            }

            scalar opposingMu = 0.0;
            scalar opposingLambda = 0.0;
            for (label node = 0; node < opposingFaceNodes; ++node)
            {
                const scalar* displacement = stk::mesh::field_data(
                    displacementField, opposingFaceNodeRels[node]);
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    opposingFaceDisplacement[node * SPATIAL_DIM + dim] =
                        displacement[dim];
                }
            }

            currentFaceME->general_shape_fcn(
                1, ip->currentIsoParCoords_.data(), currentShape.data());
            opposingFaceME->general_shape_fcn(
                1, ip->opposingIsoParCoords_.data(), opposingShape.data());

            for (label node = 0; node < currentFaceNodes; ++node)
            {
                const scalar E =
                    *stk::mesh::field_data(EField, currentFaceNodeRels[node]);
                const scalar nu =
                    *stk::mesh::field_data(nuField, currentFaceNodeRels[node]);
                const scalar shape = currentShape[node];
                currentMu += shape * E / (2.0 * (1.0 + nu));
                currentLambda +=
                    shape *
                    (planeStress
                         ? nu * E / ((1.0 + nu) * (1.0 - nu))
                         : nu * E /
                               ((1.0 + nu) * (1.0 - 2.0 * nu)));
            }

            for (label node = 0; node < opposingFaceNodes; ++node)
            {
                const scalar E =
                    *stk::mesh::field_data(EField, opposingFaceNodeRels[node]);
                const scalar nu =
                    *stk::mesh::field_data(nuField,
                                           opposingFaceNodeRels[node]);
                const scalar shape = opposingShape[node];
                opposingMu += shape * E / (2.0 * (1.0 + nu));
                opposingLambda +=
                    shape *
                    (planeStress
                         ? nu * E / ((1.0 + nu) * (1.0 - nu))
                         : nu * E /
                               ((1.0 + nu) * (1.0 - 2.0 * nu)));
            }

            const stk::mesh::Entity* currentElementNodeRels =
                bulkData.begin_nodes(ip->currentElement_);
            const stk::mesh::Entity* opposingElementNodeRels =
                bulkData.begin_nodes(ip->opposingElement_);
            STK_ThrowAssert(bulkData.num_nodes(ip->currentElement_) ==
                            currentElementNodes);
            STK_ThrowAssert(bulkData.num_nodes(ip->opposingElement_) ==
                            opposingElementNodes);

            for (label node = 0; node < currentElementNodes; ++node)
            {
                connectedNodes[node] = currentElementNodeRels[node];
                const scalar* coordinates = stk::mesh::field_data(
                    coordinatesField, currentElementNodeRels[node]);
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    currentElementCoordinates[node * SPATIAL_DIM + dim] =
                        coordinates[dim];
                }
            }

            for (label node = 0; node < opposingElementNodes; ++node)
            {
                connectedNodes[currentElementNodes + node] =
                    opposingElementNodeRels[node];
                const scalar* coordinates = stk::mesh::field_data(
                    coordinatesField, opposingElementNodeRels[node]);
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    opposingElementCoordinates[node * SPATIAL_DIM + dim] =
                        coordinates[dim];
                }
            }
            interfaceSideInfoPtr->transformCoordinateList(
                opposingElementCoordinates, opposingElementNodes);

            currentElementME->sidePcoords_to_elemPcoords(
                ip->currentFaceOrdinal_,
                1,
                ip->currentIsoParCoords_.data(),
                currentElementIsoParCoords.data());
            opposingElementME->sidePcoords_to_elemPcoords(
                ip->opposingFaceOrdinal_,
                1,
                ip->opposingIsoParCoords_.data(),
                opposingElementIsoParCoords.data());

            scalar gradientError = 0.0;
            currentElementME->general_face_grad_op(
                ip->currentFaceOrdinal_,
                currentElementIsoParCoords.data(),
                currentElementCoordinates.data(),
                currentDndx.data(),
                currentDetJ.data(),
                &gradientError);
            opposingElementME->general_face_grad_op(
                ip->opposingFaceOrdinal_,
                opposingElementIsoParCoords.data(),
                opposingElementCoordinates.data(),
                opposingDndx.data(),
                opposingDetJ.data(),
                &gradientError);
            STK_ThrowRequireMsg(
                gradientError == 0.0,
                "Invalid element geometry at FEM solid interface");

            const scalar* areaVector =
                stk::mesh::field_data(areaVectorField, ip->currentFace_);
            scalar areaMagnitudeSquared = 0.0;
            std::vector<scalar> currentNormal(SPATIAL_DIM);
            for (label dim = 0; dim < SPATIAL_DIM; ++dim)
            {
                const scalar area =
                    areaVector[ip->currentGaussPointId_ * SPATIAL_DIM + dim];
                areaMagnitudeSquared += area * area;
                currentNormal[dim] = area;
            }
            const scalar areaMagnitude =
                std::sqrt(areaMagnitudeSquared);
            STK_ThrowRequireMsg(
                areaMagnitude > std::numeric_limits<scalar>::epsilon(),
                "Zero-area integration point at FEM solid interface");
            for (scalar& component : currentNormal)
                component /= areaMagnitude;

            const label* currentFaceOrdinals =
                currentElementME->side_node_ordinals(
                    ip->currentFaceOrdinal_);
            const label* opposingFaceOrdinals =
                opposingElementME->side_node_ordinals(
                    ip->opposingFaceOrdinal_);

            scalar currentInverseLength = 0.0;
            for (label node = 0; node < currentFaceNodes; ++node)
            {
                const label elementNode = currentFaceOrdinals[node];
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    currentInverseLength +=
                        currentDndx[elementNode * SPATIAL_DIM + dim] *
                        currentNormal[dim];
                }
            }

            std::vector<scalar> opposingNormal(currentNormal);
            for (scalar& component : opposingNormal)
                component = -component;

            scalar opposingInverseLength = 0.0;
            for (label node = 0; node < opposingFaceNodes; ++node)
            {
                const label elementNode = opposingFaceOrdinals[node];
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    opposingInverseLength +=
                        opposingDndx[elementNode * SPATIAL_DIM + dim] *
                        opposingNormal[dim];
                }
            }

            const scalar currentStiffness =
                currentLambda + 2.0 * currentMu / 3.0;
            const scalar opposingStiffness =
                opposingLambda + 2.0 * opposingMu / 3.0;
            const scalar penalty =
                penaltyFactor *
                (currentStiffness * currentInverseLength +
                 opposingStiffness * opposingInverseLength) /
                2.0;
            const scalar penaltyMeasure =
                penalty * areaMagnitude * ip->areaFraction_;

            const bool interfaceStatus =
                assembleFEMInterfacePenalty(
                    currentFaceNodes,
                    opposingFaceNodes,
                    currentShape.data(),
                    opposingShape.data(),
                    opposingToCurrent,
                    penaltyMeasure,
                    currentFaceDisplacement.data(),
                    opposingFaceDisplacement.data(),
                    currentResidual.data(),
                    currentCurrentTangent.data(),
                    currentOpposingTangent.data());
            STK_ThrowRequireMsg(
                interfaceStatus,
                "OpenAccel FEM solid interface kernel failed");

            for (label currentTest = 0; currentTest < currentFaceNodes;
                 ++currentTest)
            {
                const label currentTestElementNode =
                    currentFaceOrdinals[currentTest];
                for (label rowComponent = 0;
                     rowComponent < SPATIAL_DIM;
                     ++rowComponent)
                {
                    const label faceRow =
                        currentTest * SPATIAL_DIM + rowComponent;
                    const label elementRow =
                        currentTestElementNode * SPATIAL_DIM + rowComponent;
                    rhs[elementRow] -= currentResidual[faceRow];

                    for (label currentTrial = 0;
                         currentTrial < currentFaceNodes;
                         ++currentTrial)
                    {
                        const label currentTrialElementNode =
                            currentFaceOrdinals[currentTrial];
                        for (label columnComponent = 0;
                             columnComponent < SPATIAL_DIM;
                             ++columnComponent)
                        {
                            const label faceColumn =
                                currentTrial * SPATIAL_DIM + columnComponent;
                            const label elementColumn =
                                currentTrialElementNode * SPATIAL_DIM +
                                columnComponent;
                            lhs[elementRow * totalDofs + elementColumn] +=
                                currentCurrentTangent
                                    [faceRow * currentFaceDofs + faceColumn];
                        }
                    }

                    for (label opposingTrial = 0;
                         opposingTrial < opposingFaceNodes;
                         ++opposingTrial)
                    {
                        const label opposingTrialElementNode =
                            opposingFaceOrdinals[opposingTrial];
                        for (label columnComponent = 0;
                             columnComponent < SPATIAL_DIM;
                             ++columnComponent)
                        {
                            const label faceColumn =
                                opposingTrial * SPATIAL_DIM + columnComponent;
                            const label elementColumn =
                                (currentElementNodes +
                                 opposingTrialElementNode) *
                                    SPATIAL_DIM +
                                columnComponent;
                            lhs[elementRow * totalDofs + elementColumn] +=
                                currentOpposingTangent
                                    [faceRow * opposingFaceDofs + faceColumn];
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
