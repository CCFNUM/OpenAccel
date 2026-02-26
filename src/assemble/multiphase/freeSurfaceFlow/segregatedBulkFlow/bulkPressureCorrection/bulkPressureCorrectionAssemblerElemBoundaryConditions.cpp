// File       : bulkPressureCorrectionAssemblerElemBoundaryConditions.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "bulkPressureCorrectionAssembler.h"
#include "freeSurfaceFlowModel.h"
#include "zoneTransformation.h"

namespace accel
{

void bulkPressureCorrectionAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
{
    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        const boundaryConditionType pBCType =
            model_->pRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();
        const boundaryConditionType UBCType =
            model_->URef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                {
                    assembleElemTermsBoundarySymmetry_(domain, boundary, ctx);
                }
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::slip:
                            {
                                assembleElemTermsBoundarySymmetry_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::noSlip:
                            {
                                assembleElemTermsBoundaryWallNoSlip_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("invalid velocity boundary "
                                     "condition at wall");
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (pBCType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::specifiedValue:
                                    case boundaryConditionType::normalSpeed:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedVelocity_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    case boundaryConditionType::massFlowRate:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedMassFlowRate_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        case boundaryConditionType::staticPressure:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::specifiedValue:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    case boundaryConditionType::
                                        specifiedDirection:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        case boundaryConditionType::totalPressure:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::
                                        specifiedDirection:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        default:
                            errorMsg("invalid pressure boundary "
                                     "condition at inlet");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (pBCType)
                    {
                        case boundaryConditionType::staticPressure:
                        case boundaryConditionType::averageStaticPressure:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::zeroGradient:
                                        {
                                            assembleElemTermsBoundaryOutletSpecifiedPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at outlet");
                                }
                            }
                            break;

                        case boundaryConditionType::massFlowRate:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::massFlowRate:
                                        {
                                            assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at outlet");
                                }
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            {
                                switch (UBCType)
                                {
                                    case boundaryConditionType::zeroGradient:
                                        {
                                            assembleElemTermsBoundaryOutletOutflow_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid velocity boundary "
                                                 "condition at outlet");
                                }
                            }
                            break;

                        default:
                            errorMsg("invalid pressure boundary "
                                     "condition at outlet");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    switch (pBCType)
                    {
                        case boundaryConditionType::staticPressure:
                        case boundaryConditionType::totalPressure:
                            {
                                assembleElemTermsBoundaryOpening_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("invalid pressure boundary "
                                     "condition at outlet");
                    }
                    break;
                }
        }
    }
}

void bulkPressureCorrectionAssembler::assembleElemTermsBoundarySymmetry_(
    const domain* domain,
    const boundary* boundary,
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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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
                const label offSetSF_face = ip * nodesPerSide;
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

                // form mDot; rho*uj*Aj = 0
                scalar mDot = 0.0;

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

void bulkPressureCorrectionAssembler::assembleElemTermsBoundaryWallNoSlip_(
    const domain* domain,
    const boundary* boundary,
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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryInletSpecifiedVelocity_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const bool compressible = domain->isMaterialCompressible();
    const scalar comp = compressible ? 1.0 : 0.0;

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
    std::vector<scalar> ws_psi;
    std::vector<scalar> ws_alpha;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& sideUSTKFieldRef = model_->URef().sideFieldRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef(phaseIndex_).sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    const auto* psiSTKFieldPtr =
        compressible ? model_->psiRef().stkFieldPtr() : nullptr;

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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
        ws_psi.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_psi = &ws_psi[0];
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
                p_psi[ni] = compressible
                                ? *stk::mesh::field_data(*psiSTKFieldPtr, node)
                                : 0;

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
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

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
                const scalar tmDot =
                    (stk::mesh::field_data(mDotSideSTKFieldRef, side))[ip];

                if (rfflag[ip] == 0)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // zero out vector quantities; form aMag
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_umBip[j] = 0.0;
                        p_coordBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar rhoBip = 0;
                    scalar psiBip = 0;
                    scalar alphaBip = 0;
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];
                        const scalar r_coord =
                            p_coordinate_face_shape_function[offSetSF_face +
                                                             ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        psiBip += r_vel * p_psi[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[ic * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho*uspecj*Aj
                    scalar mDot = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        mDot += (rhoBip * UbcVec[ip * SPATIAL_DIM + j]) * axj;
                    }

                    //================================
                    // Compressibility contribution at inlet
                    // Newton-Raphson: ∂(ρU_bc)/∂p = U_bc * ∂ρ/∂p = U_bc * ψ
                    // LHS coefficient = mDot / ρ_bip * ψ_bip
                    // where ρ_bip and ψ_bip are interpolated to boundary ip
                    //================================
                    label rowR = nearestNode * nodesPerElement;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        p_lhs[rowR + inn] += tmDot * r_vel * psiBip / rhoBip *
                                             comp * alphaBip / densityScale;
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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
                else
                {
                    // slip-wall

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
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_vel * p_coordinates[inn * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho*uj*Aj = 0
                    scalar mDot = 0.0;

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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryInletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_alpha;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;

    // Get fields
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef(phaseIndex_).sideFieldRef().stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();

    scalar densityScale = model_->rhoRef(phaseIndex_).scale();

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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
        ws_alpha.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

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

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];

            // get element; its face ordinal number and populate
            // face_node_ordinals
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = face_elem_ords[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

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
            }

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
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const scalar tmDot =
                    (stk::mesh::field_data(mDotSideSTKFieldRef, side))[ip];

                const label nearestNode = ipNodeMap[ip];

                // interpolate to bip
                scalar alphaBip = 0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF_face + ic];

                    // use velocity shape functions
                    alphaBip += r_vel * p_alpha[ic];
                }

                // residual
                p_rhs[nearestNode] -= tmDot * alphaBip / densityScale;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
        const domain* domain,
        const boundary* boundary,
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
    const auto& sideRhoSTKFieldRef =
        model_->rhoRef(phaseIndex_).sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            const scalar* rhoBc =
                stk::mesh::field_data(sideRhoSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

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

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                if (rfflag[ip] == 1)
                {
                    errorMsg("invalid outflow at velocity-inlet");
                }
                else
                {
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
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];
                        const scalar r_coord =
                            p_coordinate_face_shape_function[offSetSF_face +
                                                             ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[ic * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho_spec*u_spec*Aj
                    scalar mDot = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        mDot +=
                            (rhoBc[ip] * UbcVec[ip * SPATIAL_DIM + j]) * axj;
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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryOutletSpecifiedPressure_(const domain* domain,
                                                      const boundary* boundary,
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
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> umBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> duRhsBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_umBip = &umBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_duRhsBip = &duRhsBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_duRhs;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_bcMultiplier;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& nodalSidePSTKFieldRef =
        model_->pRef().nodeSideFieldRef().stkFieldRef();
    const auto& gradPSTKFieldRef = model_->pRef().gradRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Pressure diffusivity: duTilde for LHS (SIMPLEC), du for RHS
    const bool consistent = model_->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    const auto& duLhsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        consistent ? flowModel::duTilde_ID : flowModel::du_ID);
    const auto& duRhsSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get body force fields for buoyancy pressure stabilization
    const auto& FSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto& FOrigSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = model_->pRef().isGradientShifted();

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
        const scalar f = 1.0 / static_cast<scalar>(nodesPerSide);

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_duRhs.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_duRhs = &ws_duRhs[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_bcMultiplier[ni] = 1.0;

                // gather vectors
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du =
                    stk::mesh::field_data(duLhsSTKFieldRef, node);
                const scalar* duRhs =
                    stk::mesh::field_data(duRhsSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                    p_duRhs[offSet + j] = duRhs[j];
                }

                // overwrite pressure value at boundary with dirichlet value
                p_p[ic] = *stk::mesh::field_data(nodalSidePSTKFieldRef, node);

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

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(FSTKFieldRef, node);
                const scalar* FOrig =
                    stk::mesh::field_data(FOrigSTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }
            }

            // compute dndx for residual
            scalar scs_error = 0.0;
            if (isPGradientShifted)
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

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                if (rfflag[ip] == 0)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // zero out vector quantities; form aMag
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] = 0.0;
                        p_uBip[j] = 0.0;
                        p_umBip[j] = 0.0;
                        p_GpdxBip[j] = 0.0;
                        p_dpdxBip[j] = 0.0;
                        p_duBip[j] = 0.0;
                        p_duRhsBip[j] = 0.0;
                        p_FBip[j] = 0.0;
                        p_FOrigBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar rhoBip = 0;
                    scalar alphaBip = 0;
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];
                        const scalar r_coord =
                            p_coordinate_face_shape_function[offSetSF_face +
                                                             ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_uBip[j] += r_vel * p_U[icNdim + j];
                            p_duBip[j] += r_vel * p_du[icNdim + j];
                            p_duRhsBip[j] += r_vel * p_duRhs[icNdim + j];
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                            // face-centre average for original body force
                            p_FOrigBip[j] += r_vel * p_FOrig[icNdim + j];

                            // arithmetic interpolation
                            p_GpdxBip[j] += f * p_Gpdx[icNdim + j];

                            // interpolate redistributed body force using shape
                            // functions
                            p_FBip[j] += f * p_F[icNdim + j];
                        }
                    }

                    // form dpdxBip
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;
                        const scalar p = p_p[ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dpdxBip[j] += p_dndx[offSetDnDx + j] * p;
                        }
                    }

                    // lhs for pressure system
                    label rowR = nearestNode * nodesPerElement;

                    //================================
                    // laplacian: -rhoip*Dip*Gpn.Sip
                    //================================

                    // element-based gradient: here we add contributions from
                    // all nodes. However, due to nature of boundary
                    // (dirichlet), boundary nodes must not contribute to the
                    // lhs, so we multiply by 0 the lhs when the current node
                    // belongs to the boundary
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;

                        scalar lhsFac = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            lhsFac += -rhoBip * p_duBip[j] *
                                      p_dndx[offSetDnDx + j] *
                                      areaVec[ip * SPATIAL_DIM + j];
                        }
                        p_lhs[rowR + ic] += lhsFac * p_bcMultiplier[ic] *
                                            alphaBip / densityScale;
                    }

                    //================================
                    // All-in-all: form mDot rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj
                    // + rho*du*(FOrigj - Fj)*Aj (buoyancy
                    // stabilization)
                    //================================
                    scalar mDot = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                        // divergence + pressure Rhie-Chow (D for RHS)
                        mDot += (rhoBip * p_uBip[j] -
                                 rhoBip * p_duRhsBip[j] *
                                     (p_dpdxBip[j] - p_GpdxBip[j])) *
                                axj;

                        // // buoyancy stabilization: +rho*D*(F_orig - F)·S
                        // mDot += rhoBip * p_duRhsBip[j] *
                        //         (p_FOrigBip[j] - p_FBip[j]) * axj;
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
                        mDot -= rhoBip * p_umBip[i] *
                                areaVec[ip * SPATIAL_DIM + i] * alphaBip /
                                densityScale;
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i];
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
                else
                {
                    // slip-wall

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
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_vel * p_coordinates[inn * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho*uj*Aj = 0
                    scalar mDot = 0.0;

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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::assembleElemTermsBoundaryOutletOutflow_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const bool compressible = domain->isMaterialCompressible();
    const scalar comp = compressible ? 1.0 : 0.0;

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
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> umBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_umBip = &umBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_psi;
    std::vector<scalar> ws_alpha;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef(phaseIndex_).sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

    const auto* psiSTKFieldPtr =
        compressible ? model_->psiRef().stkFieldPtr() : nullptr;

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerSide);
        ws_psi.resize(nodesPerSide);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_psi = &ws_psi[0];
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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

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
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
                p_psi[ni] = compressible
                                ? *stk::mesh::field_data(*psiSTKFieldPtr, node)
                                : 0;

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                }

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

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const scalar tmDot =
                    (stk::mesh::field_data(mDotSideSTKFieldRef, side))[ip];

                if (rfflag[ip] == 0)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // zero out vector quantities; form aMag
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] = 0.0;
                        p_uBip[j] = 0.0;
                        p_umBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar rhoBip = 0;
                    scalar psiBip = 0;
                    scalar alphaBip = 0;
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];
                        const scalar r_coord =
                            p_coordinate_face_shape_function[offSetSF_face +
                                                             ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        psiBip += r_vel * p_psi[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_uBip[j] += r_vel * p_U[icNdim + j];
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_coord * p_coordinates[inn * SPATIAL_DIM + j];
                        }
                    }

                    // lhs for pressure system
                    label rowR = nearestNode * nodesPerElement;

                    //================================
                    // All-in-all: form mDot rho*uj*Aj
                    //================================
                    scalar mDot = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        mDot += rhoBip * p_uBip[j] * axj;
                    }

                    //================================
                    // Compressibility contribution at open boundary
                    // Newton-Raphson: ∂(ρU)/∂p = U * ∂ρ/∂p = U * ψ
                    // LHS coefficient = mDot / ρ_bip * ψ_bip
                    // where U is extrapolated from interior, ρ_bip and ψ_bip
                    // are interpolated to boundary ip
                    //================================
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        p_lhs[rowR + inn] += tmDot * r_vel * psiBip / rhoBip *
                                             comp * alphaBip / densityScale;
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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
                else
                {
                    // slip-wall

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
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_vel * p_coordinates[inn * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho*uj*Aj = 0
                    scalar mDot = 0.0;

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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::assembleElemTermsBoundaryOpening_(
    const domain* domain,
    const boundary* boundary,
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
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> umBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> duRhsBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_umBip = &umBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_duRhsBip = &duRhsBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_duRhs;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_bcMultiplier;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& nodalSidePSTKFieldRef =
        model_->pRef().nodeSideFieldRef().stkFieldRef();
    const auto& gradPSTKFieldRef = model_->pRef().gradRef().stkFieldRef();

    // Pressure diffusivity: duTilde for LHS (SIMPLEC), du for RHS
    const bool consistent = model_->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    const auto& duLhsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        consistent ? flowModel::duTilde_ID : flowModel::du_ID);
    const auto& duRhsSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get body force fields for buoyancy pressure stabilization
    const auto& FSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto& FOrigSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = model_->pRef().isGradientShifted();

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
        const scalar f = 1.0 / static_cast<scalar>(nodesPerSide);

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_duRhs.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_duRhs = &ws_duRhs[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_bcMultiplier[ni] = 1.0;

                // gather vectors
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
                p_p[ic] = *stk::mesh::field_data(nodalSidePSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du =
                    stk::mesh::field_data(duLhsSTKFieldRef, node);
                const scalar* duRhs =
                    stk::mesh::field_data(duRhsSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                    p_duRhs[offSet + j] = duRhs[j];
                }

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = Um[j];
                    }
                }
                else
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = 0.0;
                    }
                }

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(FSTKFieldRef, node);
                const scalar* FOrig =
                    stk::mesh::field_data(FOrigSTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }
            }

            // compute dndx for residual
            scalar scs_error = 0.0;
            if (isPGradientShifted)
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

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // zero out vector quantities; form aMag
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordBip[j] = 0.0;
                    p_uBip[j] = 0.0;
                    p_umBip[j] = 0.0;
                    p_GpdxBip[j] = 0.0;
                    p_dpdxBip[j] = 0.0;
                    p_duBip[j] = 0.0;
                    p_duRhsBip[j] = 0.0;
                    p_FBip[j] = 0.0;
                    p_FOrigBip[j] = 0.0;
                }

                // interpolate to bip
                scalar rhoBip = 0;
                scalar alphaBip = 0;
                const label offSetSF_face = ip * nodesPerSide;
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
                        // use velocity shape functions
                        p_uBip[j] += r_vel * p_U[icNdim + j];
                        p_duBip[j] += r_vel * p_du[icNdim + j];
                        p_duRhsBip[j] += r_vel * p_duRhs[icNdim + j];
                        p_umBip[j] += r_vel * p_Um[icNdim + j];

                        // use coordinates
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                        // face-centre average for original body force
                        p_FOrigBip[j] += r_vel * p_FOrig[icNdim + j];

                        // arithmetic interpolation
                        p_GpdxBip[j] += f * p_Gpdx[icNdim + j];

                        // interpolate redistributed body force using shape
                        // functions
                        p_FBip[j] += f * p_F[icNdim + j];
                    }
                }

                // form dpdxBip
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    const scalar p = p_p[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_dpdxBip[j] += p_dndx[offSetDnDx + j] * p;
                    }
                }

                // lhs for pressure system
                label rowR = nearestNode * nodesPerElement;

                //================================
                // laplacian: -rhoip*Dip*Gpn.Sip
                //================================

                // element-based gradient: here we add contributions from
                // all nodes. However, due to nature of boundary
                // (dirichlet), boundary nodes must not contribute to the
                // lhs, so we multiply by 0 the lhs when the current node
                // belongs to the boundary
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;

                    scalar lhsFac = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        lhsFac += -rhoBip * p_duBip[j] *
                                  p_dndx[offSetDnDx + j] *
                                  areaVec[ip * SPATIAL_DIM + j];
                    }
                    p_lhs[rowR + ic] +=
                        lhsFac * p_bcMultiplier[ic] * alphaBip / densityScale;
                }

                //================================
                // All-in-all: form mDot rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj
                // + rho*du*(FOrigj - Fj)*Aj (buoyancy
                // stabilization)
                //================================
                scalar mDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                    // divergence + pressure Rhie-Chow (D for RHS)
                    mDot += (rhoBip * p_uBip[j] -
                             rhoBip * p_duRhsBip[j] *
                                 (p_dpdxBip[j] - p_GpdxBip[j])) *
                            axj;

                    // // buoyancy stabilization: +rho*D*(F_orig - F)·S
                    // mDot += rhoBip * p_duRhsBip[j] *
                    //         (p_FOrigBip[j] - p_FBip[j]) * axj;
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

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryInletSpecifiedPressure_(const domain* domain,
                                                     const boundary* boundary,
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
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> umBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> duRhsBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_umBip = &umBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_duRhsBip = &duRhsBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_duRhs;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_bcMultiplier;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& nodalSidePSTKFieldRef =
        model_->pRef().nodeSideFieldRef().stkFieldRef();
    const auto& gradPSTKFieldRef = model_->pRef().gradRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Pressure diffusivity: duTilde for LHS (SIMPLEC), du for RHS
    const bool consistent = model_->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    const auto& duLhsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        consistent ? flowModel::duTilde_ID : flowModel::du_ID);
    const auto& duRhsSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get body force fields for buoyancy pressure stabilization
    const auto& FSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto& FOrigSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

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
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = model_->pRef().isGradientShifted();

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
        const scalar f = 1.0 / static_cast<scalar>(nodesPerSide);

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_duRhs.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_duRhs = &ws_duRhs[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_bcMultiplier[ni] = 1.0;

                // gather vectors
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
                p_p[ic] = *stk::mesh::field_data(nodalSidePSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du =
                    stk::mesh::field_data(duLhsSTKFieldRef, node);
                const scalar* duRhs =
                    stk::mesh::field_data(duRhsSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                    p_duRhs[offSet + j] = duRhs[j];
                }

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = Um[j];
                    }
                }
                else
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = 0.0;
                    }
                }

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(FSTKFieldRef, node);
                const scalar* FOrig =
                    stk::mesh::field_data(FOrigSTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }
            }

            // compute dndx for residual
            scalar scs_error = 0.0;
            if (isPGradientShifted)
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

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                if (rfflag[ip] == 0)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // zero out vector quantities; form aMag
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] = 0.0;
                        p_uBip[j] = 0.0;
                        p_umBip[j] = 0.0;
                        p_GpdxBip[j] = 0.0;
                        p_dpdxBip[j] = 0.0;
                        p_duBip[j] = 0.0;
                        p_duRhsBip[j] = 0.0;
                        p_FBip[j] = 0.0;
                        p_FOrigBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar rhoBip = 0;
                    scalar alphaBip = 0;
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];
                        const scalar r_coord =
                            p_coordinate_face_shape_function[offSetSF_face +
                                                             ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_uBip[j] += r_vel * p_U[icNdim + j];
                            p_umBip[j] += r_vel * p_Um[icNdim + j];
                            p_duBip[j] += r_vel * p_du[icNdim + j];
                            p_duRhsBip[j] += r_vel * p_duRhs[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                            // face-centre average for original body force
                            p_FOrigBip[j] += r_vel * p_FOrig[icNdim + j];

                            // arithmetic interpolation
                            p_GpdxBip[j] += f * p_Gpdx[icNdim + j];

                            // interpolate redistributed body force using shape
                            // functions
                            p_FBip[j] += f * p_F[icNdim + j];
                        }
                    }

                    // form dpdxBip
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;
                        const scalar p = p_p[ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dpdxBip[j] += p_dndx[offSetDnDx + j] * p;
                        }
                    }

                    // lhs for pressure system
                    label rowR = nearestNode * nodesPerElement;

                    //================================
                    // laplacian: -rhoip*Dip*Gpn.Sip
                    //================================

                    // element-based gradient: here we add contributions from
                    // all nodes. However, due to nature of boundary
                    // (dirichlet), boundary nodes must not contribute to the
                    // lhs, so we multiply by 0 the lhs when the current node
                    // belongs to the boundary
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;

                        scalar lhsFac = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            lhsFac += -rhoBip * p_duBip[j] *
                                      p_dndx[offSetDnDx + j] *
                                      areaVec[ip * SPATIAL_DIM + j];
                        }
                        p_lhs[rowR + ic] += lhsFac * p_bcMultiplier[ic] *
                                            alphaBip / densityScale;
                    }

                    //================================
                    // All-in-all: form mDot rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj
                    // + rho*du*(FOrigj - Fj)*Aj (buoyancy
                    // stabilization)
                    //================================
                    scalar mDot = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                        // divergence + pressure Rhie-Chow (D for RHS)
                        mDot += (rhoBip * p_uBip[j] -
                                 rhoBip * p_duRhsBip[j] *
                                     (p_dpdxBip[j] - p_GpdxBip[j])) *
                                axj;

                        // // buoyancy stabilization: +rho*D*(F_orig - F)·S
                        // mDot += rhoBip * p_duRhsBip[j] *
                        //         (p_FOrigBip[j] - p_FBip[j]) * axj;
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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
                else
                {
                    // slip-wall

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
                    const label offSetSF_face = ip * nodesPerSide;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF_face + ic];

                        // use velocity shape functions
                        rhoBip += r_vel * p_rho[ic];
                        alphaBip += r_vel * p_alpha[ic];

                        const label icNdim = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            // use velocity shape functions
                            p_umBip[j] += r_vel * p_Um[icNdim + j];

                            // use coordinates shape functions
                            p_coordBip[j] +=
                                r_vel * p_coordinates[inn * SPATIAL_DIM + j];
                        }
                    }

                    // form mDot; rho*uj*Aj = 0
                    scalar mDot = 0.0;

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
                        mDot -=
                            rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                    }

#ifndef NDEBUG
                    // SCL check: accumulate grid flux to nearest node
                    if (sclCheckSTKFieldPtr)
                    {
                        scalar gridFlux = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            gridFlux += rhoBip * p_umBip[i] *
                                        areaVec[ip * SPATIAL_DIM + i] *
                                        alphaBip / densityScale;
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
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void bulkPressureCorrectionAssembler::
    assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_alpha;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;

    // Get fields
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef(phaseIndex_).sideFieldRef().stkFieldRef();
    const auto& alphaSTKFieldRef = model_->alphaRef(phaseIndex_).stkFieldRef();

    scalar densityScale = model_->rhoRef(phaseIndex_).scale();

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

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
        ws_alpha.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

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

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];

            // get element; its face ordinal number and populate
            // face_node_ordinals
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = face_elem_ords[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

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
            }

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
                p_alpha[ni] = *stk::mesh::field_data(alphaSTKFieldRef, node);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const scalar tmDot =
                    (stk::mesh::field_data(mDotSideSTKFieldRef, side))[ip];

                const label nearestNode = ipNodeMap[ip];

                // interpolate to bip
                scalar alphaBip = 0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF_face + ic];

                    // use velocity shape functions
                    alphaBip += r_vel * p_alpha[ic];
                }

                // residual
                p_rhs[nearestNode] -= tmDot * alphaBip / densityScale;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
