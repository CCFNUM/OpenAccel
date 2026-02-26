// File       : thermalTemperatureAssemblerElemBoundaryConditions.cpp
// Created    : Thu Apr 14 2024 8:36:38 (+0100)
// Author     : Fabian Wermelinger
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef WITH_THERMAL_TEMPERATURE

#include "thermalTemperatureAssembler.h"

namespace accel
{

void thermalTemperatureAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
{
    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            model_->TRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                assembleElemTermsBoundaryWallFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::specifiedFlux:
                            {
                                assembleElemTermsBoundaryWallSpecifiedFlux_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            break;

                        case boundaryConditionType::mixed:
                            {
                                assembleElemTermsBoundaryWallMixed_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                assembleElemTermsBoundaryInletFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                assembleElemTermsBoundaryOutletZeroGradient_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("boundary condition invalid");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    assembleElemTermsBoundaryOpening_(domain, boundary, ctx);
                }
                break;

            default:
                break;
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryWallFixedValue_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    if ((domain->type() == domainType::solid) ||
        (domain->type() == domainType::fluid &&
         domain->turbulence_.option_ == turbulenceOption::laminar))
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

        // nodal fields to gather
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_lambdaEff;
        std::vector<scalar> ws_cp;
        std::vector<scalar> ws_T;
        std::vector<scalar> ws_bcMultiplier;

        // master element
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
        const auto& nodalSideTSTKFieldRef =
            model_->TRef().nodeSideFieldRef().stkFieldRef();

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
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for field?
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
            ws_T.resize(nodesPerElement);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_lambdaEff.resize(nodesPerSide);
            ws_cp.resize(nodesPerSide);
            ws_bcMultiplier.resize(nodesPerElement);
            ws_shape_function.resize(numScsBip * nodesPerSide);
            ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
            ws_det_j.resize(numScsBip);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_T = &ws_T[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_lambdaEff = &ws_lambdaEff[0];
            scalar* p_cp = &ws_cp[0];
            scalar* p_bcMultiplier = &ws_bcMultiplier[0];
            scalar* p_shape_function = &ws_shape_function[0];
            scalar* p_dndx = &ws_dndx[0];

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

                // pointer to face data
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
                const label faceOrdinal =
                    bulkData.begin_element_ordinals(side)[0];

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

                    p_bcMultiplier[iNode] = 1.0;

                    // gather scalars
                    p_T[iNode] = *stk::mesh::field_data(TSTKFieldRef, node);

                    // gather vectors
                    scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label offSet = iNode * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordinates[offSet + j] = coords[j];
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

                // Fill gamma + override boundary nodes of phi with the
                // dirichlet values
                for (label iNode = 0; iNode < numSideNodes; ++iNode)
                {
                    stk::mesh::Entity node = sideNodeRels[iNode];

                    const label ic = faceNodeOrdinals[iNode];

                    // gather scalars
                    p_lambdaEff[iNode] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                    p_cp[iNode] = *stk::mesh::field_data(cpSTKFieldRef, node);
                    p_T[ic] =
                        *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                    // set 0 the boundary nodes
                    p_bcMultiplier[ic] = 0.0;
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
                    const label offSetSF_face = ip * nodesPerSide;

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];
                    }

                    //================================
                    // Diffusion
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

                            // matrix entries
                            label indexR = nearestNode;
                            label rowR = indexR * nodesPerElement;

                            const scalar Ti = p_T[ic];

                            // -lambda*dT/dxj*Aj
                            scalar lhsfac = -lambdaEffBip * dndxj * axj;
                            p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic];
                            p_rhs[indexR] -= lhsfac * Ti;
                        }
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
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

        // space for LHS/RHS; nodesPerElement*nodesPerElement and
        // nodesPerElement
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // nodal fields to gather
        std::vector<scalar> ws_T;

        // master element
        std::vector<scalar> ws_shape_function;

        // Get transport fields/side fields
        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& sideTSTKFieldRef =
            model_->TRef().sideFieldRef().stkFieldRef();
        const auto* TWallCoeffsSTKFieldPtr =
            model_->TWallCoeffsRef().stkFieldPtr();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for field?
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
            const label lhsSize = nodesPerSide * nodesPerSide;
            const label rhsSize = nodesPerSide;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerSide);

            // algorithm related; element
            ws_T.resize(nodesPerSide);
            ws_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_T = &ws_T[0];
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
                const scalar* Tbc =
                    stk::mesh::field_data(sideTSTKFieldRef, side);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);

                // Fill scratch vectors
                for (label ic = 0; ic < numSideNodes; ++ic)
                {
                    stk::mesh::Entity node = sideNodeRels[ic];

                    // set connected nodes
                    connectedNodes[ic] = node;

                    // gather scalars
                    p_T[ic] = *stk::mesh::field_data(TSTKFieldRef, node);
                }

                // loop over side ip's
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = faceIpNodeMap[ip];

                    const label offSetSF_face = ip * nodesPerSide;

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

                    scalar TBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];
                        TBip += r * p_T[ic];
                    }

                    scalar TWallCoeffBip = (stk::mesh::field_data(
                        *TWallCoeffsSTKFieldPtr, side))[ip];
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        // matrix entries
                        label rowR = nearestNode * nodesPerSide;

                        p_lhs[rowR + ic] += TWallCoeffBip * amag * r;
                    }

                    p_rhs[nearestNode] +=
                        TWallCoeffBip * amag * (Tbc[ip] - TBip);
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryWallSpecifiedFlux_(
    const domain* domain,
    const boundary* boundary,
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

    // master element
    std::vector<scalar> ws_shape_function;

    // Get transport fields/side fields
    const auto& sideTFluxSTKFieldRef =
        model_->TRef().sideFluxFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
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

        // mapping from ip to nodes for this ordinal; face perspective (use with
        // face_node_relations)
        const label* faceIpNodeMap = meFC->ipNodeMap();

        // resize some things; matrix related
        const label lhsSize = nodesPerSide * nodesPerSide;
        const label rhsSize = nodesPerSide;
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
            const scalar* TFluxVec =
                stk::mesh::field_data(sideTFluxSTKFieldRef, sideBucket, iSide);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
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

                // flux value is stored into domain: multiply by -1
                p_rhs[nearestNode] -= (-TFluxVec[ip]) * amag;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryWallMixed_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const scalar* fixedValueFlag =
        ((phi_->boundaryConditionRef(domain->index(), boundary->index()))
             .template data<1>("fixed_value_flag"))
            .value();

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
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get transport fields/side fields
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideTFluxSTKFieldRef =
        model_->TRef().sideFluxFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
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
        ws_T.resize(nodesPerElement);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_lambdaEff.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_T = &ws_T[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_lambdaEff = &ws_lambdaEff[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isShifted)
        {
            meFC->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_shape_function[0]);
        }

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
            const scalar* TFluxVec =
                stk::mesh::field_data(sideTFluxSTKFieldRef, sideBucket, iSide);

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

                p_bcMultiplier[iNode] = 1.0;

                // gather scalars
                p_T[iNode] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = iNode * SPATIAL_DIM;
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

            // Fill gamma + overwrite phi values at boundary with dirichlet
            // values
            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                const label ic = faceNodeOrdinals[iNode];

                stk::mesh::Entity node = sideNodeRels[iNode];

                // gather scalars
                p_lambdaEff[iNode] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);

                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;
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
                const label offSetSF_face = ip * nodesPerSide;

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

                // 1) specified value
                scalar lambdaEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];
                    lambdaEffBip += r * p_lambdaEff[ic];
                }

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

                        const scalar phi = p_T[ic];

                        // -lambda*dT/dxj*Aj
                        scalar lhsfac = -lambdaEffBip * dndxj * axj;
                        p_lhs[rowR + ic] +=
                            lhsfac * p_bcMultiplier[ic] * fixedValueFlag[0];
                        p_rhs[indexR] -= lhsfac * phi * fixedValueFlag[0];
                    }
                }

                // 2) specified flux

                // matrix entries
                label indexR = nearestNode;

                // flux value is stored into domain: multiply by -1
                p_rhs[indexR] -=
                    (-TFluxVec[ip]) * amag * (1.0 - fixedValueFlag[0]);
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryInletFixedValue_(
    const domain* domain,
    const boundary* boundary,
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
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideTSTKFieldRef = model_->TRef().sideFieldRef().stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
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
        ws_T.resize(nodesPerElement);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_lambdaEff.resize(nodesPerSide);
        ws_cp.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_T = &ws_T[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_lambdaEff = &ws_lambdaEff[0];
        scalar* p_cp = &ws_cp[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isShifted)
        {
            meFC->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_shape_function[0]);
        }

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

            // pointer to face data
            const scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);
            const scalar* Tbc =
                stk::mesh::field_data(sideTSTKFieldRef, sideBucket, iSide);

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

                p_bcMultiplier[iNode] = 1.0;

                // gather scalars
                p_T[iNode] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = iNode * SPATIAL_DIM;
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

            // Fill gamma + override boundary nodes of phi with the dirichlet
            // values
            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                stk::mesh::Entity node = sideNodeRels[iNode];

                const label ic = faceNodeOrdinals[iNode];

                // gather scalars
                p_lambdaEff[iNode] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[iNode] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;
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
                const label trfflag =
                    (stk::mesh::field_data(*reversalFlagSTKFieldPtr, side))[ip];

                if (trfflag == 1)
                {
                    // reversed flow .. treat as slip wall, that is, no
                    // advection and no diffusion
                }
                else
                {
                    const label nearestNode = ipNodeMap[ip];

                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];
                    }

                    //================================
                    // Advection: entrainment; advect
                    // in from specified value
                    //================================
                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const scalar aflux = tmDot * cpBip * Tbc[ip];
                    p_rhs[nearestNode] -= aflux;

                    //================================
                    // Diffusion
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

                            // matrix entries
                            label indexR = nearestNode;
                            label rowR = indexR * nodesPerElement;

                            const scalar Ti = p_T[ic];

                            // -lambda*dT/dxj*Aj
                            scalar lhsfac = -lambdaEffBip * dndxj * axj;
                            p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic];
                            p_rhs[indexR] -= lhsfac * Ti;
                        }
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryOutletZeroGradient_(
    const domain* domain,
    const boundary* boundary,
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
    std::vector<scalar> ws_cp;

    // master element
    std::vector<scalar> ws_shape_function;

    const auto& TSTKFieldRef = phi_->stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
    const bool isShifted = phi_->isShifted();

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
        ws_cp.resize(nodesPerSide);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_cp = &ws_cp[0];
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
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);

            // Fill gamma + override boundary nodes of phi with the dirichlet
            // values
            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                stk::mesh::Entity node = sideNodeRels[iNode];

                p_cp[iNode] = *stk::mesh::field_data(cpSTKFieldRef, node);
            }

            // loop over side ip's
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label trfflag =
                    (stk::mesh::field_data(*reversalFlagSTKFieldPtr, side))[ip];

                if (trfflag == 0)
                {
                    //================================
                    // Advection: only
                    //================================

                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const label nearestNode = ipNodeMap[ip];
                    const label offSetSF_face = ip * nodesPerSide;

                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        cpBip += r * p_cp[ic];
                    }

                    // right node; right is on the face
                    stk::mesh::Entity nodeR = elemNodeRels[nearestNode];

                    const scalar TIpUpw =
                        *stk::mesh::field_data(TSTKFieldRef, nodeR);

                    const label rowR = nearestNode * nodesPerElement;

                    // total advection
                    const scalar aflux = tmDot * cpBip * TIpUpw;

                    p_rhs[nearestNode] -= aflux;

                    // upwind lhs
                    p_lhs[rowR + nearestNode] += tmDot * cpBip;

                    // no diffusion .. zero-gradient
                }
                else
                {
                    // reversed flow .. treat as slip wall, that is,
                    // zero-gradient
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleElemTermsBoundaryOpening_(
    const domain* domain,
    const boundary* boundary,
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
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_Cp;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideTSTKFieldRef = model_->TRef().sideFieldRef().stkFieldRef();
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
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
        ws_T.resize(nodesPerElement);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_lambdaEff.resize(nodesPerSide);
        ws_Cp.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_T = &ws_T[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_lambdaEff = &ws_lambdaEff[0];
        scalar* p_cp = &ws_Cp[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isShifted)
        {
            meFC->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_shape_function[0]);
        }

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

            // pointer to face data
            const scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);
            const scalar* Tbc =
                stk::mesh::field_data(sideTSTKFieldRef, sideBucket, iSide);

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

                // gather scalars
                p_T[iNode] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = iNode * SPATIAL_DIM;
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

            // Fill gamma + override boundary nodes of phi with the dirichlet
            // values
            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                stk::mesh::Entity node = sideNodeRels[iNode];

                const label ic = faceNodeOrdinals[iNode];

                // gather scalars
                p_lambdaEff[iNode] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[iNode] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;
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
                const label trfflag =
                    (stk::mesh::field_data(*reversalFlagSTKFieldPtr, side))[ip];

                if (trfflag == 0)
                {
                    //================================
                    // Advection: only
                    //================================

                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const label nearestNode = ipNodeMap[ip];
                    const label offSetSF_face = ip * nodesPerSide;

                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        cpBip += r * p_cp[ic];
                    }

                    // right node; right is on the face
                    stk::mesh::Entity nodeR = elemNodeRels[nearestNode];

                    const scalar TIpUpw =
                        *stk::mesh::field_data(TSTKFieldRef, nodeR);

                    const label indexR = nearestNode;
                    const label rowR = indexR * nodesPerElement;

                    // total advection
                    const scalar aflux = tmDot * cpBip * TIpUpw;

                    p_rhs[indexR] -= aflux;

                    // upwind lhs
                    p_lhs[rowR + nearestNode] += tmDot * cpBip;

                    // no diffusion .. zero-gradient
                }
                else
                {
                    const label nearestNode = ipNodeMap[ip];

                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];
                    }

                    //================================
                    // Advection: entrainment; advect
                    // in from specified value
                    //================================
                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const scalar aflux = tmDot * cpBip * Tbc[ip];
                    p_rhs[nearestNode] -= aflux;

                    //================================
                    // Diffusion
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

                            // matrix entries
                            label indexR = nearestNode;
                            label rowR = indexR * nodesPerElement;

                            const scalar Ti = p_T[ic];

                            // -lambda*dT/dxj*Aj
                            scalar lhsfac = -lambdaEffBip * dndxj * axj;
                            p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic];
                            p_rhs[indexR] -= lhsfac * Ti;
                        }
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif
