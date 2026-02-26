// File       : totalEnergyAssemblerElemBoundaryConditions.cpp
// Created    : Thu Apr 02 2025 15:40:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "totalEnergyAssembler.h"

namespace accel
{

void totalEnergyAssembler::assembleElemTermsBoundary_(const domain* domain,
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
                {
                    assembleElemTermsBoundarySymmetry_(domain, boundary, ctx);
                }
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
                            {
                                assembleElemTermsBoundaryWallZeroGradient_(
                                    domain, boundary, ctx);
                            }
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

void totalEnergyAssembler::assembleElemTermsBoundarySymmetry_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and
    // nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_vw;

    // master element
    std::vector<scalar> ws_shape_function;

    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

    // Get geometric fields
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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_vw = &ws_vw[0];

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

            // Fill gamma + override boundary nodes of phi with the
            // dirichlet values
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
            }

            // loop over side ip's
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF_face = ip * nodesPerSide;

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_vwBip[j] = 0.0;
                }

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                    }
                }

                //================================
                // Viscous Work
                //================================
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // matrix entries
                    label indexR = nearestNode;

                    const scalar axj = areaVec[faceOffSet + j];
                    p_rhs[indexR] += vwBip[j] * axj;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryWallZeroGradient_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and
    // nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_vw;

    // master element
    std::vector<scalar> ws_shape_function;

    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

    // Get geometric fields
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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_vw = &ws_vw[0];
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

            // pointer to face data
            const scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);

            //==========================================
            // gather nodal data off of element; n/a
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
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);

            // Fill gamma + override boundary nodes of phi with the
            // dirichlet values
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
            }

            // loop over side ip's
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF_face = ip * nodesPerSide;

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_vwBip[j] = 0.0;
                }

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                    }
                }

                //================================
                // Viscous Work
                //================================
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // matrix entries
                    label indexR = nearestNode;

                    const scalar axj = areaVec[faceOffSet + j];
                    p_rhs[indexR] += vwBip[j] * axj;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryWallFixedValue_(
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

        const bool includeAdv = domain->type() == domainType::fluid;
        const bool includeViscousWork =
            domain->heatTransfer_.includeViscousWork_ && includeAdv;

        // space for LHS/RHS; nodesPerElement*nodesPerElement and
        // nodesPerElement
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values
        std::vector<scalar> vwBip(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_vwBip = &vwBip[0];

        // nodal fields to gather
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_lambdaEff;
        std::vector<scalar> ws_cp;
        std::vector<scalar> ws_T;
        std::vector<scalar> ws_vw;

        std::vector<scalar> ws_bcMultiplier;

        // master element
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
        const auto& nodalSideTSTKFieldRef =
            model_->TRef().nodeSideFieldRef().stkFieldRef();
        const auto* USTKFieldPtr =
            includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
        const auto* gradUSTKFieldPtr =
            includeViscousWork ? model_->URef().gradRef().stkFieldPtr()
                               : nullptr;
        const auto* muEffSTKFieldPtr =
            includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
            ws_vw.resize(nodesPerSide * SPATIAL_DIM);

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
            scalar* p_vw = &ws_vw[0];

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
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = elemNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    p_bcMultiplier[ni] = 1.0;

                    // gather scalars
                    p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                    // gather vectors
                    scalar* coords =
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
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);

                // Fill gamma + override boundary nodes of phi with the
                // dirichlet values
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    const label ic = faceNodeOrdinals[ni];

                    // gather scalars
                    p_lambdaEff[ni] =
                        *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                    p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                    p_T[ic] =
                        *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                    // set 0 the boundary nodes
                    p_bcMultiplier[ic] = 0.0;

                    if (includeViscousWork)
                    {
                        const scalar muEff =
                            *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                        const scalar* U =
                            stk::mesh::field_data(*USTKFieldPtr, node);
                        const scalar* dudx =
                            stk::mesh::field_data(*gradUSTKFieldPtr, node);

                        // calculate divergence of velocity
                        scalar divU = 0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            divU += dudx[i * SPATIAL_DIM + i];
                        }

                        // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                        // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vw[ni * SPATIAL_DIM + i] =
                                -2.0 / 3.0 * muEff * divU * U[i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_vw[ni * SPATIAL_DIM + i] +=
                                    muEff *
                                    (dudx[i * SPATIAL_DIM + j] +
                                     dudx[j * SPATIAL_DIM + i]) *
                                    U[j];
                            }
                        }
                    }
                    else
                    {
                        // set viscous work to 0
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vw[ni * SPATIAL_DIM + i] = 0;
                        }
                    }
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

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
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

                            const scalar T = p_T[ic];

                            // -lambda*dT/dxj*Aj
                            scalar lhsfac = -lambdaEffBip * dndxj * axj;
                            p_lhs[rowR + ic] +=
                                lhsfac * p_bcMultiplier[ic] / cpBip;
                            p_rhs[indexR] -= lhsfac * T;
                        }
                    }

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // matrix entries
                        label indexR = nearestNode;

                        const scalar axj = areaVec[faceOffSet + j];
                        p_rhs[indexR] += vwBip[j] * axj;
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

        const bool includeAdv = domain->type() == domainType::fluid;
        const bool includeViscousWork =
            domain->heatTransfer_.includeViscousWork_ && includeAdv;

        // space for LHS/RHS; nodesPerElement*nodesPerElement and
        // nodesPerElement
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values
        std::vector<scalar> vwBip(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_vwBip = &vwBip[0];

        // nodal fields to gather
        std::vector<scalar> ws_cp;
        std::vector<scalar> ws_T;
        std::vector<scalar> ws_vw;

        // master element
        std::vector<scalar> ws_shape_function;

        // Get transport fields/side fields
        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
        const auto& sideTSTKFieldRef =
            model_->TRef().sideFieldRef().stkFieldRef();
        const auto* TWallCoeffsSTKFieldPtr =
            model_->TWallCoeffsRef().stkFieldPtr();
        const auto* USTKFieldPtr =
            includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
        const auto* gradUSTKFieldPtr =
            includeViscousWork ? model_->URef().gradRef().stkFieldPtr()
                               : nullptr;
        const auto* muEffSTKFieldPtr =
            includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
            ws_cp.resize(nodesPerSide);
            ws_vw.resize(nodesPerSide * SPATIAL_DIM);
            ws_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_T = &ws_T[0];
            scalar* p_cp = &ws_cp[0];
            scalar* p_vw = &ws_vw[0];
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
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                    p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                    if (includeViscousWork)
                    {
                        const scalar muEff =
                            *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                        const scalar* U =
                            stk::mesh::field_data(*USTKFieldPtr, node);
                        const scalar* dudx =
                            stk::mesh::field_data(*gradUSTKFieldPtr, node);

                        // calculate divergence of velocity
                        scalar divU = 0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            divU += dudx[i * SPATIAL_DIM + i];
                        }

                        // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                        // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vw[ni * SPATIAL_DIM + i] =
                                -2.0 / 3.0 * muEff * divU * U[i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_vw[ni * SPATIAL_DIM + i] +=
                                    muEff *
                                    (dudx[i * SPATIAL_DIM + j] +
                                     dudx[j * SPATIAL_DIM + i]) *
                                    U[j];
                            }
                        }
                    }
                    else
                    {
                        // set viscous work to 0
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vw[ni * SPATIAL_DIM + i] = 0;
                        }
                    }
                }

                // loop over side ip's
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = faceIpNodeMap[ip];

                    const label offSetSF_face = ip * nodesPerSide;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

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

                    scalar cpBip = 0.0;
                    scalar TBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];
                        cpBip += r * p_cp[ic];
                        TBip += r * p_T[ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    scalar TWallCoeffBip = (stk::mesh::field_data(
                        *TWallCoeffsSTKFieldPtr, side))[ip];
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        // matrix entries
                        label rowR = nearestNode * nodesPerSide;

                        p_lhs[rowR + ic] += TWallCoeffBip * amag * r / cpBip;
                    }

                    p_rhs[nearestNode] +=
                        TWallCoeffBip * amag * (Tbc[ip] - TBip);

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[nearestNode] += vwBip[j] * axj;
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryWallSpecifiedFlux_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_vw;

    // master element
    std::vector<scalar> ws_shape_function;

    // Get transport fields/side fields
    const auto& sideTFluxSTKFieldRef =
        model_->TRef().sideFluxFieldRef().stkFieldRef();
    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_vw = &ws_vw[0];
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
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
            }

            // loop over side ip's
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = faceIpNodeMap[ip];

                const label offSetSF_face = ip * nodesPerSide;

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_vwBip[j] = 0.0;
                }

                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                    }
                }

                //================================
                // Diffusion: only diffusion
                //================================

                // matrix entries
                label indexR = nearestNode;

                // flux value is stored into domain: multiply by -1
                p_rhs[indexR] -= (-TFluxVec[ip]) * amag;

                //================================
                // Viscous Work
                //================================
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                    p_rhs[indexR] += vwBip[j] * axj;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryWallMixed_(
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

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_bcMultiplier;
    std::vector<scalar> ws_vw;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get transport fields/side fields
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideTFluxSTKFieldRef =
        model_->TRef().sideFluxFieldRef().stkFieldRef();
    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);
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
        scalar* p_vw = &ws_vw[0];
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
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                p_bcMultiplier[ni] = 1.0;

                // gather scalars
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
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

            // Fill gamma + overwrite phi values at boundary with dirichlet
            // values
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
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

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_vwBip[j] = 0.0;
                }

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
                scalar cpBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];
                    lambdaEffBip += r * p_lambdaEff[ic];
                    cpBip += r * p_cp[ic];

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                    }
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
                        p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic] *
                                            fixedValueFlag[0] / cpBip;
                        p_rhs[indexR] -= lhsfac * phi * fixedValueFlag[0];
                    }
                }

                // 2) specified flux

                // matrix entries
                label indexR = nearestNode;

                // flux value is stored into domain: multiply by -1
                p_rhs[indexR] -=
                    (-TFluxVec[ip]) * amag * (1.0 - fixedValueFlag[0]);

                //================================
                // Viscous Work
                //================================
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                    p_rhs[indexR] += vwBip[j] * axj;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryInletFixedValue_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_vw;

    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideH0STKFieldRef =
        model_->h0Ref().sideFieldRef().stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();
    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);

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
        scalar* p_vw = &ws_vw[0];
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
            const scalar* h0bcVec =
                stk::mesh::field_data(sideH0STKFieldRef, sideBucket, iSide);

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
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                p_bcMultiplier[ni] = 1.0;

                // gather scalars
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
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

            // Fill gamma + override boundary nodes of phi with the dirichlet
            // values
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                const label ic = faceNodeOrdinals[ni];

                // gather scalars
                p_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
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

                    const label nearestNode = ipNodeMap[ip];

                    const label offSetSF_face = ip * nodesPerSide;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[nearestNode] += vwBip[j] * axj;
                    }
                }
                else
                {
                    const label nearestNode = ipNodeMap[ip];

                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    //================================
                    // Advection: entrainment; advect
                    // in from specified value
                    //================================
                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const scalar aflux = tmDot * h0bcVec[ip];
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

                            const scalar T = p_T[ic];

                            // -lambda*dT/dxj*Aj
                            scalar lhsfac = -lambdaEffBip * dndxj * axj;
                            p_lhs[rowR + ic] +=
                                lhsfac * p_bcMultiplier[ic] / cpBip;
                            p_rhs[indexR] -= lhsfac * T;
                        }
                    }

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // matrix entries
                        label indexR = nearestNode;

                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[indexR] += vwBip[j] * axj;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryOutletZeroGradient_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_vw;

    // master element
    std::vector<scalar> ws_shape_function;

    const auto& h0STKFieldRef = phi_->stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();
    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_vw = &ws_vw[0];
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

            // Fill scratch vectors
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
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

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    // right node; right is on the face
                    stk::mesh::Entity nodeR = elemNodeRels[nearestNode];

                    const scalar h0IpUpw =
                        *stk::mesh::field_data(h0STKFieldRef, nodeR);

                    const label indexR = nearestNode;
                    const label rowR = indexR * nodesPerElement;

                    // total advection
                    const scalar aflux = tmDot * h0IpUpw;

                    p_rhs[indexR] -= aflux;

                    // upwind lhs
                    p_lhs[rowR + nearestNode] += tmDot;

                    // no diffusion .. zero-gradient

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[indexR] += vwBip[j] * axj;
                    }
                }
                else
                {
                    // reversed flow .. treat as slip wall, that is,
                    // zero-gradient

                    const label nearestNode = ipNodeMap[ip];

                    const label offSetSF_face = ip * nodesPerSide;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // matrix entries
                        label indexR = nearestNode;

                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[indexR] += vwBip[j] * axj;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleElemTermsBoundaryOpening_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;

    // space for LHS/RHS; nodesPerElement*nodesPerElement and nodesPerElement
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values
    std::vector<scalar> vwBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_vwBip = &vwBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_vw;

    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    const auto& h0STKFieldRef = phi_->stkFieldRef();
    const auto& nodalSideTSTKFieldRef =
        model_->TRef().nodeSideFieldRef().stkFieldRef();
    const auto& sideH0STKFieldRef =
        model_->h0Ref().sideFieldRef().stkFieldRef();
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* reversalFlagSTKFieldPtr =
        model_->URef().reversalFlagRef().stkFieldPtr();
    const auto* USTKFieldPtr =
        includeViscousWork ? model_->URef().stkFieldPtr() : nullptr;
    const auto* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        ws_vw.resize(nodesPerSide * SPATIAL_DIM);

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
        scalar* p_vw = &ws_vw[0];

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
            const scalar* h0bcVec =
                stk::mesh::field_data(sideH0STKFieldRef, sideBucket, iSide);

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
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_bcMultiplier[ni] = 1.0;

                // gather scalars
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
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

            // Fill gamma + override boundary nodes of phi with the dirichlet
            // values
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                const label ic = faceNodeOrdinals[ni];

                // gather scalars
                p_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);

                p_T[ic] = *stk::mesh::field_data(nodalSideTSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                if (includeViscousWork)
                {
                    const scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    const scalar* dudx =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);

                    // calculate divergence of velocity
                    scalar divU = 0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += dudx[i * SPATIAL_DIM + i];
                    }

                    // calculate viscous work: VW_i = Σⱼ τ_ji U_j
                    // where τ_ji = μ(∂U_j/∂x_i + ∂U_i/∂x_j - 2/3 δ_ji ∇·U)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] =
                            -2.0 / 3.0 * muEff * divU * U[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vw[ni * SPATIAL_DIM + i] +=
                                muEff *
                                (dudx[i * SPATIAL_DIM + j] +
                                 dudx[j * SPATIAL_DIM + i]) *
                                U[j];
                        }
                    }
                }
                else
                {
                    // set viscous work to 0
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vw[ni * SPATIAL_DIM + i] = 0;
                    }
                }
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

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    // right node; right is on the face
                    stk::mesh::Entity nodeR = elemNodeRels[nearestNode];

                    const scalar h0IpUpw =
                        *stk::mesh::field_data(h0STKFieldRef, nodeR);

                    const label indexR = nearestNode;
                    const label rowR = indexR * nodesPerElement;

                    // total advection
                    const scalar aflux = tmDot * h0IpUpw;

                    p_rhs[indexR] -= aflux;

                    // upwind lhs
                    p_lhs[rowR + nearestNode] += tmDot;

                    // no diffusion .. zero-gradient

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[indexR] += vwBip[j] * axj;
                    }
                }
                else
                {
                    const label nearestNode = ipNodeMap[ip];

                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF_face = ip * nodesPerSide;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_vwBip[j] = 0.0;
                    }

                    scalar lambdaEffBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];

                        lambdaEffBip += r * p_lambdaEff[ic];
                        cpBip += r * p_cp[ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_vwBip[i] += r * p_vw[ic * SPATIAL_DIM + i];
                        }
                    }

                    //================================
                    // Advection: entrainment; advect
                    // in from specified value
                    //================================
                    const scalar tmDot = (stk::mesh::field_data(
                        *mDotSideSTKFieldPtr_, side))[ip];

                    const scalar aflux = tmDot * h0bcVec[ip];
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
                            p_lhs[rowR + ic] +=
                                lhsfac * p_bcMultiplier[ic] / cpBip;
                            p_rhs[indexR] -= lhsfac * Ti;
                        }
                    }

                    //================================
                    // Viscous Work
                    //================================
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                        p_rhs[nearestNode] += vwBip[j] * axj;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
