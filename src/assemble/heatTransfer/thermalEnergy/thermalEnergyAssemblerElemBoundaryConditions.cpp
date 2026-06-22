// File       : thermalEnergyAssemblerElemBoundaryConditions.cpp
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef WITH_THERMAL_TEMPERATURE

#include "thermalEnergyAssembler.h"

namespace accel
{

void thermalEnergyAssembler::assembleElemTermsBoundary_(const domain* domain,
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

void thermalEnergyAssembler::assembleElemTermsBoundaryWallFixedValue_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    if ((domain->type() == domainType::solid) ||
        (domain->type() == domainType::fluid &&
         domain->turbulence_.option_ == turbulenceOption::laminar))
    {
        // resolved-gradient Dirichlet treatment
        Base::assembleElemTermsBoundaryWallFixedValue_(domain, boundary, ctx);
    }
    else
    {
        // turbulent fluid: wall-function flux q = TWallCoeffs*(Tbc - T)
        const auto& mesh = field_broker_->meshRef();
        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // space for LHS/RHS; nodesPerSide*nodesPerSide and nodesPerSide
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // nodal fields to gather
        std::vector<scalar> ws_cp;
        std::vector<scalar> ws_T;

        // master element
        std::vector<scalar> ws_shape_function;

        // Get transport fields/side fields
        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
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

            // algorithm related; face
            ws_T.resize(nodesPerSide);
            ws_cp.resize(nodesPerSide);
            ws_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_T = &ws_T[0];
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

                    scalar cpBip = 0.0;
                    scalar TBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF_face + ic];
                        cpBip += r * p_cp[ic];
                        TBip += r * p_T[ic];
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
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void thermalEnergyAssembler::assembleElemTermsBoundaryWallMixed_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    // Robin wall q = U*(Tout - Tnw); turbulent fluid: 1/U = 1/hExt + 1/hWf
    // with hWf = rho*cp*uStar/TPlus, otherwise U = hExt
    const auto& bcRef =
        model_->TRef().boundaryConditionRef(domain->index(), boundary->index());
    const scalar Tout = (bcRef.template data<1>("value")).value()[0];
    const scalar hExt =
        (bcRef.template data<1>("transfer_coefficient")).value()[0];

    const bool turbulentWall =
        domain->type() == domainType::fluid &&
        domain->turbulence_.option_ != turbulenceOption::laminar;

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerSide*nodesPerSide and nodesPerSide
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_T;

    // master element
    std::vector<scalar> ws_shape_function;

    // Get transport fields/side fields
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* TWallCoeffsSTKFieldPtr =
        turbulentWall ? model_->TWallCoeffsRef().stkFieldPtr() : nullptr;

    // Get geometric fields
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

        // algorithm related; face
        ws_T.resize(nodesPerSide);
        ws_cp.resize(nodesPerSide);
        ws_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_T = &ws_T[0];
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

            // face data
            const scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);
            const scalar* TWallCoeffs =
                turbulentWall
                    ? stk::mesh::field_data(*TWallCoeffsSTKFieldPtr, side)
                    : nullptr;

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

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
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

                scalar cpBip = 0.0;
                scalar TBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];
                    cpBip += r * p_cp[ic];
                    TBip += r * p_T[ic];
                }

                // series combination with the log-law wall conductance
                scalar UBip = hExt;
                if (turbulentWall)
                {
                    const scalar hWf = TWallCoeffs[ip];
                    UBip = hExt * hWf / std::max(hExt + hWf, SMALL);
                }

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF_face + ic];

                    // matrix entries
                    label rowR = nearestNode * nodesPerSide;

                    p_lhs[rowR + ic] += UBip * amag * r / cpBip;
                }

                p_rhs[nearestNode] += UBip * amag * (Tout - TBip);
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalEnergyAssembler::assembleElemTermsBoundaryWallSpecifiedFlux_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // side flux field comes from TEMPERATURE (where the heat_flux BC stores
    // it), not from phi_ (enthalpy)
    const auto& sidePhiFluxSTKFieldRef =
        model_->TRef().sideFluxFieldRef().stkFieldRef();

    // geometric field
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;
        const label* faceIpNodeMap = meFC->ipNodeMap();

        const label lhsSize = nodesPerSide * nodesPerSide;
        const label rhsSize = nodesPerSide;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerSide);

        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];

        const stk::mesh::Bucket::size_type nBoundarySides = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundarySides;
             ++iSide)
        {
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            stk::mesh::Entity side = sideBucket[iSide];

            const scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);
            const scalar* phiFluxVec = stk::mesh::field_data(
                sidePhiFluxSTKFieldRef, sideBucket, iSide);

            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);
            STK_ThrowAssert(numSideNodes == nodesPerSide);

            for (label iNode = 0; iNode < numSideNodes; ++iNode)
            {
                connectedNodes[iNode] = sideNodeRels[iNode];
            }

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

                // flux value is stored into the domain: multiply by -1
                p_rhs[nearestNode] -= (-phiFluxVec[ip]) * amag;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif
