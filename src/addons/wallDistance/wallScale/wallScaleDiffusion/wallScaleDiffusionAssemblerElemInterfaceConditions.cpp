// File       : wallScaleDiffusionAssemblerElemInterfaceConditions.cpp
// Created    : Thu Jun 19 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "wallScaleDiffusionAssembler.h"

namespace accel
{

void wallScaleDiffusionAssembler::assembleElemTermsInterfaces_(
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
                // treat as a no-slip wall
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

void wallScaleDiffusionAssembler::assembleElemTermsInterfaceSideNoSlipWall_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
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
    std::vector<scalar> ws_yScale;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get transport fields/side fields
    const auto& yScaleSTKFieldRef = phi_->stkFieldRef();
    const auto& nodalSideYScaleSTKFieldRef =
        phi_->nodeSideFieldRef().stkFieldRef();

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
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

    // shifted ip's for field
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
        ws_yScale.resize(nodesPerElement);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_yScale = &ws_yScale[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_dndx = &ws_dndx[0];

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

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = iNode * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }

                // gather 1-dim fields
                p_yScale[iNode] =
                    *stk::mesh::field_data(yScaleSTKFieldRef, node);
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

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather; scalar
                p_yScale[ic] =
                    *stk::mesh::field_data(nodalSideYScaleSTKFieldRef, node);
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

                //================================
                // Diffusion: only diffusion
                //================================

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

                        const scalar yScale = p_yScale[ic];

                        // -dphi/dxj*Aj
                        scalar lhsfac = -dndxj * axj;
                        p_lhs[rowR + ic] += lhsfac * p_bcMultiplier[ic];
                        p_rhs[indexR] -= lhsfac * yScale;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
