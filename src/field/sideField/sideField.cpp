// File       : sideField.cpp
// Created    : Tue Sep 02 2025 17:52:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "sideField.h"
#include "boundary.h"
#include "controls.h"
#include "nodeSideField.h"

namespace accel
{

template <>
void sideField<scalar, 1>::interpolate(const nodeSideField<scalar, 1>& nsf,
                                       label iZone,
                                       label iBoundary,
                                       bool shifted)
{
    stk::mesh::MetaData& metaData = this->meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshPtr_->bulkDataRef();

    // zero-out side field on the current boundary parts
    this->setToValue(
        {0.0}, this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& nsSTKFieldRef = nsf.stkFieldRef();

    // master element
    std::vector<scalar> ws_face_shape_function;

    // nodal fields to gather
    std::vector<scalar> ws_phi;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

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

        // algorithm related; element/face
        ws_phi.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_phi = &ws_phi[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (shifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            scalar* sfvalue = stk::mesh::field_data(stkFieldRef, side);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerSide);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather phi
                p_phi[ni] = *stk::mesh::field_data(nsSTKFieldRef, node);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    sfvalue[ip] += r * p_phi[ic];
                }
            }
        }
    }
}

template <>
void sideField<scalar, SPATIAL_DIM>::interpolate(
    const nodeSideField<scalar, SPATIAL_DIM>& nsf,
    label iZone,
    label iBoundary,
    bool shifted)
{
    stk::mesh::MetaData& metaData = this->meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshPtr_->bulkDataRef();

    // zero-out side field on the current boundary parts
    this->setToValue(
        std::vector<scalar>(SPATIAL_DIM, 0.0),
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& nsSTKFieldRef = nsf.stkFieldRef();

    // master element
    std::vector<scalar> ws_face_shape_function;

    // nodal fields to gather
    std::vector<scalar> ws_phi;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

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

        // algorithm related; element/face
        ws_phi.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_phi = &ws_phi[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (shifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            scalar* sfvalue = stk::mesh::field_data(stkFieldRef, side);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerSide);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather phi
                scalar* phi = stk::mesh::field_data(nsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_phi[offSet + j] = phi[j];
                }
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        sfvalue[SPATIAL_DIM * ip + j] += r * p_phi[icNdim + j];
                    }
                }
            }
        }
    }
}

#ifdef HAS_INTERFACE
template <>
void sideField<scalar, 1>::interpolate(const nodeSideField<scalar, 1>& nsf,
                                       label iInterface,
                                       bool master,
                                       bool shifted)
{
    stk::mesh::MetaData& metaData = this->meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshPtr_->bulkDataRef();

    // get interface info ptr
    auto* interfaceSideInfoPtr =
        master ? this->meshPtr()->interfaceRef(iInterface).masterInfoPtr()
               : this->meshPtr()->interfaceRef(iInterface).slaveInfoPtr();

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& nsSTKFieldRef = nsf.stkFieldRef();

    // master element
    std::vector<scalar> ws_face_shape_function;

    // nodal fields to gather
    std::vector<scalar> ws_phi;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

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

        // algorithm related; element/face
        ws_phi.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_phi = &ws_phi[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (shifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            scalar* sfvalue = stk::mesh::field_data(stkFieldRef, side);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerSide);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather phi
                p_phi[ni] = *stk::mesh::field_data(nsSTKFieldRef, node);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // zero-out
                sfvalue[ip] = 0;

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    sfvalue[ip] += r * p_phi[ic];
                }
            }
        }
    }
}

template <>
void sideField<scalar, SPATIAL_DIM>::interpolate(
    const nodeSideField<scalar, SPATIAL_DIM>& nsf,
    label iInterface,
    bool master,
    bool shifted)
{
    stk::mesh::MetaData& metaData = this->meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshPtr_->bulkDataRef();

    // get interface info ptr
    auto* interfaceSideInfoPtr =
        master ? this->meshPtr()->interfaceRef(iInterface).masterInfoPtr()
               : this->meshPtr()->interfaceRef(iInterface).slaveInfoPtr();

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& nsSTKFieldRef = nsf.stkFieldRef();

    // master element
    std::vector<scalar> ws_face_shape_function;

    // nodal fields to gather
    std::vector<scalar> ws_phi;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

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

        // algorithm related; element/face
        ws_phi.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_phi = &ws_phi[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (shifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            scalar* sfvalue = stk::mesh::field_data(stkFieldRef, side);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerSide);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather phi
                scalar* phi = stk::mesh::field_data(nsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_phi[offSet + j] = phi[j];
                }
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // zero-out
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    sfvalue[SPATIAL_DIM * ip + j] = 0;
                }

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        sfvalue[SPATIAL_DIM * ip + j] += r * p_phi[icNdim + j];
                    }
                }
            }
        }
    }
}

template <>
void sideField<scalar, 1>::transfer(label iInterface,
                                    bool reverse,
                                    bool shifted)
{
    bool verbose =
        this->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.advancedOptions_.interfaceTransfer_.verbose_ > 0;

    assert(iInterface < this->meshRef().nInterfaces());

    const interface& interf = this->meshRef().interfaceRef(iInterface);

    // check if field is defined on the first interface side
    if (!definedOn(interf.masterInfoPtr()->currentPartVec_))
    {
        errorMsg("Field " + this->name() +
                 " is not defined on one or more master parts of interface " +
                 interf.name());
    }

    // check if field is defined on the first interface side
    if (!this->definedOn(interf.masterInfoPtr()->opposingPartVec_))
    {
        errorMsg("Field " + this->name() +
                 " is not defined on one or more slave parts of interface " +
                 interf.name());
    }

    // get interface side that is sitting in the solid domain
    const interfaceSideInfo* interfaceSideInfoPtr =
        reverse ? interf.masterInfoPtr() : interf.slaveInfoPtr();

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // get fields
    auto& sideSTKFieldRef = this->stkFieldRef();

    if (interf.isConformal())
    {
        if (verbose && messager::master())
        {
            const interfaceSideInfo* opposingSide =
                reverse ? interf.slaveInfoPtr() : interf.masterInfoPtr();

            std::cout << "\n[Transfer: " << this->name() << "]"
                      << "\n\tConformal: " << opposingSide->name() << " -> "
                      << interfaceSideInfoPtr->name() << std::endl;
        }

        // extract vector of dgInfo
        const std::vector<std::vector<dgInfo*>>& dgInfoVec =
            interfaceSideInfoPtr->dgInfoVec_;

        for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size());
             iSide++)
        {
            const std::vector<dgInfo*>& faceDgInfoVec = dgInfoVec[iSide];

            for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
            {
                dgInfo* dgInfo = faceDgInfoVec[k];

                if (dgInfo->gaussPointExposed_)
                    continue;

                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;

                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                const label opposingGaussPointId =
                    dgInfo->opposingGaussPointId_;

                scalar* c_field =
                    stk::mesh::field_data(sideSTKFieldRef, currentFace);
                const scalar* o_field =
                    stk::mesh::field_data(sideSTKFieldRef, opposingFace);

                c_field[currentGaussPointId] = o_field[opposingGaussPointId];
            }
        }
    }
    else
    {
        if (verbose && messager::master())
        {
            const interfaceSideInfo* opposingSide =
                reverse ? interf.slaveInfoPtr() : interf.masterInfoPtr();

            std::cout << "\n[Transfer: " << this->name() << "]"
                      << "\n\tNon-conformal: " << opposingSide->name() << " -> "
                      << interfaceSideInfoPtr->name() << std::endl;
        }

        std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

        // Workspace for ip-to-node conversion.
        // M[ip][node] = N_node(xi_ip) is the shape-function matrix at the
        // integration points.  Inverting it recovers pseudo-nodal values
        // from ip values so that interpolatePoint (node-based) can be reused.
        std::vector<scalar> ws_M_inv;
        std::vector<scalar> ws_pseudo_nodal;
        MasterElement* prevMeFCOpposing = nullptr;

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

                // extract current/opposing face
                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;
                MasterElement* meFCOpposing = dgInfo->meFCOpposing_;

                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

                const label npe = meFCOpposing->nodesPerElement_;

                // Recompute M^{-1} when the opposing topology changes
                if (meFCOpposing != prevMeFCOpposing)
                {
                    prevMeFCOpposing = meFCOpposing;

                    std::vector<scalar> ws_shape_fcn(npe * npe);
                    if (shifted)
                        meFCOpposing->shifted_shape_fcn(&ws_shape_fcn[0]);
                    else
                        meFCOpposing->shape_fcn(&ws_shape_fcn[0]);

                    // Gauss-Jordan elimination with partial pivoting
                    std::vector<scalar> M(ws_shape_fcn);
                    ws_M_inv.assign(npe * npe, 0.0);
                    for (label i = 0; i < npe; ++i)
                        ws_M_inv[i * npe + i] = 1.0;

                    for (label col = 0; col < npe; ++col)
                    {
                        label pivotRow = col;
                        scalar pivotVal = std::abs(M[col * npe + col]);
                        for (label row = col + 1; row < npe; ++row)
                        {
                            const scalar val = std::abs(M[row * npe + col]);
                            if (val > pivotVal)
                            {
                                pivotVal = val;
                                pivotRow = row;
                            }
                        }
                        if (pivotRow != col)
                        {
                            for (label j = 0; j < npe; ++j)
                            {
                                std::swap(M[col * npe + j],
                                          M[pivotRow * npe + j]);
                                std::swap(ws_M_inv[col * npe + j],
                                          ws_M_inv[pivotRow * npe + j]);
                            }
                        }

                        const scalar scale = 1.0 / M[col * npe + col];
                        for (label j = 0; j < npe; ++j)
                        {
                            M[col * npe + j] *= scale;
                            ws_M_inv[col * npe + j] *= scale;
                        }

                        for (label row = 0; row < npe; ++row)
                        {
                            if (row == col)
                                continue;
                            const scalar factor = M[row * npe + col];
                            for (label j = 0; j < npe; ++j)
                            {
                                M[row * npe + j] -= factor * M[col * npe + j];
                                ws_M_inv[row * npe + j] -=
                                    factor * ws_M_inv[col * npe + j];
                            }
                        }
                    }

                    ws_pseudo_nodal.resize(npe);
                }

                scalar* c_field =
                    stk::mesh::field_data(sideSTKFieldRef, currentFace);
                const scalar* o_field =
                    stk::mesh::field_data(sideSTKFieldRef, opposingFace);

                // Convert ip values to pseudo-nodal values:
                // c_node = sum_ip( M_inv[node][ip] * f_ip )
                for (label node = 0; node < npe; ++node)
                {
                    ws_pseudo_nodal[node] = 0.0;
                    for (label ip = 0; ip < npe; ++ip)
                    {
                        ws_pseudo_nodal[node] +=
                            ws_M_inv[node * npe + ip] * o_field[ip];
                    }
                }

                // Interpolate from pseudo-nodal values to the projection point
                scalar interpValue = 0.0;
                meFCOpposing->interpolatePoint(1,
                                               &opposingIsoParCoords[0],
                                               &ws_pseudo_nodal[0],
                                               &interpValue);

                c_field[currentGaussPointId] = interpValue;
            }
        }
    } // else (non-conformal)

    if (verbose)
    {
        scalar sumCurrent = 0;
        scalar sumOpposing = 0;

        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), mesh::exposed_area_vector_ID);

        // current total sum
        {
            // define some common selectors
            stk::mesh::Selector selOwnedSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nBoundaryFaces =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nBoundaryFaces;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                    const scalar* field =
                        stk::mesh::field_data(sideSTKFieldRef, side);

                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        // zero out vector quantities; squeeze in
                        // aMag
                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        sumCurrent += aMag * field[ip];
                    }
                }
            }
        }

        // opposing total sum
        {
            // define some common selectors
            stk::mesh::Selector selOwnedSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(interfaceSideInfoPtr->opposingPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nBoundaryFaces =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nBoundaryFaces;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                    const scalar* field =
                        stk::mesh::field_data(sideSTKFieldRef, side);

                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        // zero out vector quantities; squeeze in
                        // aMag
                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        sumOpposing += aMag * field[ip];
                    }
                }
            }
        }

        if (messager::parallel())
        {
            messager::sumReduce(sumCurrent);
            messager::sumReduce(sumOpposing);
        }

        if (messager::master())
        {
            const scalar diff = sumCurrent - sumOpposing;
            std::cout << "\n\tTransfer integral check (" << this->name()
                      << "):" << "\n\t  " << std::setw(18) << std::left
                      << "Current" << std::setw(18) << "Opposing"
                      << "Difference" << "\n\t  " << std::scientific
                      << std::setprecision(8) << std::setw(18) << sumCurrent
                      << std::setw(18) << sumOpposing << diff << std::endl;
        }
    }
}

template <>
void sideField<scalar, SPATIAL_DIM>::transfer(label iInterface,
                                              bool reverse,
                                              bool shifted)
{
    bool verbose =
        this->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.advancedOptions_.interfaceTransfer_.verbose_ > 0;

    assert(iInterface < this->meshRef().nInterfaces());

    const interface& interf = this->meshRef().interfaceRef(iInterface);

    // check if field is defined on the first interface side
    if (!definedOn(interf.masterInfoPtr()->currentPartVec_))
    {
        errorMsg("Field " + this->name() +
                 " is not defined on one or more master parts of interface " +
                 interf.name());
    }

    // check if field is defined on the first interface side
    if (!this->definedOn(interf.masterInfoPtr()->opposingPartVec_))
    {
        errorMsg("Field " + this->name() +
                 " is not defined on one or more slave parts of interface " +
                 interf.name());
    }

    // get interface side that is sitting in the solid domain
    const interfaceSideInfo* interfaceSideInfoPtr =
        reverse ? interf.masterInfoPtr() : interf.slaveInfoPtr();

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    auto& sideSTKFieldRef = this->stkFieldRef();

    if (interf.isConformal())
    {
        if (verbose && messager::master())
        {
            const interfaceSideInfo* opposingSide =
                reverse ? interf.slaveInfoPtr() : interf.masterInfoPtr();

            std::cout << "\n[Transfer: " << this->name() << "]"
                      << "\n\tConformal: " << opposingSide->name() << " -> "
                      << interfaceSideInfoPtr->name() << std::endl;
        }

        // extract vector of dgInfo
        const std::vector<std::vector<dgInfo*>>& dgInfoVec =
            interfaceSideInfoPtr->dgInfoVec_;

        for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size());
             iSide++)
        {
            const std::vector<dgInfo*>& faceDgInfoVec = dgInfoVec[iSide];

            for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
            {
                dgInfo* dgInfo = faceDgInfoVec[k];

                if (dgInfo->gaussPointExposed_)
                    continue;

                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;

                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                const label opposingGaussPointId =
                    dgInfo->opposingGaussPointId_;

                scalar* c_field =
                    stk::mesh::field_data(sideSTKFieldRef, currentFace);
                const scalar* o_field =
                    stk::mesh::field_data(sideSTKFieldRef, opposingFace);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    c_field[currentGaussPointId * SPATIAL_DIM + j] =
                        o_field[opposingGaussPointId * SPATIAL_DIM + j];
                }
            }
        }
    }
    else
    {
        if (verbose && messager::master())
        {
            const interfaceSideInfo* opposingSide =
                reverse ? interf.slaveInfoPtr() : interf.masterInfoPtr();

            std::cout << "\n[Transfer: " << this->name() << "]"
                      << "\n\tNon-conformal: " << opposingSide->name() << " -> "
                      << interfaceSideInfoPtr->name() << std::endl;
        }

        std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

        // Workspace for ip-to-node conversion and transposed field.
        // M[ip][node] = N_node(xi_ip) is the shape-function matrix at the
        // integration points.  Inverting it recovers pseudo-nodal values
        // from ip values so that interpolatePoint (node-based) can be reused.
        std::vector<scalar> ws_M_inv;
        std::vector<scalar> ws_pseudo_nodal;
        MasterElement* prevMeFCOpposing = nullptr;

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

                // extract current/opposing face
                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;
                MasterElement* meFCOpposing = dgInfo->meFCOpposing_;

                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

                const label npe = meFCOpposing->nodesPerElement_;

                // Recompute M^{-1} when the opposing topology changes
                if (meFCOpposing != prevMeFCOpposing)
                {
                    prevMeFCOpposing = meFCOpposing;

                    std::vector<scalar> ws_shape_fcn(npe * npe);
                    if (shifted)
                        meFCOpposing->shifted_shape_fcn(&ws_shape_fcn[0]);
                    else
                        meFCOpposing->shape_fcn(&ws_shape_fcn[0]);

                    // Gauss-Jordan elimination with partial pivoting
                    std::vector<scalar> M(ws_shape_fcn);
                    ws_M_inv.assign(npe * npe, 0.0);
                    for (label i = 0; i < npe; ++i)
                        ws_M_inv[i * npe + i] = 1.0;

                    for (label col = 0; col < npe; ++col)
                    {
                        label pivotRow = col;
                        scalar pivotVal = std::abs(M[col * npe + col]);
                        for (label row = col + 1; row < npe; ++row)
                        {
                            const scalar val = std::abs(M[row * npe + col]);
                            if (val > pivotVal)
                            {
                                pivotVal = val;
                                pivotRow = row;
                            }
                        }
                        if (pivotRow != col)
                        {
                            for (label j = 0; j < npe; ++j)
                            {
                                std::swap(M[col * npe + j],
                                          M[pivotRow * npe + j]);
                                std::swap(ws_M_inv[col * npe + j],
                                          ws_M_inv[pivotRow * npe + j]);
                            }
                        }

                        const scalar scale = 1.0 / M[col * npe + col];
                        for (label j = 0; j < npe; ++j)
                        {
                            M[col * npe + j] *= scale;
                            ws_M_inv[col * npe + j] *= scale;
                        }

                        for (label row = 0; row < npe; ++row)
                        {
                            if (row == col)
                                continue;
                            const scalar factor = M[row * npe + col];
                            for (label j = 0; j < npe; ++j)
                            {
                                M[row * npe + j] -= factor * M[col * npe + j];
                                ws_M_inv[row * npe + j] -=
                                    factor * ws_M_inv[col * npe + j];
                            }
                        }
                    }

                    ws_pseudo_nodal.resize(npe * SPATIAL_DIM);
                }

                scalar* c_field =
                    stk::mesh::field_data(sideSTKFieldRef, currentFace);
                const scalar* o_field =
                    stk::mesh::field_data(sideSTKFieldRef, opposingFace);

                // Convert ip values to pseudo-nodal values in
                // component-major layout for interpolatePoint:
                // c[comp * npe + node] = sum_ip( M_inv[node][ip] *
                //                                f_ip[ip * NDIM + comp] )
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    for (label node = 0; node < npe; ++node)
                    {
                        ws_pseudo_nodal[j * npe + node] = 0.0;
                        for (label ip = 0; ip < npe; ++ip)
                        {
                            ws_pseudo_nodal[j * npe + node] +=
                                ws_M_inv[node * npe + ip] *
                                o_field[ip * SPATIAL_DIM + j];
                        }
                    }
                }

                // Interpolate from pseudo-nodal values to the projection point
                scalar interpValue[SPATIAL_DIM] = {0};
                meFCOpposing->interpolatePoint(SPATIAL_DIM,
                                               &opposingIsoParCoords[0],
                                               &ws_pseudo_nodal[0],
                                               interpValue);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    c_field[currentGaussPointId * SPATIAL_DIM + j] =
                        interpValue[j];
                }
            }
        }
    } // else (non-conformal)

    if (verbose)
    {
        std::vector<scalar> sumCurrent(SPATIAL_DIM, 0);
        std::vector<scalar> sumOpposing(SPATIAL_DIM, 0);

        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), mesh::exposed_area_vector_ID);

        // current total sum
        {
            // define some common selectors
            stk::mesh::Selector selOwnedSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nBoundaryFaces =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nBoundaryFaces;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                    const scalar* field =
                        stk::mesh::field_data(sideSTKFieldRef, side);

                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        // zero out vector quantities; squeeze in
                        // aMag
                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            sumCurrent[j] += aMag * field[ip * SPATIAL_DIM + j];
                        }
                    }
                }
            }
        }

        // opposing total sum
        {
            // define some common selectors
            stk::mesh::Selector selOwnedSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(interfaceSideInfoPtr->opposingPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nBoundaryFaces =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nBoundaryFaces;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                    const scalar* field =
                        stk::mesh::field_data(sideSTKFieldRef, side);

                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        // zero out vector quantities; squeeze in
                        // aMag
                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            sumOpposing[j] +=
                                aMag * field[ip * SPATIAL_DIM + j];
                        }
                    }
                }
            }
        }

        if (messager::parallel())
        {
            messager::sumReduce(sumCurrent);
            messager::sumReduce(sumOpposing);
        }

        if (messager::master())
        {
            const char* compLabel[] = {"x", "y", "z"};
            std::cout << "\n\tTransfer integral check:" << "\n\t  "
                      << std::setw(4) << " " << std::setw(18) << std::left
                      << "Current" << std::setw(18) << "Opposing"
                      << "Difference" << std::endl;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar diff = sumCurrent[j] - sumOpposing[j];
                std::cout << "\t  " << std::setw(4) << std::left << compLabel[j]
                          << std::scientific << std::setprecision(8)
                          << std::setw(18) << sumCurrent[j] << std::setw(18)
                          << sumOpposing[j] << diff << std::endl;
            }
        }
    }
}
#endif /* HAS_INTERFACE */

} // namespace accel
