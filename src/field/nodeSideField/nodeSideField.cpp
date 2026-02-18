// File : nodeSideField.cpp
// Created : Tue Sep 02 2025 17:41:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "nodeSideField.h"
#include "boundary.h"
#include "sideField.h"

namespace accel
{

template <>
void nodeSideField<scalar, 1>::interpolate(const sideField<scalar, 1>& sf,
                                           label iZone,
                                           label iBoundary)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // zero-out side field on the current boundary parts
    this->setToValue(
        {0.0}, this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // get assembled area field
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(iZone));

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& sfSTKFieldRef = sf.stkFieldRef();

    // flag whether to calculate assembled area field. This is a must in two
    // cases: 1) the first time the code visits this function 2) in case the
    // zone is deforming
    bool areaCalc = false;

    // create and fill the moving area field
    if (localAssembledAreaFieldPtr_ == nullptr)
    {
        // create field
        localAssembledAreaFieldPtr_ = std::make_unique<simpleScalarField>(
            this->meshPtr(),
            stk::topology::NODE_RANK,
            this->name() + "_assembled_area");

        // set flag to true
        areaCalc = true;
    }

    // ensure the local area field is created
    assert(localAssembledAreaFieldPtr_);

    // get stk pointer
    auto* localAssembledAreaSTKFieldPtr =
        localAssembledAreaFieldPtr_->stkFieldPtr();

    // define field on boundary parts if not yet
    for (const stk::mesh::Part* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        if (!(localAssembledAreaSTKFieldPtr->defined_on(*part)))
        {
            // put on all parts which the current node-side field is defined on
            scalar initialValue = 0.0;
            stk::mesh::put_field_on_mesh(
                *localAssembledAreaSTKFieldPtr, *part, 1, &initialValue);

            // set flag to true
            areaCalc = true;
        }
    }

    if (this->meshPtr()->zonePtr(iZone)->meshDeforming())
    {
        areaCalc = true;
    }

    // calculate assembled area for the current patch
    if (areaCalc)
    {
        // zero first
        {
            stk::mesh::Selector selAllNodes =
                metaData.universal_part() &
                stk::mesh::selectUnion(this->meshPtr()
                                           ->zonePtr(iZone)
                                           ->boundaryRef(iBoundary)
                                           .parts());
            stk::mesh::BucketVector const& sideNodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

            for (stk::mesh::BucketVector::const_iterator ib =
                     sideNodeBuckets.begin();
                 ib != sideNodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideNodeBucket = **ib;
                const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                    sideNodeBucket.size();

                scalar* valueb = stk::mesh::field_data(
                    *localAssembledAreaSTKFieldPtr, sideNodeBucket);

                for (stk::mesh::Bucket::size_type iSideNode = 0;
                     iSideNode < nSideNodesPerBucket;
                     ++iSideNode)
                {
                    valueb[iSideNode] = 0.0;
                }
            }
        }

        // fill in with assembled area
        {
            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    accel::MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                // mapping from ip to nodes for this ordinal; face perspective
                // (use with face_node_relations)
                const label* faceIpNodeMap = meFC->ipNodeMap();

                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nSidesPerBucket =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                    //======================================
                    // gather nodal data off of face
                    //======================================
                    stk::mesh::Entity const* sideNodeRels =
                        bulkData.begin_nodes(side);

                    // loop over boundary ips
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        const label offSetAveraVec = ip * SPATIAL_DIM;

                        const label localFaceNode = faceIpNodeMap[ip];

                        // left and right nodes; right is on the face; left is
                        // the opposing node
                        stk::mesh::Entity node = sideNodeRels[localFaceNode];

                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[offSetAveraVec + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        scalar* ma = stk::mesh::field_data(
                            *localAssembledAreaSTKFieldPtr, node);
                        (*ma) += aMag;
                    }
                }
            }
        }

        // communicate if parallel
        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData,
                                              {localAssembledAreaSTKFieldPtr});
        }
    }

    // accumulate
    {
        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());

            // mapping from ip to nodes for this ordinal; face perspective
            // (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            const label numScsBip = meFC->numIntPoints_;

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);

                const scalar* svalue =
                    stk::mesh::field_data(sfSTKFieldRef, side);

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on the face; left is
                    // the opposing node
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    scalar* value = stk::mesh::field_data(stkFieldRef, node);
                    const scalar aa = *stk::mesh::field_data(
                        *localAssembledAreaSTKFieldPtr, node);

                    (*value) += aMag * svalue[ip] / aa;
                }
            }
        }
    }

    // communicate if parallel
    if (messager::parallel())
    {
        this->synchronizeGhostedEntities(iZone, iBoundary);
    }
}

template <>
void nodeSideField<scalar, SPATIAL_DIM>::interpolate(
    const sideField<scalar, SPATIAL_DIM>& sf,
    label iZone,
    label iBoundary)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // zero-out side field on the current boundary parts
    this->setToValue(
        std::vector<scalar>(SPATIAL_DIM, 0.0),
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    // get assembled area field
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(iZone));

    // get field
    auto& stkFieldRef = this->stkFieldRef();
    const auto& sfSTKFieldRef = sf.stkFieldRef();

    // flag whether to calculate assembled area field. This is a must in two
    // cases: 1) the first time the code visits this function 2) in case the
    // zone is deforming
    bool areaCalc = false;

    // create and fill the moving area field
    if (localAssembledAreaFieldPtr_ == nullptr)
    {
        // create field
        localAssembledAreaFieldPtr_ = std::make_unique<simpleScalarField>(
            this->meshPtr(),
            stk::topology::NODE_RANK,
            this->name() + "_assembled_area");

        // set flag to true
        areaCalc = true;
    }

    // ensure the local area field is created
    assert(localAssembledAreaFieldPtr_);

    // get stk pointer
    auto* localAssembledAreaSTKFieldPtr =
        localAssembledAreaFieldPtr_->stkFieldPtr();

    // define field on boundary parts if not yet
    for (const stk::mesh::Part* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        if (!(localAssembledAreaSTKFieldPtr->defined_on(*part)))
        {
            // put on all parts which the current node-side field is defined on
            scalar initialValue = 0.0;
            stk::mesh::put_field_on_mesh(
                *localAssembledAreaSTKFieldPtr, *part, 1, &initialValue);

            // set flag to true
            areaCalc = true;
        }
    }

    if (this->meshPtr()->zonePtr(iZone)->meshDeforming())
    {
        areaCalc = true;
    }

    // calculate assembled area for the current patch
    if (areaCalc)
    {
        // zero first
        {
            stk::mesh::Selector selAllNodes =
                metaData.universal_part() &
                stk::mesh::selectUnion(this->meshPtr()
                                           ->zonePtr(iZone)
                                           ->boundaryRef(iBoundary)
                                           .parts());
            stk::mesh::BucketVector const& sideNodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

            for (stk::mesh::BucketVector::const_iterator ib =
                     sideNodeBuckets.begin();
                 ib != sideNodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideNodeBucket = **ib;
                const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                    sideNodeBucket.size();

                scalar* valueb = stk::mesh::field_data(
                    *localAssembledAreaSTKFieldPtr, sideNodeBucket);

                for (stk::mesh::Bucket::size_type iSideNode = 0;
                     iSideNode < nSideNodesPerBucket;
                     ++iSideNode)
                {
                    valueb[iSideNode] = 0.0;
                }
            }
        }

        // fill in with assembled area
        {
            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                // face master element
                MasterElement* meFC =
                    accel::MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                // mapping from ip to nodes for this ordinal; face perspective
                // (use with face_node_relations)
                const label* faceIpNodeMap = meFC->ipNodeMap();

                const label numScsBip = meFC->numIntPoints_;

                const stk::mesh::Bucket::size_type nSidesPerBucket =
                    sideBucket.size();

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // pointer to face data
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                    //======================================
                    // gather nodal data off of face
                    //======================================
                    stk::mesh::Entity const* sideNodeRels =
                        bulkData.begin_nodes(side);

                    // loop over boundary ips
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        const label offSetAveraVec = ip * SPATIAL_DIM;

                        const label localFaceNode = faceIpNodeMap[ip];

                        // left and right nodes; right is on the face; left is
                        // the opposing node
                        stk::mesh::Entity node = sideNodeRels[localFaceNode];

                        scalar aMag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[offSetAveraVec + j];
                            aMag += axj * axj;
                        }
                        aMag = std::sqrt(aMag);

                        scalar* ma = stk::mesh::field_data(
                            *localAssembledAreaSTKFieldPtr, node);
                        (*ma) += aMag;
                    }
                }
            }
        }

        // communicate if parallel
        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData,
                                              {localAssembledAreaSTKFieldPtr});
        }
    }

    // accumulate
    {
        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());

            // mapping from ip to nodes for this ordinal; face perspective
            // (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            const label numScsBip = meFC->numIntPoints_;

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);

                const scalar* svalue =
                    stk::mesh::field_data(sfSTKFieldRef, side);

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on the face; left is
                    // the opposing node
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    scalar* value = stk::mesh::field_data(stkFieldRef, node);
                    const scalar aa = *stk::mesh::field_data(
                        *localAssembledAreaSTKFieldPtr, node);

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        value[i] += aMag * svalue[SPATIAL_DIM * ip + i] / aa;
                    }
                }
            }
        }
    }

    // communicate if parallel
    if (messager::parallel())
    {
        this->synchronizeGhostedEntities(iZone, iBoundary);
    }
}

} // namespace accel
