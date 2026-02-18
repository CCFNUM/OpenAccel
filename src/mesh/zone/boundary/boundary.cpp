// File : boundary.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "boundary.h"
#include "messager.h"

namespace accel
{

boundary::boundary(zone* zonePtr, std::string name, label index)
    : zonePtr_(zonePtr), name_(name), index_(index)
{
}

void boundary::read(const YAML::Node& node)
{
    // store boundary type
    physicalType_ = convertBoundaryPhysicalTypeFromString(
        node["type"].template as<std::string>());

    // store parts
    for (const auto& locationName :
         node["location"].template as<std::vector<std::string>>())
    {
        stk::mesh::Part* part =
            zonePtr_->meshRef().metaDataRef().get_part(locationName);

        if (!part)
        {
            errorMsg("Part " + locationName +
                     " not found in the mesh for boundary " + name_);
        }

        parts_.push_back(part);
    }

    // do some checks
    if (zonePtr_->frameRotating() || zonePtr_->meshTransforming())
    {
        // store frame type: only in case of domain motion
        if (node["frame_type"])
        {
            frameType_ = convertBoundaryRelativeFrameTypeFromString(
                node["frame_type"].template as<std::string>());
        }

        // get motion information class
        const auto& motion = zonePtr_->transformationRef();

        switch (motion.type())
        {
            case meshMotionType::stationary:
                assert(frameType_ == boundaryRelativeFrameType::absolute);
                break;

            case meshMotionType::translating:
                // can be absolute or relative
                break;

            case meshMotionType::rotating:
                {
                    if (physicalType_ == boundaryPhysicalType::wall)
                    {
                        assert(frameType_ ==
                               boundaryRelativeFrameType::relative);
                    }
                }
                break;
        }
    }
}

// Methods

void boundary::setup()
{
}

void boundary::initialize()
{
    // compute stats: area and centroid
    computeStats0_();
}

void boundary::update()
{
    computeStats_();
}

void boundary::computeStats0_()
{
    const auto& metaData = zonePtr_->meshPtr()->metaDataRef();
    const auto& bulkData = zonePtr_->meshPtr()->bulkDataRef();

    // set to 0
    stats_.area_ = 0.0;

    // initialze boundary area and centroid
    // extract coordinates field
    STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    // setup for buckets; union parts and ask for owned part
    stk::mesh::BucketVector const& sideBuckets = bulkData.get_buckets(
        metaData.side_rank(),
        metaData.locally_owned_part() & stk::mesh::selectUnion(this->parts()));

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());

        // extract master element specifics
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // define scratch field
        std::vector<scalar> ws_coordinates(nodesPerSide * SPATIAL_DIM);
        std::vector<scalar> ws_scs_areav(numScsBip * SPATIAL_DIM);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // face node relations for nodal gather
            stk::mesh::Entity const* sideNodeRels =
                sideBucket.begin_nodes(iSide);

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            label num_nodes = sideBucket.num_nodes(iSide);
            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];
                scalar* coords = stk::mesh::field_data(*coordinates, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_coordinates[offSet + j] = coords[j];
                }
            }

            // compute scs integration point areavec
            scalar scs_error = 0.0;
            meFC->determinant(
                1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

            // scarrer to area vector
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                scalar amag = 0.0;
                const label offSet = ip * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    amag += ws_scs_areav[offSet + j] * ws_scs_areav[offSet + j];
                }
                amag = std::sqrt(amag);

                // accumulate
                stats_.area_ += amag;
            }
        }
    }

    // sync
    if (messager::parallel())
    {
        messager::sumReduce(stats_.area_);
    }

    // Calculate centroid
    const scalar largeNumber = 1.0e16;
    scalar minCoord[3] = {largeNumber, largeNumber, largeNumber};
    scalar maxCoord[3] = {-largeNumber, -largeNumber, -largeNumber};

    // model coords are fine in this case
    auto* modelCoords = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        this->zonePtr()->meshPtr()->getCoordinateFieldName());

    // select all nodes
    stk::mesh::Selector s_all_nodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(this->parts());

    // select all locally owned nodes for bounding box
    stk::mesh::BucketVector const& node_buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, s_all_nodes);
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
         ib != node_buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* mCoord = stk::mesh::field_data(*modelCoords, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            minCoord[0] = std::min(minCoord[0], mCoord[k * SPATIAL_DIM + 0]);
            maxCoord[0] = std::max(maxCoord[0], mCoord[k * SPATIAL_DIM + 0]);
            minCoord[1] = std::min(minCoord[1], mCoord[k * SPATIAL_DIM + 1]);
            maxCoord[1] = std::max(maxCoord[1], mCoord[k * SPATIAL_DIM + 1]);
            if (SPATIAL_DIM == 3)
            {
                minCoord[2] =
                    std::min(minCoord[2], mCoord[k * SPATIAL_DIM + 2]);
                maxCoord[2] =
                    std::max(maxCoord[2], mCoord[k * SPATIAL_DIM + 2]);
            }
        }
    }

    // parallel reduction on min/max
    scalar g_minCoord[3] = {};
    scalar g_maxCoord[3] = {};
    stk::ParallelMachine comm = bulkData.parallel();
    stk::all_reduce_min(comm, minCoord, g_minCoord, 3);
    stk::all_reduce_max(comm, maxCoord, g_maxCoord, 3);
    for (label j = 0; j < SPATIAL_DIM; ++j)
        stats_.centroid_[j] = 0.5 * (g_maxCoord[j] + g_minCoord[j]);

    // store initial stats
    stats_.area0_ = stats_.area_;
    stats_.centroid0_ = stats_.centroid_;
}

void boundary::computeStats_()
{
    const auto& metaData = zonePtr_->meshPtr()->metaDataRef();
    const auto& bulkData = zonePtr_->meshPtr()->bulkDataRef();

    // set to 0
    stats_.area_ = 0.0;

    // Now we can make use of exposed area field. Get the area field
    STKScalarField& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);

    // setup for buckets; union parts and ask for locally owned part
    stk::mesh::BucketVector const& sideBuckets = bulkData.get_buckets(
        metaData.side_rank(),
        metaData.locally_owned_part() & stk::mesh::selectUnion(this->parts()));

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());

        // extract master element specifics
        const label numScsIp = meFC->numIntPoints_;

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // face data
            scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);

            // scarrer to area vector
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                const label faceOffSet = ip * SPATIAL_DIM;

                scalar localArea = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    localArea +=
                        areaVec[faceOffSet + i] * areaVec[faceOffSet + i];
                }
                localArea = std::sqrt(localArea);

                stats_.area_ += localArea;
            }
        }
    }

    // sync
    if (messager::parallel())
    {
        messager::sumReduce(stats_.area_);
    }

    // Calculate centroid
    const scalar largeNumber = 1.0e16;
    scalar minCoord[3] = {largeNumber, largeNumber, largeNumber};
    scalar maxCoord[3] = {-largeNumber, -largeNumber, -largeNumber};

    // model coords are fine in this case
    auto* modelCoords = metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                   mesh::coordinates_ID);

    // select all nodes
    stk::mesh::Selector s_all_nodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(this->parts());

    // select all locally owned nodes for bounding box
    stk::mesh::BucketVector const& node_buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, s_all_nodes);
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
         ib != node_buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* mCoord = stk::mesh::field_data(*modelCoords, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            minCoord[0] = std::min(minCoord[0], mCoord[k * SPATIAL_DIM + 0]);
            maxCoord[0] = std::max(maxCoord[0], mCoord[k * SPATIAL_DIM + 0]);
            minCoord[1] = std::min(minCoord[1], mCoord[k * SPATIAL_DIM + 1]);
            maxCoord[1] = std::max(maxCoord[1], mCoord[k * SPATIAL_DIM + 1]);
            if (SPATIAL_DIM == 3)
            {
                minCoord[2] =
                    std::min(minCoord[2], mCoord[k * SPATIAL_DIM + 2]);
                maxCoord[2] =
                    std::max(maxCoord[2], mCoord[k * SPATIAL_DIM + 2]);
            }
        }
    }

    // parallel reduction on min/max
    scalar g_minCoord[3] = {};
    scalar g_maxCoord[3] = {};
    stk::ParallelMachine comm = bulkData.parallel();
    stk::all_reduce_min(comm, minCoord, g_minCoord, 3);
    stk::all_reduce_max(comm, maxCoord, g_maxCoord, 3);
    for (label j = 0; j < SPATIAL_DIM; ++j)
        stats_.centroid_[j] = 0.5 * (g_maxCoord[j] + g_minCoord[j]);
}

} // namespace accel
