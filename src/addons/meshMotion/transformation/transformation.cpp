// File       : transformation.cpp
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "transformation.h"
#include "meshMotion.h"

namespace accel
{

transformation::transformation(meshMotion* meshMotionPtr)
    : meshMotionPtr_(meshMotionPtr)
{
}

void transformation::setup()
{
}

void transformation::initialize()
{
}

void transformation::update()
{
    // calculate displacement due to zone motion
    for (label iZone = 0; iZone < meshMotionPtr_->meshRef().nZones(); iZone++)
    {
        const zone* zonePtr = meshMotionPtr_->meshRef().zonePtr(iZone);

        if (zonePtr->meshTransforming())
        {
#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "updating total displacement field due to "
                             "transformation in zone: "
                          << zonePtr->name() << std::endl;
            }
#endif

            auto& mesh = meshMotionPtr_->meshRef();
            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
            stk::mesh::MetaData& metaData = mesh.metaDataRef();

            // current time
            const scalar t = meshMotionPtr_->meshRef().controlsRef().time;

            // rotation and translation details
            const scalar& omega =
                zonePtr->transformationRef().rotation().omega_;
            const rvectorD& origin =
                zonePtr->transformationRef().rotation().origin_;
// axis is not used in a 2D case
#if SPATIAL_DIM == 3
            const rvectorD& axis =
                zonePtr->transformationRef().rotation().axis_;
#else
            const rvectorD& axis = rvectorD::Zero();
#endif
            const rvectorD& v = zonePtr->transformationRef().translation().v_;

            // --- build the rigid-body transformation once per zone ---
            const rtensorD R = utils::getRotationMatrix(omega * t, axis);
            // translation displacement: v * t
            const rvectorD tra = v * t;

            // Get original/model coords field
            const auto& orgCoordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::original_coordinates_ID);

            // Get fields
            const auto& DtSTKFieldRef = meshMotionPtr_->DtRef().stkFieldRef();

            // get interior parts of the zone
            const stk::mesh::PartVector& intPartVec = zonePtr->interiorParts();

            // get stationary parts of the zone
            const stk::mesh::PartVector& statPartVec =
                zonePtr->stationaryParts();

            // define some common selectors; select all nodes
            stk::mesh::Selector selUniversalNodes =
                metaData.universal_part() & stk::mesh::selectUnion(intPartVec) &
                !stk::mesh::selectUnion(statPartVec);

            stk::mesh::BucketVector const& nodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK, selUniversalNodes);

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                scalar* Dtb = stk::mesh::field_data(DtSTKFieldRef, nodeBucket);
                const scalar* orgCoordsb =
                    stk::mesh::field_data(orgCoordsSTKFieldRef, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    // map flat coordinate arrays into Eigen vectors
                    // (zero-copy views over the STK field data)
                    const rvectorDViewC mcX(orgCoordsb + iNode * SPATIAL_DIM);
                    rvectorDView dt(Dtb + iNode * SPATIAL_DIM);

                    // total displacement:
                    //   rot_disp = R * (mcX - origin) + origin - mcX
                    //   tra      = v * t
                    dt += R * (mcX - origin) + origin - mcX + tra;
                }
            }
        }
    }
}

} /* namespace accel */
