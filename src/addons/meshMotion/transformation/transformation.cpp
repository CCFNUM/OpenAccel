// File       : transformation.cpp
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

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
            scalar t = meshMotionPtr_->meshRef().controlsRef().time;

            // rotation details: extract pointers if required
            const auto& omega = zonePtr->transformationRef().rotation().omega_;
            const auto& axis = zonePtr->transformationRef().rotation().axis_;
            const auto& origin =
                zonePtr->transformationRef().rotation().origin_;

            const scalar* p_axis = axis.data();
            const scalar* p_ori = origin.data();

            // translation details
            const auto& v = zonePtr->transformationRef().translation().v_;

            const scalar* p_v = v.data();

            // local space; current coords and rotated coords; generalized for
            // 2D and 3D. Do also for translation.
            scalar mcX[3] = {0.0, 0.0, 0.0};
            scalar rcX[3] = {0.0, 0.0, 0.0};
            scalar tra[3] = {0.0, 0.0, 0.0};

            // Get coords field
            const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::coordinates_ID);
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
                    // load the current and model coords
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        mcX[i] = orgCoordsb[iNode * SPATIAL_DIM + i];
                    }

                    const scalar cX = mcX[0] - p_ori[0];
                    const scalar cY = mcX[1] - p_ori[1];
                    const scalar cZ = mcX[2] - p_ori[2];

                    const scalar sinOTby2 = sin(omega * t * 0.5);
                    const scalar cosOTby2 = cos(omega * t * 0.5);

                    const scalar q0 = cosOTby2;
                    const scalar q1 = sinOTby2 * p_axis[0];
                    const scalar q2 = sinOTby2 * p_axis[1];
                    const scalar q3 = sinOTby2 * p_axis[2];

                    // rotated model coordinates; converted to displacement; add
                    // back in centroid
                    rcX[0] = (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) * cX +
                             2.0 * (q1 * q2 - q0 * q3) * cY +
                             2.0 * (q0 * q2 + q1 * q3) * cZ - mcX[0] + p_ori[0];
                    rcX[1] = 2.0 * (q1 * q2 + q0 * q3) * cX +
                             (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) * cY +
                             2.0 * (q2 * q3 - q0 * q1) * cZ - mcX[1] + p_ori[1];
                    rcX[2] = 2.0 * (q1 * q3 - q0 * q2) * cX +
                             2.0 * (q0 * q1 + q2 * q3) * cY +
                             (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) * cZ -
                             mcX[2] + p_ori[2];

                    tra[0] = p_v[0] * t;
                    tra[1] = p_v[1] * t;
                    tra[2] = p_v[2] * t;

                    // set displacement
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        Dtb[SPATIAL_DIM * iNode + i] += rcX[i] + tra[i];
                    }
                }
            }
        }
    }
}

} /* namespace accel */
