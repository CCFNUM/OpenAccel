// File : sideField.cpp
// Created : Tue Sep 02 2025 17:52:24 (+0100)
// Author : Mhamad Mahdi Alloush
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
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
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
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
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

} // namespace accel
