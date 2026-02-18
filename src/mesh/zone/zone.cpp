// File : zone.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "zone.h"
#include "boundary.h"
#include "messager.h"

namespace accel
{

zone::zone(mesh* meshPtr, label index, std::string name)
    : meshPtr_(meshPtr), index_(index), name_(name)
{
}

void zone::setup()
{
}

void zone::initialize()
{
    // initialize boundaries
    for (label iBoundary = 0; iBoundary < nBoundaries(); iBoundary++)
    {
        boundaryRef(iBoundary).initialize();
    }

    // compute stats: zone volume
    computeStats0_();
}

void zone::update()
{
    // update boundaries
    for (label iBoundary = 0; iBoundary < nBoundaries(); iBoundary++)
    {
        boundaryRef(iBoundary).update();
    }

    computeStats_();
}

void zone::computeStats0_()
{
    const auto& metaData = this->meshPtr()->metaDataRef();
    const auto& bulkData = this->meshPtr()->bulkDataRef();

    // reset to 0
    stats_.volume_ = 0.0;

    // extract coordinates field
    STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    // setup for element buckets; union parts and ask for locally owned part
    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK,
                             metaData.locally_owned_part() &
                                 stk::mesh::selectUnion(this->interiorParts()));

    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;

        // extract master element for this bucket's topology
        MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCV->nodesPerElement_;
        const label numScvIp = meSCV->numIntPoints_;

        // define scratch field
        std::vector<scalar> ws_coordinates(nodesPerElement * SPATIAL_DIM);
        std::vector<scalar> ws_scv_volume(numScvIp);

        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            // gather nodal coordinates
            stk::mesh::Entity const* nodeRels =
                elementBucket.begin_nodes(iElement);
            label num_nodes = elementBucket.num_nodes(iElement);

            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];
                scalar* coords = stk::mesh::field_data(*coordinates, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_coordinates[offSet + j] = coords[j];
                }
            }

            // compute sub-control volume determinants
            scalar scv_error = 0.0;
            meSCV->determinant(
                1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

            // sum SCV volumes to get element volume
            for (label ip = 0; ip < numScvIp; ++ip)
            {
                stats_.volume_ += ws_scv_volume[ip];
            }
        }
    }

    if (messager::parallel())
    {
        messager::sumReduce(stats_.volume_);
    }

    // store initial stats
    stats_.volume0_ = stats_.volume_;
}

void zone::computeStats_()
{
    const auto& metaData = this->meshPtr()->metaDataRef();
    const auto& bulkData = this->meshPtr()->bulkDataRef();

    // reset to 0
    stats_.volume_ = 0.0;

    // Now we can make use of exposed area field. Get the area field
    STKScalarField& dualVolumeSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

    // setup for buckets; union parts and ask for locally owned part
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK,
                             metaData.locally_owned_part() &
                                 stk::mesh::selectUnion(this->interiorParts()));

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        const scalar* volb =
            stk::mesh::field_data(dualVolumeSTKFieldRef, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stats_.volume_ += volb[iNode];
        }
    }

    if (messager::parallel())
    {
        messager::sumReduce(stats_.volume_);
    }
}

// Access

mesh* zone::meshPtr()
{
    return meshPtr_;
}

const mesh* zone::meshPtr() const
{
    return meshPtr_;
}

mesh& zone::meshRef()
{
    return *meshPtr_;
}

const mesh& zone::meshRef() const
{
    return *meshPtr_;
}

// Boundaries

label zone::nBoundaries() const
{
    return boundaryVector_.size();
}

boundary* zone::boundaryPtr(label iBoundary)
{
    return boundaryVector_[iBoundary].get();
}

const boundary* zone::boundaryPtr(label iBoundary) const
{
    return boundaryVector_[iBoundary].get();
}

boundary& zone::boundaryRef(label iBoundary)
{
    return *boundaryVector_[iBoundary].get();
}

const boundary& zone::boundaryRef(label iBoundary) const
{
    return *boundaryVector_[iBoundary].get();
}

// domain motion
zoneTransformation& zone::transformationRef()
{
    // lazy instantiation
    if (transformationPtr_ == nullptr)
    {
        transformationPtr_ = std::make_unique<zoneTransformation>(this);
    }
    return *transformationPtr_.get();
}

const zoneTransformation& zone::transformationRef() const
{
    // lazy instantiation
    if (transformationPtr_ == nullptr)
    {
        transformationPtr_ =
            std::make_unique<zoneTransformation>(const_cast<zone*>(this));
    }
    return *transformationPtr_.get();
}

// mesh deformation
zoneDeformation& zone::deformationRef()
{
    // lazy instantiation
    if (deformationPtr_ == nullptr)
    {
        deformationPtr_ = std::make_unique<zoneDeformation>(this);
    }
    return *deformationPtr_.get();
}

const zoneDeformation& zone::deformationRef() const
{
    // lazy instantiation
    if (deformationPtr_ == nullptr)
    {
        deformationPtr_ =
            std::make_unique<zoneDeformation>(const_cast<zone*>(this));
    }
    return *deformationPtr_.get();
}

} // namespace accel
