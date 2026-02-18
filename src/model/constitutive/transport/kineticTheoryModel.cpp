// File : kineticTheoryModel.cpp
// Created : Tue Jun 17 2025 17:05:11 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "kineticTheoryModel.h"

namespace accel
{

kineticTheoryModel::kineticTheoryModel(realm* realm) : transportModel(realm)
{
}

void kineticTheoryModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    updateDynamicViscosity(domain);
}

void kineticTheoryModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateDynamicViscosity(domain, iPhase);
}

void kineticTheoryModel::initializeThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    updateThermalConductivity(domain);
}

void kineticTheoryModel::initializeThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateThermalConductivity(domain, iPhase);
}

void kineticTheoryModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    errorMsg("Not implemented yet");
}

void kineticTheoryModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    errorMsg("Not implemented yet");
}

void kineticTheoryModel::updateThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* lambdaSTKFieldPtr = lambdaRef().stkFieldPtr();

    scalar c1 =
        domain->materialRef().transportProperties_.thermalConductivity_.c1_;
    scalar c2 =
        domain->materialRef().transportProperties_.thermalConductivity_.c2_;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* lambdab = stk::mesh::field_data(*lambdaSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc lambda
            lambdab[iNode] = c1 + c2 * T;
        }
    }
}

void kineticTheoryModel::updateThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* lambdaSTKFieldPtr = lambdaRef(iPhase).stkFieldPtr();

    scalar c1 = domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
                    .transportProperties_.thermalConductivity_.c1_;
    scalar c2 = domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
                    .transportProperties_.thermalConductivity_.c2_;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* lambdab = stk::mesh::field_data(*lambdaSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc lambda
            lambdab[iNode] = c1 + c2 * T;
        }
    }
}

} /* namespace accel */
