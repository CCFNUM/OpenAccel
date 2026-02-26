// File       : sutherlandsFormulaModel.cpp
// Created    : Thu Apr 03 2025 17:05:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "sutherlandsFormulaModel.h"

namespace accel
{

sutherlandsFormulaModel::sutherlandsFormulaModel(realm* realm)
    : transportModel(realm)
{
}

void sutherlandsFormulaModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    updateDynamicViscosity(domain);
}

void sutherlandsFormulaModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateDynamicViscosity(domain, iPhase);
}

void sutherlandsFormulaModel::initializeThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    updateThermalConductivity(domain);
}

void sutherlandsFormulaModel::initializeThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateThermalConductivity(domain, iPhase);
}

void sutherlandsFormulaModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();

    scalar T_ref =
        domain->materialRef()
            .transportProperties_.dynamicViscosity_.referenceTemperature_;
    scalar mu_ref =
        domain->materialRef()
            .transportProperties_.dynamicViscosity_.referenceViscosity_;
    scalar S = domain->materialRef()
                   .transportProperties_.dynamicViscosity_.sutherlandsConstant_;
    scalar n = domain->materialRef()
                   .transportProperties_.dynamicViscosity_.temperatureExponent_;

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
        scalar* mub = stk::mesh::field_data(*muSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc rho
            mub[iNode] =
                mu_ref * std::pow(T / T_ref, n) * (T_ref + S) / (T + S);
        }
    }
}

void sutherlandsFormulaModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* muSTKFieldPtr = muRef(iPhase).stkFieldPtr();

    scalar T_ref =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .transportProperties_.dynamicViscosity_.referenceTemperature_;
    scalar mu_ref =
        domain->materialRef()
            .transportProperties_.dynamicViscosity_.referenceViscosity_;
    scalar S = domain->materialRef()
                   .transportProperties_.dynamicViscosity_.sutherlandsConstant_;
    scalar n = domain->materialRef()
                   .transportProperties_.dynamicViscosity_.temperatureExponent_;

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
        scalar* mub = stk::mesh::field_data(*muSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc mu
            mub[iNode] =
                mu_ref * std::pow(T / T_ref, n) * (T_ref + S) / (T + S);
        }
    }
}

void sutherlandsFormulaModel::updateThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* lambdaSTKFieldPtr = lambdaRef().stkFieldPtr();

    scalar T_ref =
        domain->materialRef()
            .transportProperties_.thermalConductivity_.referenceTemperature_;
    scalar lambda_ref = domain->materialRef()
                            .transportProperties_.thermalConductivity_
                            .referenceThermalConductivity_;
    scalar S =
        domain->materialRef()
            .transportProperties_.thermalConductivity_.sutherlandsConstant_;
    scalar n =
        domain->materialRef()
            .transportProperties_.thermalConductivity_.temperatureExponent_;

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
            lambdab[iNode] =
                lambda_ref * std::pow(T / T_ref, n) * (T_ref + S) / (T + S);
        }
    }
}

void sutherlandsFormulaModel::updateThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* lambdaSTKFieldPtr = lambdaRef(iPhase).stkFieldPtr();

    scalar T_ref =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .transportProperties_.thermalConductivity_.referenceTemperature_;
    scalar lambda_ref = domain->materialRef()
                            .transportProperties_.thermalConductivity_
                            .referenceThermalConductivity_;
    scalar S =
        domain->materialRef()
            .transportProperties_.thermalConductivity_.sutherlandsConstant_;
    scalar n =
        domain->materialRef()
            .transportProperties_.thermalConductivity_.temperatureExponent_;

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

            // calc rho
            lambdab[iNode] =
                lambda_ref * std::pow(T / T_ref, n) * (T_ref + S) / (T + S);
        }
    }
}

} /* namespace accel */
