// File : idealGasModel.cpp
// Created : Thu Apr 03 2025 17:05:11 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "idealGasModel.h"
#include "polynomial.h"

namespace accel
{

idealGasModel::idealGasModel(realm* realm) : thermoModel(realm)
{
}

void idealGasModel::initializeDensity(const std::shared_ptr<domain> domain)
{
    updateDensity(domain);
}

void idealGasModel::initializeDensity(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    updateDensity(domain, iPhase);
}

void idealGasModel::initializeCompressibility(
    const std::shared_ptr<domain> domain)
{
    updateCompressibility(domain);
}

void idealGasModel::initializeCompressibility(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateCompressibility(domain, iPhase);
}

void idealGasModel::initializeSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    updateSpecificHeatCapacity(domain);
}

void idealGasModel::initializeSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateSpecificHeatCapacity(domain, iPhase);
}

void idealGasModel::initializeSpecificEnthalpy(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();

    const auto& coeffs =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.coeffs_;

    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

    // Precompute enthalpy polynomial coefficients: h_i = Rs * a_i / (i + 1)
    std::vector<scalar> hCoeffs(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i)
    {
        hCoeffs[i] = Rs * coeffs[i] / static_cast<scalar>(i + 1);
    }

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Select owned nodes in the domain
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    const stk::mesh::BucketVector& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const stk::mesh::Bucket* b : nodeBuckets)
    {
        const stk::mesh::Bucket& nodeBucket = *b;
        const auto nNodesPerBucket = nodeBucket.size();

        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);

        for (size_t iNode = 0; iNode < nNodesPerBucket; ++iNode)
        {
            scalar T = Tb[iNode];
            scalar Ti = T;
            scalar h = 0.0;

            for (size_t i = 0; i < hCoeffs.size(); ++i)
            {
                h += hCoeffs[i] * Ti;
                Ti *= T; // Ti = T^(i+2) in next iteration
            }

            hb[iNode] = h;
        }
    }
}

void idealGasModel::initializeSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();
    const STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
    STKScalarField* h0STKFieldPtr = h0Ref().stkFieldPtr();

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
        const scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
        const scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
        scalar* h0b = stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // caluclate UmagSqr
            scalar UmagSqr = 0.0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                UmagSqr +=
                    Ub[SPATIAL_DIM * iNode + i] * Ub[SPATIAL_DIM * iNode + i];
            }

            scalar h = hb[iNode];

            h0b[iNode] = h + 0.5 * UmagSqr;
        }
    }
}

void idealGasModel::initializeTemperature(const std::shared_ptr<domain> domain)
{
    errorMsg("Must not reach here");
}

void idealGasModel::updateDensity(const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr = pRef().stkFieldPtr();
    STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();

    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        const scalar* pb = stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar p = pb[iNode];
            scalar T = Tb[iNode];

            // calc rho
            rhob[iNode] = p / (Rs * T);
        }
    }

    // add contribution from hydrostatic pressure in case of buoyancy
    if (domain->buoyancy_.option_ == buoyancyOption::buoyant)
    {
        switch (domain->buoyancy_.model_)
        {
            case buoyancyModel::full:
                {
                    // get gravity
                    const auto& gravity = domain->buoyancy_.gravity_;

                    // get geometric fields
                    STKScalarField* coordsSTKFieldPtr =
                        metaData.get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(domain));

                    const auto& buoyancyReferenceDensity =
                        domain->buoyancy_.referenceDensity_;
                    const auto& referenceLocation =
                        domain->buoyancy_.referenceLocation_;

                    for (stk::mesh::BucketVector::const_iterator ib =
                             nodeBuckets.begin();
                         ib != nodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& nodeBucket = **ib;

                        const stk::mesh::Bucket::size_type nNodesPerBucket =
                            nodeBucket.size();

                        // field chunks in bucket
                        const scalar* Tb =
                            stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
                        scalar* rhob =
                            stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);

                        for (stk::mesh::Bucket::size_type iNode = 0;
                             iNode < nNodesPerBucket;
                             ++iNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                *coordsSTKFieldPtr, nodeBucket, iNode);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                rhob[iNode] +=
                                    (buoyancyReferenceDensity * gravity[i] *
                                     (coords[i] - referenceLocation[i])) /
                                    (Rs * Tb[iNode]);
                            }
                        }
                    }
                }
                break;

            case buoyancyModel::boussinesq:
                {
                    errorMsg("Must not reach here");
                }
                break;
        }
    }
}

void idealGasModel::updateDensity(const std::shared_ptr<domain> domain,
                                  label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr = pRef().stkFieldPtr();
    STKScalarField* rhoSTKFieldPtr = rhoRef(iPhase).stkFieldPtr();

    scalar molarMass =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        const scalar* pb = stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar p = pb[iNode];
            scalar T = Tb[iNode];

            // calc rho
            rhob[iNode] = p / (Rs * T);
        }
    }

    // add contribution from hydrostatic pressure in case of buoyancy
    if (domain->buoyancy_.option_ == buoyancyOption::buoyant)
    {
        switch (domain->buoyancy_.model_)
        {
            case buoyancyModel::full:
                {
                    // get gravity
                    const auto& gravity = domain->buoyancy_.gravity_;

                    // get geometric fields
                    STKScalarField* coordsSTKFieldPtr =
                        metaData.get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(domain));

                    const auto& buoyancyReferenceDensity =
                        domain->buoyancy_.referenceDensity_;
                    const auto& referenceLocation =
                        domain->buoyancy_.referenceLocation_;

                    for (stk::mesh::BucketVector::const_iterator ib =
                             nodeBuckets.begin();
                         ib != nodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& nodeBucket = **ib;

                        const stk::mesh::Bucket::size_type nNodesPerBucket =
                            nodeBucket.size();

                        // field chunks in bucket
                        const scalar* Tb =
                            stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
                        scalar* rhob =
                            stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);

                        for (stk::mesh::Bucket::size_type iNode = 0;
                             iNode < nNodesPerBucket;
                             ++iNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                *coordsSTKFieldPtr, nodeBucket, iNode);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                rhob[iNode] +=
                                    (buoyancyReferenceDensity * gravity[i] *
                                     (coords[i] - referenceLocation[i])) /
                                    (Rs * Tb[iNode]);
                            }
                        }
                    }
                }
                break;

            case buoyancyModel::boussinesq:
                {
                    errorMsg("Must not reach here");
                }
                break;
        }
    }
}

void idealGasModel::updateCompressibility(const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* psiSTKFieldPtr = psiRef().stkFieldPtr();

    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

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
        scalar* psib = stk::mesh::field_data(*psiSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc psi
            psib[iNode] = 1.0 / (Rs * T);
        }
    }
}

void idealGasModel::updateCompressibility(const std::shared_ptr<domain> domain,
                                          label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* psiSTKFieldPtr = psiRef(iPhase).stkFieldPtr();

    scalar molarMass =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

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
        scalar* psib = stk::mesh::field_data(*psiSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            scalar T = Tb[iNode];

            // calc psi
            psib[iNode] = 1.0 / (Rs * T);
        }
    }
}

void idealGasModel::updateSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

    const auto& coeffs =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.coeffs_;

    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

    // Get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Define selectors: select owned nodes in the interior
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    const stk::mesh::BucketVector& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const stk::mesh::Bucket* b : nodeBuckets)
    {
        const stk::mesh::Bucket& nodeBucket = *b;
        const auto nNodesPerBucket = nodeBucket.size();

        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* cpb = stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);

        for (size_t iNode = 0; iNode < nNodesPerBucket; ++iNode)
        {
            scalar T = Tb[iNode];
            scalar Ti = 1.0;
            scalar cp = 0.0;

            for (size_t i = 0; i < coeffs.size(); ++i)
            {
                cp += coeffs[i] * Ti;
                Ti *= T;
            }

            cpb[iNode] = Rs * cp;
        }
    }
}

void idealGasModel::updateSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

    const auto& coeffs =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.coeffs_;

    label localPhaseIndexInDomain = domain->globalToLocalMaterialIndex(iPhase);

    scalar molarMass =
        domain->materialRef(localPhaseIndexInDomain)
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = universalGasConstant_ / molarMass;

    // Get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Define selectors: select owned nodes in the interior
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    const stk::mesh::BucketVector& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const stk::mesh::Bucket* b : nodeBuckets)
    {
        const stk::mesh::Bucket& nodeBucket = *b;
        const auto nNodesPerBucket = nodeBucket.size();

        const scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* cpb = stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);

        for (size_t iNode = 0; iNode < nNodesPerBucket; ++iNode)
        {
            scalar T = Tb[iNode];
            scalar Ti = 1.0;
            scalar cp = 0.0;

            for (size_t i = 0; i < coeffs.size(); ++i)
            {
                cp += coeffs[i] * Ti;
                Ti *= T;
            }

            cpb[iNode] = Rs * cp;
        }
    }
}

void idealGasModel::updateSpecificEnthalpy(const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();
    const STKScalarField* h0STKFieldPtr = h0Ref().stkFieldPtr();
    STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();

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
        const scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
        const scalar* h0b = stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);
        scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // caluclate UmagSqr
            scalar UmagSqr = 0.0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                UmagSqr +=
                    Ub[SPATIAL_DIM * iNode + i] * Ub[SPATIAL_DIM * iNode + i];
            }

            scalar h0 = h0b[iNode];

            // calc h
            hb[iNode] = h0 - 0.5 * UmagSqr;
        }
    }
}

void idealGasModel::updateSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    errorMsg("Must not reach here");
}

void idealGasModel::updateTemperature(const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Get fields
    const STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
    STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();

    scalar Tmin = TRef().minAcceptedValue();
    scalar Tmax = TRef().maxAcceptedValue();

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Define selector
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    const stk::mesh::BucketVector& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    auto option = domain->materialRef()
                      .thermodynamicProperties_.specificHeatCapacity_.option_;

    switch (option)
    {
        case specificHeatCapacityOption::value:
            {
                const STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

                for (const stk::mesh::Bucket* b : nodeBuckets)
                {
                    const stk::mesh::Bucket& nodeBucket = *b;
                    const auto nNodesPerBucket = nodeBucket.size();

                    const scalar* cpb =
                        stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
                    const scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
                    scalar* Tb =
                        stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);

                    for (size_t iNode = 0; iNode < nNodesPerBucket; ++iNode)
                    {
                        scalar h = hb[iNode];
                        scalar cp = cpb[iNode];
                        scalar T = h / cp;

                        if (T < 0.9 * Tb[iNode])
                            T = 0.9 * Tb[iNode];

                        T = std::max(std::min(T, Tmax), Tmin);
                        Tb[iNode] = T;
                    }
                }
            }
            break;

        case specificHeatCapacityOption::zeroPressurePolynomial:
            {
                const auto& coeffs =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
                scalar molarMass =
                    domain->materialRef()
                        .thermodynamicProperties_.equationOfState_.molarMass_;
                scalar Rs = universalGasConstant_ / molarMass;

                std::array<scalar, 9> hCoeffs = {0.0};
                for (size_t i = 0; i < 8; ++i)
                {
                    hCoeffs[i + 1] =
                        Rs * coeffs[i] / static_cast<scalar>(i + 1);
                }

                utils::polynomial<8> p8(hCoeffs);

                for (const stk::mesh::Bucket* b : nodeBuckets)
                {
                    const stk::mesh::Bucket& nodeBucket = *b;
                    const auto nNodesPerBucket = nodeBucket.size();

                    const scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
                    scalar* Tb =
                        stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);

                    for (size_t iNode = 0; iNode < nNodesPerBucket; ++iNode)
                    {
                        scalar h = hb[iNode];

                        scalar T = p8.solve(h, Tb[iNode]);

                        if (T < 0.9 * Tb[iNode])
                            T = 0.9 * Tb[iNode];

                        T = std::max(std::min(T, Tmax), Tmin);
                        Tb[iNode] = T;
                    }
                }
            }
            break;
    }
}

} /* namespace accel */
