// File       : heatTransferModel.cpp
// Created    : Thu Feb 01 2024 16:41:11 (+0100)
// Author     : Fabian Wermelinger
// Description: Heat transfer base class implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "heatTransferModel.h"
#include "domain.h"
#include "idealGasModel.h"
#include "kineticTheoryModel.h"
#include "polynomial.h"
#include "realm.h"
#include "simulation.h"
#include "sutherlandsFormulaModel.h"

namespace accel
{

heatTransferModel::heatTransferModel(realm* realm) : model(realm)
{
    // create field instances
    rhoRef();
    cpRef();
    lambdaRef();
    lambdaEffRef(); // always required even if not turbulent
    TRef();
    hRef();
    h0Ref();
    T0Ref();
    TWallCoeffsRef();
    qDotRef();
    URef();    // always required (for solid will be 0)
    pRef();    // always required (for solid will be 0)
    mDotRef(); // always required (for solid will be 0)

    // instantiate compressible fields
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->type() == domainType::fluid &&
            domain->isMaterialCompressible())
        {
            psiRef();
            break;
        }
    }

    // instantiate buoyancy-related fields (i.e. thermal expansivity)
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->type() == domainType::fluid &&
            !domain->isMaterialCompressible() &&
            domain->buoyancy_.option_ == buoyancyOption::buoyant &&
            domain->nMaterials() == 1)
        {
            betaRef();
            break;
        }
    }

    // instantiate phasic property fields in case of a multiphase
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->type() == domainType::fluid &&
            domain->multiphase_.homogeneous_ &&
            domain->multiphase_.freeSurfaceModel_.option_ ==
                freeSurfaceModelOption::standard)
        {
            for (const auto& material : domain->materialVector())
            {
                std::string materialName = material.name_;

                // global material index
                label phaseIndex =
                    realm->simulationRef().materialIndex(materialName);

                rhoRef(phaseIndex);
                cpRef(phaseIndex);
                lambdaRef(phaseIndex);
            }
        }
    }

    constexpr label count = 7; // number of fields in BoundaryData struct
    label block_length[count] = {1, 1, 1, 1, 1, 32, 32};
    MPI_Aint block_disp[count] = {offsetof(HeatBoundaryData, in),
                                  offsetof(HeatBoundaryData, out),
                                  offsetof(HeatBoundaryData, in_area),
                                  offsetof(HeatBoundaryData, out_area),
                                  offsetof(HeatBoundaryData, total_area),
                                  offsetof(HeatBoundaryData, name),
                                  offsetof(HeatBoundaryData, type)};
    MPI_Datatype block_type[count] = {MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_CHAR,
                                      MPI_CHAR};
    MPI_Datatype TempType;
    MPI_Type_create_struct(
        count, block_length, block_disp, block_type, &TempType);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(TempType, &lb, &extent);
    MPI_Type_create_resized(TempType, lb, extent, &MPIHeatBoundaryData);
    MPI_Type_commit(&MPIHeatBoundaryData);

    MPI_Op_create(
        heatTransferModel::sumHeatBoundaryData_, 1, &MPIHeatBoundaryData_SUM);
}

void heatTransferModel::updateTotalTemperatureField_(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get fields
    STKScalarField* T0STKFieldPtr = T0Ref().stkFieldPtr();
    const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();
    const STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

    auto eosOption =
        domain->materialRef().thermodynamicProperties_.equationOfState_.option_;

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    const stk::mesh::BucketVector& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    switch (eosOption)
    {
        case equationOfStateOption::value:
            {
                for (const stk::mesh::Bucket* b : nodeBuckets)
                {
                    const auto& bucket = *b;
                    const size_t n = bucket.size();

                    scalar* T0b = stk::mesh::field_data(*T0STKFieldPtr, bucket);
                    const scalar* Tb =
                        stk::mesh::field_data(*TSTKFieldPtr, bucket);
                    const scalar* Ub =
                        stk::mesh::field_data(*USTKFieldPtr, bucket);
                    const scalar* cpb =
                        stk::mesh::field_data(*cpSTKFieldPtr, bucket);

                    for (size_t iNode = 0; iNode < n; ++iNode)
                    {
                        scalar T = Tb[iNode];
                        scalar cp = cpb[iNode];
                        scalar UmagSqr = 0.0;

                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            scalar u = Ub[SPATIAL_DIM * iNode + d];
                            UmagSqr += u * u;
                        }

                        T0b[iNode] = T + 0.5 * UmagSqr / cp;
                    }
                }
            }
            break;

        case equationOfStateOption::idealGas:
            {
                auto cpOption =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                switch (cpOption)
                {
                    case specificHeatCapacityOption::value:
                        {
                            const STKScalarField* MaSTKFieldPtr =
                                MaRef().stkFieldPtr();
                            const STKScalarField* cpSTKFieldPtr =
                                cpRef().stkFieldPtr();

                            scalar M = domain->materialRef()
                                           .thermodynamicProperties_
                                           .equationOfState_.molarMass_;
                            scalar Rs = thermoModel::universalGasConstant_ / M;

                            for (const stk::mesh::Bucket* b : nodeBuckets)
                            {
                                const auto& bucket = *b;
                                const size_t n = bucket.size();

                                scalar* T0b = stk::mesh::field_data(
                                    *T0STKFieldPtr, bucket);
                                const scalar* Tb = stk::mesh::field_data(
                                    *TSTKFieldPtr, bucket);
                                const scalar* cpb = stk::mesh::field_data(
                                    *cpSTKFieldPtr, bucket);
                                const scalar* Mab = stk::mesh::field_data(
                                    *MaSTKFieldPtr, bucket);

                                for (size_t iNode = 0; iNode < n; ++iNode)
                                {
                                    scalar T = Tb[iNode];
                                    scalar cp = cpb[iNode];
                                    scalar Ma = Mab[iNode];

                                    scalar gamma = cp / (cp - Rs);
                                    T0b[iNode] = T * (1.0 + (gamma - 1.0) *
                                                                0.5 * Ma * Ma);
                                }
                            }
                        }
                        break;

                    case specificHeatCapacityOption::zeroPressurePolynomial:
                        {
                            const auto& coeffs =
                                domain->materialRef()
                                    .thermodynamicProperties_
                                    .specificHeatCapacity_.coeffs_;

                            scalar M = domain->materialRef()
                                           .thermodynamicProperties_
                                           .equationOfState_.molarMass_;
                            scalar Rs = thermoModel::universalGasConstant_ / M;

                            std::array<scalar, 9> hCoeffs = {0.0};
                            for (size_t i = 0; i < 8; ++i)
                            {
                                hCoeffs[i + 1] =
                                    Rs * coeffs[i] / static_cast<scalar>(i + 1);
                            }
                            utils::polynomial<8> p8(hCoeffs);

                            for (const stk::mesh::Bucket* b : nodeBuckets)
                            {
                                const auto& bucket = *b;
                                const size_t n = bucket.size();

                                scalar* T0b = stk::mesh::field_data(
                                    *T0STKFieldPtr, bucket);
                                const scalar* Tb = stk::mesh::field_data(
                                    *TSTKFieldPtr, bucket);
                                const scalar* Ub = stk::mesh::field_data(
                                    *USTKFieldPtr, bucket);

                                for (size_t iNode = 0; iNode < n; ++iNode)
                                {
                                    scalar T = Tb[iNode];
                                    scalar UmagSqr = 0.0;

                                    for (label d = 0; d < SPATIAL_DIM; ++d)
                                    {
                                        scalar u = Ub[SPATIAL_DIM * iNode + d];
                                        UmagSqr += u * u;
                                    }

                                    // Compute h(T) using progressive powers
                                    scalar h = 0.0;
                                    scalar Ti = T;
                                    for (size_t i = 1; i < hCoeffs.size(); ++i)
                                    {
                                        h += hCoeffs[i] * Ti;
                                        Ti *= T;
                                    }

                                    scalar h0 = h + 0.5 * UmagSqr;
                                    T0b[iNode] = p8.solve(h0, T);
                                }
                            }
                        }
                        break;
                }
            }
            break;
    }
}

void heatTransferModel::setupSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    // phasic
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        for (const auto& material : domain->materialVector())
        {
            std::string materialName = material.name_;

            // global material index
            label phaseIndex =
                this->realmRef().simulationRef().materialIndex(materialName);

            fieldBroker::setupSpecificHeatCapacity(domain, phaseIndex);
        }

        // bulk
        if (cpRef().isZoneUnset(domain->index()))
        {
            cpRef().setZone(domain->index());
        }
    }
    else
    {
        fieldBroker::setupSpecificHeatCapacity(domain);
    }
}

void heatTransferModel::setupThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    // phasic
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        for (const auto& material : domain->materialVector())
        {
            std::string materialName = material.name_;

            // global material index
            label phaseIndex =
                this->realmRef().simulationRef().materialIndex(materialName);

            fieldBroker::setupThermalConductivity(domain, phaseIndex);
        }

        // bulk
        if (lambdaRef().isZoneUnset(domain->index()))
        {
            lambdaRef().setZone(domain->index());
        }
    }
    else
    {
        fieldBroker::setupThermalConductivity(domain);
    }
}

void heatTransferModel::initializeSpecificEnthalpy(
    const std::shared_ptr<domain> domain)
{
    auto option =
        domain->materialRef().thermodynamicProperties_.equationOfState_.option_;

    switch (option)
    {
        case equationOfStateOption::value:
            {
                auto& mesh = this->meshRef();
                stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
                stk::mesh::MetaData& metaData = mesh.metaDataRef();

                // Get fields
                STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
                const STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
                const STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() & stk::mesh::selectUnion(partVec);

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* cpb =
                        stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
                    const scalar* Tb =
                        stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
                    scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        scalar cp = cpb[iNode];
                        scalar T = Tb[iNode];

                        // calc h
                        hb[iNode] = cp * T;
                    }
                }
            }
            break;

        case equationOfStateOption::idealGas:
            {
                auto option =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                switch (option)
                {
                    case specificHeatCapacityOption::value:
                        {
                            auto& mesh = this->meshRef();
                            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
                            stk::mesh::MetaData& metaData = mesh.metaDataRef();

                            // Get fields
                            STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
                            const STKScalarField* TSTKFieldPtr =
                                TRef().stkFieldPtr();
                            const STKScalarField* cpSTKFieldPtr =
                                cpRef().stkFieldPtr();

                            // get interior parts the domain is defined on
                            const stk::mesh::PartVector& partVec =
                                domain->zonePtr()->interiorParts();

                            // define some common selectors; select owned nodes
                            stk::mesh::Selector selUniversalNodes =
                                metaData.universal_part() &
                                stk::mesh::selectUnion(partVec);

                            stk::mesh::BucketVector const& nodeBuckets =
                                bulkData.get_buckets(stk::topology::NODE_RANK,
                                                     selUniversalNodes);
                            for (stk::mesh::BucketVector::const_iterator ib =
                                     nodeBuckets.begin();
                                 ib != nodeBuckets.end();
                                 ++ib)
                            {
                                stk::mesh::Bucket& nodeBucket = **ib;

                                const stk::mesh::Bucket::size_type
                                    nNodesPerBucket = nodeBucket.size();

                                // field chunks in bucket
                                const scalar* cpb = stk::mesh::field_data(
                                    *cpSTKFieldPtr, nodeBucket);
                                const scalar* Tb = stk::mesh::field_data(
                                    *TSTKFieldPtr, nodeBucket);
                                scalar* hb = stk::mesh::field_data(
                                    *hSTKFieldPtr, nodeBucket);

                                for (stk::mesh::Bucket::size_type iNode = 0;
                                     iNode < nNodesPerBucket;
                                     ++iNode)
                                {
                                    scalar cp = cpb[iNode];
                                    scalar T = Tb[iNode];

                                    // calc h
                                    hb[iNode] = cp * T;
                                }
                            }
                        }
                        break;

                    case specificHeatCapacityOption::zeroPressurePolynomial:
                        {
                            std::unique_ptr<idealGasModel> model =
                                std::make_unique<idealGasModel>(
                                    this->realmPtr_);
                            model->initializeSpecificEnthalpy(domain);
                        }
                        break;
                }
            }
            break;
    }

    // model-based side updates
    updateSpecificEnthalpySideFields_(domain);
}

void heatTransferModel::initializeSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get fields
    STKScalarField* h0STKFieldPtr = h0Ref().stkFieldPtr();
    const STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();

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
            scalar UmagSqr = 0.0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                UmagSqr +=
                    Ub[SPATIAL_DIM * iNode + i] * Ub[SPATIAL_DIM * iNode + i];
            }

            // calc h0
            h0b[iNode] = hb[iNode] + 0.5 * UmagSqr;
        }
    }

    // model-based side updates
    updateSpecificTotalEnthalpySideFields_(domain);
}

void heatTransferModel::initializeSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        // update phasic electrical conductivity
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            auto option =
                domain->materialRef(iPhase)
                    .thermodynamicProperties_.specificHeatCapacity_.option_;

            switch (option)
            {
                case specificHeatCapacityOption::value:
                    {
                        fieldBroker::initializeSpecificHeatCapacity(domain,
                                                                    phaseIndex);
                    }
                    break;

                case specificHeatCapacityOption::zeroPressurePolynomial:
                    {
                        std::unique_ptr<idealGasModel> model =
                            std::make_unique<idealGasModel>(this->realmPtr_);
                        model->initializeSpecificHeatCapacity(domain,
                                                              phaseIndex);
                    }
                    break;
            }
        }

        // update bulk specific heat capacity
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // Get pointer to global density field to be used
        const auto* bulkCpSTKFieldPtr = cpRef().stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Initialize global density to zero before every update
        cpRef().setToValue({0}, partVec);

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            // Get fields for a given iPhase
            const STKScalarField* cpSTKFieldPtr =
                this->cpRef(phaseIndex).stkFieldPtr();
            const STKScalarField* alphaSTKFieldPtr =
                this->alphaRef(phaseIndex).stkFieldPtr();

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* cpb =
                    stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
                const scalar* alphab =
                    stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
                scalar* bulkCpb =
                    stk::mesh::field_data(*bulkCpSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    bulkCpb[iNode] += cpb[iNode] * alphab[iNode];
                }
            }
        }
    }
    else
    {
        auto option =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.option_;

        switch (option)
        {
            case specificHeatCapacityOption::value:
                fieldBroker::initializeSpecificHeatCapacity(domain);
                break;

            case specificHeatCapacityOption::zeroPressurePolynomial:
                {
                    std::unique_ptr<idealGasModel> model =
                        std::make_unique<idealGasModel>(this->realmPtr_);
                    model->initializeSpecificHeatCapacity(domain);
                }
                break;
        }
    }
}

void heatTransferModel::initializeThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        // update phasic electrical conductivity
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            auto option =
                domain->materialRef(iPhase)
                    .transportProperties_.thermalConductivity_.option_;

            switch (option)
            {
                case thermalConductivityOption::value:
                    fieldBroker::initializeThermalConductivity(domain,
                                                               phaseIndex);
                    break;

                case thermalConductivityOption::kineticTheoryModel:
                    {
                        assert(domain->type() == domainType::fluid);

                        errorMsg("Not yet implemented");
                    }
                    break;

                case thermalConductivityOption::sutherlandsFormula:
                    {
                        assert(domain->type() == domainType::fluid);

                        errorMsg("Not yet implemented");
                    }
                    break;
            }
        }

        // update bulk electrical conductivity
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // Get pointer to global density field to be used
        const auto* bulkLambdaSTKFieldPtr = lambdaRef().stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Initialize global density to zero before every update
        lambdaRef().setToValue({0}, partVec);

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            // Get fields for a given iPhase
            const STKScalarField* lambdaSTKFieldPtr =
                this->lambdaRef(phaseIndex).stkFieldPtr();
            const STKScalarField* alphaSTKFieldPtr =
                this->alphaRef(phaseIndex).stkFieldPtr();

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* lambdab =
                    stk::mesh::field_data(*lambdaSTKFieldPtr, nodeBucket);
                const scalar* alphab =
                    stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
                scalar* bulkLambdab =
                    stk::mesh::field_data(*bulkLambdaSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    bulkLambdab[iNode] += lambdab[iNode] * alphab[iNode];
                }
            }
        }
    }
    else
    {
        auto option = domain->materialRef()
                          .transportProperties_.thermalConductivity_.option_;

        switch (option)
        {
            case thermalConductivityOption::value:
                fieldBroker::initializeThermalConductivity(domain);
                break;

            case thermalConductivityOption::kineticTheoryModel:
                {
                    std::unique_ptr<kineticTheoryModel> model =
                        std::make_unique<kineticTheoryModel>(this->realmPtr_);
                    model->initializeThermalConductivity(domain);
                }
                break;

            case thermalConductivityOption::sutherlandsFormula:
                {
                    assert(domain->type() == domainType::fluid);

                    std::unique_ptr<sutherlandsFormulaModel> model =
                        std::make_unique<sutherlandsFormulaModel>(
                            this->realmPtr_);
                    model->initializeThermalConductivity(domain);
                }
                break;
        }
    }
}

void heatTransferModel::initializeCompressibility(
    const std::shared_ptr<domain> domain)
{
    auto option =
        domain->materialRef().thermodynamicProperties_.equationOfState_.option_;

    switch (option)
    {
        case equationOfStateOption::value:
            break;

        case equationOfStateOption::idealGas:
            {
                std::unique_ptr<idealGasModel> model =
                    std::make_unique<idealGasModel>(this->realmPtr_);
                model->initializeCompressibility(domain);
            }
            break;
    }
}

void heatTransferModel::updateTemperature(const std::shared_ptr<domain> domain)
{
    fieldBroker::updateTemperature(domain);

    auto option =
        domain->materialRef().thermodynamicProperties_.equationOfState_.option_;

    switch (option)
    {
        case equationOfStateOption::value:
            {
                auto& mesh = this->meshRef();
                stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
                stk::mesh::MetaData& metaData = mesh.metaDataRef();

                // Get fields
                STKScalarField* TSTKFieldPtr = TRef().stkFieldPtr();
                const STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
                const STKScalarField* cpSTKFieldPtr = cpRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() & stk::mesh::selectUnion(partVec);

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* cpb =
                        stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
                    const scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
                    scalar* Tb =
                        stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        Tb[iNode] = hb[iNode] / cpb[iNode];
                    }
                }
            }
            break;

        case equationOfStateOption::idealGas:
            {
                std::unique_ptr<idealGasModel> model =
                    std::make_unique<idealGasModel>(this->realmPtr_);
                model->updateTemperature(domain);
            }
            break;
    }

    // model-based side updates
    updateTemperatureSideFields_(domain);
}

void heatTransferModel::updateSpecificEnthalpy(
    const std::shared_ptr<domain> domain)
{
    auto heatTransferOption = domain->heatTransfer_.option_;

    switch (heatTransferOption)
    {
        case heatTransferOption::thermalEnergy:
            {
                // model-based side updates
                updateSpecificEnthalpySideFields_(domain);
            }
            break;

        case heatTransferOption::totalEnergy:
            {
                const auto& mesh = meshRef();
                const stk::mesh::MetaData& metaData = mesh.metaDataRef();
                const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

                // Get fields
                const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();
                const STKScalarField* h0STKFieldPtr = h0Ref().stkFieldPtr();
                STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() & stk::mesh::selectUnion(partVec);

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* Ub =
                        stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
                    const scalar* h0b =
                        stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);
                    scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        // caluclate UmagSqr
                        scalar UmagSqr = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            UmagSqr += Ub[SPATIAL_DIM * iNode + i] *
                                       Ub[SPATIAL_DIM * iNode + i];
                        }

                        scalar h0 = h0b[iNode];

                        // calc h
                        hb[iNode] = h0 - 0.5 * UmagSqr;
                    }
                }
            }
            break;

        default:
            break;
    }
}

void heatTransferModel::updateSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    auto heatTransferOption = domain->heatTransfer_.option_;

    switch (heatTransferOption)
    {
        case heatTransferOption::thermalEnergy:
            {
                const auto& mesh = meshRef();
                const stk::mesh::MetaData& metaData = mesh.metaDataRef();
                const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

                // Get fields
                const STKScalarField* USTKFieldPtr = URef().stkFieldPtr();
                const STKScalarField* hSTKFieldPtr = hRef().stkFieldPtr();
                STKScalarField* h0STKFieldPtr = h0Ref().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() & stk::mesh::selectUnion(partVec);

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* Ub =
                        stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
                    const scalar* hb =
                        stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
                    scalar* h0b =
                        stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        // caluclate UmagSqr
                        scalar UmagSqr = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            UmagSqr += Ub[SPATIAL_DIM * iNode + i] *
                                       Ub[SPATIAL_DIM * iNode + i];
                        }

                        scalar h = hb[iNode];

                        // calc h0
                        h0b[iNode] = h + 0.5 * UmagSqr;
                    }
                }
            }
            break;

        case heatTransferOption::totalEnergy:
            {
                // model-based side updates
                updateSpecificTotalEnthalpySideFields_(domain);
            }
            break;

        default:
            break;
    }
}

void heatTransferModel::updateSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        // update phasic electrical conductivity
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            auto option =
                domain->materialRef(iPhase)
                    .thermodynamicProperties_.specificHeatCapacity_.option_;

            switch (option)
            {
                case specificHeatCapacityOption::value:
                    {
                        fieldBroker::updateSpecificHeatCapacity(domain,
                                                                phaseIndex);
                    }
                    break;

                case specificHeatCapacityOption::zeroPressurePolynomial:
                    {
                        std::unique_ptr<idealGasModel> model =
                            std::make_unique<idealGasModel>(this->realmPtr_);
                        model->updateSpecificHeatCapacity(domain, phaseIndex);
                    }
                    break;
            }
        }

        // update bulk electrical conductivity
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // Get pointer to global density field to be used
        const auto* bulkCpSTKFieldPtr = cpRef().stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Initialize global density to zero before every update
        cpRef().setToValue({0}, partVec);

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            // Get fields for a given iPhase
            const STKScalarField* cpSTKFieldPtr =
                this->cpRef(phaseIndex).stkFieldPtr();
            const STKScalarField* alphaSTKFieldPtr =
                this->alphaRef(phaseIndex).stkFieldPtr();

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* cpb =
                    stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
                const scalar* alphab =
                    stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
                scalar* bulkCpb =
                    stk::mesh::field_data(*bulkCpSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    bulkCpb[iNode] += cpb[iNode] * alphab[iNode];
                }
            }
        }
    }
    else
    {
        auto option =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.option_;

        switch (option)
        {
            case specificHeatCapacityOption::value:
                fieldBroker::updateSpecificHeatCapacity(domain);
                break;

            case specificHeatCapacityOption::zeroPressurePolynomial:
                {
                    std::unique_ptr<idealGasModel> model =
                        std::make_unique<idealGasModel>(this->realmPtr_);
                    model->updateSpecificHeatCapacity(domain);
                }
                break;
        }
    }
}

void heatTransferModel::updateThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        // update phasic thermal conductivity
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            auto option =
                domain->materialRef(iPhase)
                    .transportProperties_.thermalConductivity_.option_;

            switch (option)
            {
                case thermalConductivityOption::value:
                    {
                        fieldBroker::updateThermalConductivity(domain,
                                                               phaseIndex);
                    }
                    break;

                case thermalConductivityOption::kineticTheoryModel:
                    {
                        assert(domain->type() == domainType::fluid);

                        std::unique_ptr<kineticTheoryModel> model =
                            std::make_unique<kineticTheoryModel>(
                                this->realmPtr_);
                        model->updateThermalConductivity(domain, phaseIndex);
                    }
                    break;

                case thermalConductivityOption::sutherlandsFormula:
                    {
                        assert(domain->type() == domainType::fluid);

                        std::unique_ptr<sutherlandsFormulaModel> model =
                            std::make_unique<sutherlandsFormulaModel>(
                                this->realmPtr_);
                        model->updateThermalConductivity(domain, phaseIndex);
                    }
                    break;
            }
        }

        // update bulk thermal conductivity
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // Get pointer to global density field to be used
        const auto* bulkLambdaSTKFieldPtr = lambdaRef().stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Initialize global density to zero before every update
        lambdaRef().setToValue({0}, partVec);

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            // Get fields for a given iPhase
            const STKScalarField* lambdaSTKFieldPtr =
                this->lambdaRef(phaseIndex).stkFieldPtr();
            const STKScalarField* alphaSTKFieldPtr =
                this->alphaRef(phaseIndex).stkFieldPtr();

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* lambdab =
                    stk::mesh::field_data(*lambdaSTKFieldPtr, nodeBucket);
                const scalar* alphab =
                    stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
                scalar* bulkLambdab =
                    stk::mesh::field_data(*bulkLambdaSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    bulkLambdab[iNode] += lambdab[iNode] * alphab[iNode];
                }
            }
        }
    }
    else
    {
        auto option = domain->materialRef()
                          .transportProperties_.thermalConductivity_.option_;

        switch (option)
        {
            case thermalConductivityOption::value:
                fieldBroker::updateThermalConductivity(domain);
                break;

            case thermalConductivityOption::kineticTheoryModel:
                {
                    std::unique_ptr<kineticTheoryModel> model =
                        std::make_unique<kineticTheoryModel>(this->realmPtr_);
                    model->updateThermalConductivity(domain);
                }
                break;

            case thermalConductivityOption::sutherlandsFormula:
                {
                    std::unique_ptr<sutherlandsFormulaModel> model =
                        std::make_unique<sutherlandsFormulaModel>(
                            this->realmPtr_);
                    model->updateThermalConductivity(domain);
                }
                break;
        }
    }
}

void heatTransferModel::updateCompressibility(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::fluid &&
        domain->multiphase_.homogeneous_ &&
        domain->multiphase_.freeSurfaceModel_.option_ ==
            freeSurfaceModelOption::standard)
    {
        // update phasic thermal conductivity
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            auto option =
                domain->materialRef(iPhase)
                    .thermodynamicProperties_.equationOfState_.option_;

            switch (option)
            {
                case equationOfStateOption::value:
                    break;

                case equationOfStateOption::idealGas:
                    {
                        std::unique_ptr<idealGasModel> model =
                            std::make_unique<idealGasModel>(this->realmPtr_);
                        model->updateCompressibility(domain, phaseIndex);
                    }
                    break;
            }
        }

        // update bulk compressibility
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // Get pointer to global density field to be used
        const auto* bulkPsiSTKFieldPtr = psiRef().stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Initialize global density to zero before every update
        psiRef().setToValue({0}, partVec);

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            // Get fields for a given iPhase
            const STKScalarField* psiSTKFieldPtr =
                this->psiRef(phaseIndex).stkFieldPtr();
            const STKScalarField* alphaSTKFieldPtr =
                this->alphaRef(phaseIndex).stkFieldPtr();

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* psib =
                    stk::mesh::field_data(*psiSTKFieldPtr, nodeBucket);
                const scalar* alphab =
                    stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
                scalar* bulkPsib =
                    stk::mesh::field_data(*bulkPsiSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    bulkPsib[iNode] += psib[iNode] * alphab[iNode];
                }
            }
        }
    }
    else
    {
        auto option = domain->materialRef()
                          .thermodynamicProperties_.equationOfState_.option_;

        switch (option)
        {
            case equationOfStateOption::value:
                break;

            case equationOfStateOption::idealGas:
                {
                    std::unique_ptr<idealGasModel> model =
                        std::make_unique<idealGasModel>(this->realmPtr_);
                    model->updateCompressibility(domain);
                }
                break;
        }
    }
}

void heatTransferModel::updateEffectiveThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (domain->type() == domainType::solid ||
        (domain->turbulence_.option_ == turbulenceOption::laminar &&
         domain->type() == domainType::fluid))
    {
        const auto& mesh = this->meshRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();

        const STKScalarField* lambdaSTKFieldPtr =
            this->lambdaRef().stkFieldPtr();
        STKScalarField* lambdaEffSTKFieldPtr =
            this->lambdaEffRef().stkFieldPtr();

        // define some common selectors
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();

            const scalar* lambdab =
                stk::mesh::field_data(*lambdaSTKFieldPtr, b);
            scalar* lambdaEffb =
                stk::mesh::field_data(*lambdaEffSTKFieldPtr, b);

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                lambdaEffb[k] = lambdab[k];
            }
        }
    }
}

void heatTransferModel::updateSpecificHeatCapacityGradientField(
    const std::shared_ptr<domain> domain)
{
    auto option = domain->materialRef()
                      .thermodynamicProperties_.specificHeatCapacity_.option_;

    switch (option)
    {
        case specificHeatCapacityOption::value:
            {
                // no grad needed: skip to save time
            }
            break;

        case specificHeatCapacityOption::zeroPressurePolynomial:
            {
                fieldBroker::updateSpecificHeatCapacityGradientField(domain);
            }
            break;
    }
}

void heatTransferModel::updateTemperatureSideFields_(
    const std::shared_ptr<domain> domain)
{
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();

        const auto& bcType =
            TRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::totalTemperature:
                            {
                                updateTemperatureBoundarySideFieldInletTotalTemperature_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::totalTemperature:
                            {
                                updateTemperatureBoundarySideFieldOpening_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            default:
                break;
        }
    }
}

void heatTransferModel::updateSpecificEnthalpySideFields_(
    const std::shared_ptr<domain> domain)
{
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();

        const auto& bcType =
            TRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                updateSpecificEnthalpyBoundarySideFieldWallSpecifiedValue_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                updateSpecificEnthalpyBoundarySideFieldInletSpecifiedValue_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                updateSpecificEnthalpyBoundarySideFieldOpening_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            default:
                break;
        }
    }
}

void heatTransferModel::updateSpecificTotalEnthalpySideFields_(
    const std::shared_ptr<domain> domain)
{
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();

        const auto& bcType =
            TRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                updateSpecificTotalEnthalpyBoundarySideFieldWallSpecifiedValue_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                updateSpecificTotalEnthalpyBoundarySideFieldInletSpecifiedValue_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::staticTemperature:
                        case boundaryConditionType::totalTemperature:
                            {
                                updateSpecificTotalEnthalpyBoundarySideFieldOpening_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            default:
                break;
        }
    }
}

void heatTransferModel::
    updateTemperatureBoundarySideFieldInletTotalTemperature_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    auto& bc = TRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& data = bc.data<1>("value");

    const bool correctedBoundaryNodeValues =
        TRef().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    const auto& coeffs =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = thermoModel::universalGasConstant_ / molarMass;

    std::array<scalar, 9> hCoeffs = {0.0};
    for (size_t i = 0; i < 8; ++i)
        hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

    utils::polynomial<8> p8(hCoeffs);

    auto computeStaticTemperatureFromTotalEnthalpy =
        [&](scalar T0, scalar UmagSqr, scalar Tinit) -> scalar
    {
        scalar Ti = T0;
        scalar h0 = 0.0;
        for (size_t i = 1; i < hCoeffs.size(); ++i)
        {
            h0 += hCoeffs[i] * Ti;
            Ti *= T0;
        }
        scalar rhs = h0 - 0.5 * UmagSqr;
        return p8.solve(rhs, Tinit);
    };

    switch (data.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
        case inputDataType::timeTable:
            {
                scalar T0Value = (data.type() == inputDataType::constant)
                                     ? *data.value()
                                     : data.interpolate(controlsRef().time)[0];

                const auto& cpOption =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                const auto& USTKFieldRef = URef().stkFieldRef();
                const auto& cpSTKFieldRef = cpRef().stkFieldRef();
                auto& nodeSideTSTKFieldRef =
                    TRef().nodeSideFieldRef().stkFieldRef();

                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                const auto& sideNodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

                for (const auto* b : sideNodeBuckets)
                {
                    const auto& bucket = *b;
                    const size_t n = bucket.size();

                    scalar* Tb =
                        stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
                    const scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, bucket);
                    const scalar* cpb =
                        (cpOption == specificHeatCapacityOption::value)
                            ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                            : nullptr;

                    for (size_t iNode = 0; iNode < n; ++iNode)
                    {
                        scalar UmagSqr = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            scalar u = Ub[iNode * SPATIAL_DIM + d];
                            UmagSqr += u * u;
                        }

                        if (cpOption == specificHeatCapacityOption::value)
                        {
                            Tb[iNode] = T0Value - 0.5 * UmagSqr / cpb[iNode];
                        }
                        else
                        {
                            Tb[iNode] =
                                computeStaticTemperatureFromTotalEnthalpy(
                                    T0Value, UmagSqr, Tb[iNode]);
                        }
                    }
                }

                // Interpolate to side field
                TRef().sideFieldRef().interpolate(TRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  TRef().isShifted());

                if (correctedBoundaryNodeValues)
                {
                    TRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        case inputDataType::expression:
            {
                auto cpOption =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                typedef exprtk::symbol_table<scalar> symbol_table_t;
                typedef exprtk::expression<scalar> expression_t;
                typedef exprtk::parser<scalar> parser_t;

                scalar t, x, y, z;
                symbol_table_t symbol_table;
                symbol_table.add_constants();
                symbol_table.add_variable("t", t);
                symbol_table.add_variable("x", x);
                symbol_table.add_variable("y", y);
                symbol_table.add_variable("z", z);

                expression_t expression;
                expression.register_symbol_table(symbol_table);

                parser_t parser;
                if (!parser.compile(data.expression()[0], expression))
                    errorMsg("Error in expression for total temperature");

                t = controlsRef().time;

                const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
                const auto& USTKFieldRef = URef().stkFieldRef();
                const auto& cpSTKFieldRef = cpRef().stkFieldRef();
                auto& nodeSideTSTKFieldRef =
                    TRef().nodeSideFieldRef().stkFieldRef();

                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                const auto& sideNodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

                for (const auto* b : sideNodeBuckets)
                {
                    const auto& bucket = *b;
                    const size_t n = bucket.size();

                    scalar* Tb =
                        stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
                    const scalar* coordsb =
                        stk::mesh::field_data(coordsSTKFieldRef, bucket);
                    const scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, bucket);
                    const scalar* cpb =
                        (cpOption == specificHeatCapacityOption::value)
                            ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                            : nullptr;

                    for (size_t iNode = 0; iNode < n; ++iNode)
                    {
                        scalar UmagSqr = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            scalar u = Ub[iNode * SPATIAL_DIM + d];
                            UmagSqr += u * u;
                        }

#if SPATIAL_DIM == 3
                        x = coordsb[iNode * SPATIAL_DIM + 0];
                        y = coordsb[iNode * SPATIAL_DIM + 1];
                        z = coordsb[iNode * SPATIAL_DIM + 2];
#elif SPATIAL_DIM == 2
                        x = coordsb[iNode * SPATIAL_DIM + 0];
                        y = coordsb[iNode * SPATIAL_DIM + 1];
#endif

                        scalar T0Value = expression.value();

                        if (cpOption == specificHeatCapacityOption::value)
                        {
                            Tb[iNode] = T0Value - 0.5 * UmagSqr / cpb[iNode];
                        }
                        else
                        {
                            Tb[iNode] =
                                computeStaticTemperatureFromTotalEnthalpy(
                                    T0Value, UmagSqr, Tb[iNode]);
                        }
                    }
                }

                // Interpolate to side field
                TRef().sideFieldRef().interpolate(TRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  TRef().isShifted());

                if (correctedBoundaryNodeValues)
                {
                    TRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        case inputDataType::profileData:
            {
                errorMsg("profile data not provided yet");
            }
            break;
    }
}

void heatTransferModel::updateTemperatureBoundarySideFieldOpening_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    auto& bc = TRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& data = bc.data<1>("value");

    const bool correctedBoundaryNodeValues =
        TRef().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    const auto& coeffs =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
    scalar molarMass =
        domain->materialRef()
            .thermodynamicProperties_.equationOfState_.molarMass_;
    scalar Rs = thermoModel::universalGasConstant_ / molarMass;

    std::array<scalar, 9> hCoeffs = {0.0};
    for (size_t i = 0; i < 8; ++i)
        hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

    utils::polynomial<8> p8(hCoeffs);

    auto computeStaticTemperatureFromTotalEnthalpy =
        [&](scalar T0, scalar UmagSqr, scalar Tinit) -> scalar
    {
        scalar Ti = T0;
        scalar h0 = 0.0;
        for (size_t i = 1; i < hCoeffs.size(); ++i)
        {
            h0 += hCoeffs[i] * Ti;
            Ti *= T0;
        }
        scalar rhs = h0 - 0.5 * UmagSqr;
        return p8.solve(rhs, Tinit);
    };

    switch (data.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
        case inputDataType::timeTable:
            {
                scalar T0Value = (data.type() == inputDataType::constant)
                                     ? *data.value()
                                     : data.interpolate(controlsRef().time)[0];

                const auto& cpOption =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                const auto& USTKFieldRef = URef().stkFieldRef();
                const auto& cpSTKFieldRef = cpRef().stkFieldRef();
                auto& nodeSideTSTKFieldRef =
                    TRef().nodeSideFieldRef().stkFieldRef();

                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                const auto& sideNodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

                for (const auto* b : sideNodeBuckets)
                {
                    const auto& bucket = *b;
                    const size_t n = bucket.size();

                    scalar* Tb =
                        stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
                    const scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, bucket);
                    const scalar* cpb =
                        (cpOption == specificHeatCapacityOption::value)
                            ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                            : nullptr;

                    for (size_t iNode = 0; iNode < n; ++iNode)
                    {
                        scalar UmagSqr = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            scalar u = Ub[iNode * SPATIAL_DIM + d];
                            UmagSqr += u * u;
                        }

                        if (cpOption == specificHeatCapacityOption::value)
                        {
                            Tb[iNode] = T0Value - 0.5 * UmagSqr / cpb[iNode];
                        }
                        else
                        {
                            Tb[iNode] =
                                computeStaticTemperatureFromTotalEnthalpy(
                                    T0Value, UmagSqr, Tb[iNode]);
                        }
                    }
                }

                // Interpolate to side field
                TRef().sideFieldRef().interpolate(TRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  TRef().isShifted());

                if (correctedBoundaryNodeValues)
                {
                    TRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        case inputDataType::expression:
            {
                auto cpOption =
                    domain->materialRef()
                        .thermodynamicProperties_.specificHeatCapacity_.option_;

                typedef exprtk::symbol_table<scalar> symbol_table_t;
                typedef exprtk::expression<scalar> expression_t;
                typedef exprtk::parser<scalar> parser_t;

                scalar t, x, y, z;
                symbol_table_t symbol_table;
                symbol_table.add_constants();
                symbol_table.add_variable("t", t);
                symbol_table.add_variable("x", x);
                symbol_table.add_variable("y", y);
                symbol_table.add_variable("z", z);

                expression_t expression;
                expression.register_symbol_table(symbol_table);

                parser_t parser;
                if (!parser.compile(data.expression()[0], expression))
                    errorMsg("Error in expression for total temperature");

                t = controlsRef().time;

                const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
                const auto& USTKFieldRef = URef().stkFieldRef();
                const auto& cpSTKFieldRef = cpRef().stkFieldRef();
                auto& nodeSideTSTKFieldRef =
                    TRef().nodeSideFieldRef().stkFieldRef();

                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                const auto& sideNodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

                for (const auto* b : sideNodeBuckets)
                {
                    const auto& bucket = *b;
                    const size_t n = bucket.size();

                    scalar* Tb =
                        stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
                    const scalar* coordsb =
                        stk::mesh::field_data(coordsSTKFieldRef, bucket);
                    const scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, bucket);
                    const scalar* cpb =
                        (cpOption == specificHeatCapacityOption::value)
                            ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                            : nullptr;

                    for (size_t iNode = 0; iNode < n; ++iNode)
                    {
                        scalar UmagSqr = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            scalar u = Ub[iNode * SPATIAL_DIM + d];
                            UmagSqr += u * u;
                        }

#if SPATIAL_DIM == 3
                        x = coordsb[iNode * SPATIAL_DIM + 0];
                        y = coordsb[iNode * SPATIAL_DIM + 1];
                        z = coordsb[iNode * SPATIAL_DIM + 2];
#elif SPATIAL_DIM == 2
                        x = coordsb[iNode * SPATIAL_DIM + 0];
                        y = coordsb[iNode * SPATIAL_DIM + 1];
#endif

                        scalar T0Value = expression.value();

                        if (cpOption == specificHeatCapacityOption::value)
                        {
                            Tb[iNode] = T0Value - 0.5 * UmagSqr / cpb[iNode];
                        }
                        else
                        {
                            Tb[iNode] =
                                computeStaticTemperatureFromTotalEnthalpy(
                                    T0Value, UmagSqr, Tb[iNode]);
                        }
                    }
                }

                // Interpolate to side field
                TRef().sideFieldRef().interpolate(TRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  TRef().isShifted());

                if (correctedBoundaryNodeValues)
                {
                    TRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        case inputDataType::profileData:
            {
                errorMsg("profile data not provided yet");
            }
            break;
    }
}

void heatTransferModel::
    updateSpecificEnthalpyBoundarySideFieldWallSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        hRef().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    std::array<scalar, 9> hCoeffs = {0.0};
    utils::polynomial<8> hPoly;

    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideHSTKFieldRef = hRef().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* hb = stk::mesh::field_data(nodeSideHSTKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            hb[i] = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
        }
    }

    // Interpolate to side field
    hRef().sideFieldRef().interpolate(hRef().nodeSideFieldRef(),
                                      domain->index(),
                                      boundary->index(),
                                      hRef().isShifted());

    if (correctedBoundaryNodeValues)
    {
        hRef().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::
    updateSpecificEnthalpyBoundarySideFieldInletSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        hRef().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    std::array<scalar, 9> hCoeffs = {0.0};
    utils::polynomial<8> hPoly;

    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideHSTKFieldRef = hRef().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* hb = stk::mesh::field_data(nodeSideHSTKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            hb[i] = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
        }
    }

    // Interpolate to side field
    hRef().sideFieldRef().interpolate(hRef().nodeSideFieldRef(),
                                      domain->index(),
                                      boundary->index(),
                                      hRef().isShifted());

    if (correctedBoundaryNodeValues)
    {
        hRef().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::updateSpecificEnthalpyBoundarySideFieldOpening_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        hRef().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    std::array<scalar, 9> hCoeffs = {0.0};
    utils::polynomial<8> hPoly;

    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideHSTKFieldRef = hRef().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* hb = stk::mesh::field_data(nodeSideHSTKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            hb[i] = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
        }
    }

    // Interpolate to side field
    hRef().sideFieldRef().interpolate(hRef().nodeSideFieldRef(),
                                      domain->index(),
                                      boundary->index(),
                                      hRef().isShifted());

    if (correctedBoundaryNodeValues)
    {
        hRef().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::
    updateSpecificTotalEnthalpyBoundarySideFieldWallSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        h0Ref().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    std::array<scalar, 9> hCoeffs = {0.0};
    utils::polynomial<8> hPoly;

    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& USTKFieldRef = domain->type() == domainType::fluid
                                   ? URef().nodeSideFieldRef().stkFieldRef()
                                   : URef().stkFieldRef(); // for solid U is 0
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideH0STKFieldRef = h0Ref().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* h0b = stk::mesh::field_data(nodeSideH0STKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* Ub = stk::mesh::field_data(USTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            scalar UmagSqr = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                scalar u = Ub[i * SPATIAL_DIM + d];
                UmagSqr += u * u;
            }

            scalar h = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
            h0b[i] = h + 0.5 * UmagSqr;
        }
    }

    // Interpolate to side field
    h0Ref().sideFieldRef().interpolate(h0Ref().nodeSideFieldRef(),
                                       domain->index(),
                                       boundary->index(),
                                       h0Ref().isShifted());

    if (correctedBoundaryNodeValues)
    {
        h0Ref().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::
    updateSpecificTotalEnthalpyBoundarySideFieldInletSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        h0Ref().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    std::array<scalar, 9> hCoeffs = {0.0};
    utils::polynomial<8> hPoly;

    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& USTKFieldRef = URef().stkFieldRef();
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideH0STKFieldRef = h0Ref().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* h0b = stk::mesh::field_data(nodeSideH0STKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* Ub = stk::mesh::field_data(USTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            scalar UmagSqr = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                scalar u = Ub[i * SPATIAL_DIM + d];
                UmagSqr += u * u;
            }

            scalar h = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
            h0b[i] = h + 0.5 * UmagSqr;
        }
    }

    // Interpolate to side field
    h0Ref().sideFieldRef().interpolate(h0Ref().nodeSideFieldRef(),
                                       domain->index(),
                                       boundary->index(),
                                       h0Ref().isShifted());

    if (correctedBoundaryNodeValues)
    {
        h0Ref().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::updateSpecificTotalEnthalpyBoundarySideFieldOpening_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    const bool correctedBoundaryNodeValues =
        h0Ref().correctedBoundaryNodeValues();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    const auto& cpOption =
        domain->materialRef()
            .thermodynamicProperties_.specificHeatCapacity_.option_;

    // Setup polynomial coefficients for h(T) if needed
    utils::polynomial<8> hPoly;
    if (cpOption == specificHeatCapacityOption::zeroPressurePolynomial)
    {
        const auto& coeffs =
            domain->materialRef()
                .thermodynamicProperties_.specificHeatCapacity_.coeffs_;
        scalar molarMass =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        scalar Rs = thermoModel::universalGasConstant_ / molarMass;

        std::array<scalar, 9> hCoeffs = {0.0};
        for (size_t i = 0; i < 8; ++i)
            hCoeffs[i + 1] = Rs * coeffs[i] / static_cast<scalar>(i + 1);

        hPoly = utils::polynomial<8>(hCoeffs);
    }

    // Define how to compute static enthalpy
    auto computeStaticEnthalpy = [&](scalar T, scalar cp) -> scalar
    {
        if (cpOption == specificHeatCapacityOption::value)
            return cp * T;
        else
            return hPoly.evaluate(T);
    };

    // Access fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& USTKFieldRef = URef().stkFieldRef();
    const auto& nodeSideTSTKFieldRef = TRef().nodeSideFieldRef().stkFieldRef();
    auto& nodeSideH0STKFieldRef = h0Ref().nodeSideFieldRef().stkFieldRef();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (const auto* b : sideNodeBuckets)
    {
        const auto& bucket = *b;
        const size_t n = bucket.size();

        scalar* h0b = stk::mesh::field_data(nodeSideH0STKFieldRef, bucket);
        const scalar* Tb = stk::mesh::field_data(nodeSideTSTKFieldRef, bucket);
        const scalar* Ub = stk::mesh::field_data(USTKFieldRef, bucket);
        const scalar* cpb = (cpOption == specificHeatCapacityOption::value)
                                ? stk::mesh::field_data(cpSTKFieldRef, bucket)
                                : nullptr;

        for (size_t i = 0; i < n; ++i)
        {
            scalar UmagSqr = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                scalar u = Ub[i * SPATIAL_DIM + d];
                UmagSqr += u * u;
            }

            scalar h = computeStaticEnthalpy(Tb[i], cpb ? cpb[i] : 0.0);
            h0b[i] = h + 0.5 * UmagSqr;
        }
    }

    // Interpolate to side field
    h0Ref().sideFieldRef().interpolate(h0Ref().nodeSideFieldRef(),
                                       domain->index(),
                                       boundary->index(),
                                       h0Ref().isShifted());

    if (correctedBoundaryNodeValues)
    {
        h0Ref().correctBoundaryNodes(domain->index(), boundary->index());
    }
}

void heatTransferModel::reportHeatData_()
{
#ifdef HAS_INTERFACE
    // Mass imbalance for interfaces
    if (this->meshRef().hasInterfaces())
    {
        scalar in = 0.0;
        scalar out = 0.0;
        if (messager::master())
        {
            for (const HeatBoundaryData& interfaceData :
                 heatInterfaceDataVector_)
            {
                in += interfaceData.in;
                out += interfaceData.out;
            }

            std::cout << std::endl;
            std::cout << "Interface Heat in:    " << std::scientific
                      << std::setprecision(3) << std::setw(10) << std::right
                      << in << '\n';
            std::cout << "Interface Heat out:   " << std::scientific
                      << std::setprecision(3) << std::setw(10) << std::right
                      << out << '\n';
            std::cout << "Interface Heat imbalance: " << std::scientific
                      << std::setprecision(3) << std::setw(10) << std::right
                      << (in + out) / (in + ::accel::SMALL) * 100 << " %\n";
        }

        // destroy interface data vector
        heatInterfaceDataVector_.clear();
    }
#endif /* HAS_INTERFACE */
}

void heatTransferModel::updateHeatImbalance_(
    const std::shared_ptr<domain> domain)
{
}

#ifdef HAS_INTERFACE
void heatTransferModel::updateInterfaceHeatImbalance_(
    const std::shared_ptr<domain> domain)
{
    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const auto& qDotSideSTKFieldRef = qDotRef().sideFieldRef().stkFieldRef();
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    std::vector<HeatBoundaryData> interfaceDataVector;

    bool noInterfaces = true;
    for (const interface* interf : domain->zonePtr()->interfacesRef())
    {
        noInterfaces = false;
        if (interf->isInternal())
        {
            // Master
            {
                // get interface side that is sitting in this domain
                const auto* interfaceSideInfoPtr = interf->masterInfoPtr();

                interfaceDataVector.push_back(HeatBoundaryData());
                HeatBoundaryData& interfaceData = interfaceDataVector.back();

                // set primary traits
                std::strncpy(interfaceData.name,
                             interfaceSideInfoPtr->name().c_str(),
                             sizeof(interfaceData.name));
                interfaceData.name[sizeof(interfaceData.name) - 1] = '\0';

                auto type = interf->type();
                std::strncpy(interfaceData.type,
                             toString(type).c_str(),
                             sizeof(interfaceData.type));
                interfaceData.type[sizeof(interfaceData.type) - 1] = '\0';

                // consider only if qDot is defined on the boundary
                if (!this->qDotRef().sideFieldRef().definedOn(
                        interfaceSideInfoPtr->currentPartVec_))
                {
                    continue;
                }

                // extract vector of dgInfo
                const std::vector<std::vector<dgInfo*>>& dgInfoVec =
                    interfaceSideInfoPtr->dgInfoVec_;

                for (label iSide = 0;
                     iSide < static_cast<label>(dgInfoVec.size());
                     iSide++)
                {
                    const std::vector<dgInfo*>& faceDgInfoVec =
                        dgInfoVec[iSide];

                    // now loop over all the DgInfo objects on this
                    // particular exposed face
                    for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
                    {
                        dgInfo* dgInfo = faceDgInfoVec[k];

                        // consider only owned sides
                        if (bulkData.parallel_owner_rank(
                                dgInfo->currentElement_) !=
                            messager::myProcNo())
                            continue;

                        // if gauss point is exposed (non-overlapping),
                        // then treat as a wall
                        if (dgInfo->gaussPointExposed_)
                        {
                            continue;
                        }

                        // extract current/opposing face/element
                        stk::mesh::Entity currentFace = dgInfo->currentFace_;

                        // local ip, ordinals, etc
                        const label currentGaussPointId =
                            dgInfo->currentGaussPointId_;

                        // pointer to face data
                        const scalar* c_areaVec = stk::mesh::field_data(
                            exposedAreaVecSTKFieldRef, currentFace);
                        const scalar* ncqDot = stk::mesh::field_data(
                            qDotSideSTKFieldRef, currentFace);

                        // area associated to ip
                        scalar c_amag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar c_axj =
                                c_areaVec[currentGaussPointId * SPATIAL_DIM +
                                          j];
                            c_amag += c_axj * c_axj;
                        }
                        c_amag = std::sqrt(c_amag);

                        const scalar heat_flow = ncqDot[currentGaussPointId];
                        if (heat_flow < 0.0) // heat-in
                        {
                            interfaceData.in += heat_flow;
                            interfaceData.in_area += c_amag;
                        }
                        else // heat-out
                        {
                            interfaceData.out += heat_flow;
                            interfaceData.out_area += c_amag;
                        }
                        interfaceData.total_area += c_amag;
                    }
                }
            }

            // Slave
            {
                // get interface side that is sitting in this domain
                const auto* interfaceSideInfoPtr = interf->slaveInfoPtr();

                interfaceDataVector.push_back(HeatBoundaryData());
                HeatBoundaryData& interfaceData = interfaceDataVector.back();

                // set primary traits
                std::strncpy(interfaceData.name,
                             interfaceSideInfoPtr->name().c_str(),
                             sizeof(interfaceData.name));
                interfaceData.name[sizeof(interfaceData.name) - 1] = '\0';

                auto type = interf->type();
                std::strncpy(interfaceData.type,
                             toString(type).c_str(),
                             sizeof(interfaceData.type));
                interfaceData.type[sizeof(interfaceData.type) - 1] = '\0';

                // consider only if qDot is defined on the boundary
                if (!this->qDotRef().sideFieldRef().definedOn(
                        interfaceSideInfoPtr->currentPartVec_))
                {
                    continue;
                }

                // extract vector of dgInfo
                const std::vector<std::vector<dgInfo*>>& dgInfoVec =
                    interfaceSideInfoPtr->dgInfoVec_;

                for (label iSide = 0;
                     iSide < static_cast<label>(dgInfoVec.size());
                     iSide++)
                {
                    const std::vector<dgInfo*>& faceDgInfoVec =
                        dgInfoVec[iSide];

                    // now loop over all the DgInfo objects on this
                    // particular exposed face
                    for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
                    {
                        dgInfo* dgInfo = faceDgInfoVec[k];

                        // consider only owned sides
                        if (bulkData.parallel_owner_rank(
                                dgInfo->currentElement_) !=
                            messager::myProcNo())
                            continue;

                        // if gauss point is exposed (non-overlapping),
                        // then treat as a wall
                        if (dgInfo->gaussPointExposed_)
                        {
                            continue;
                        }

                        // extract current/opposing face/element
                        stk::mesh::Entity currentFace = dgInfo->currentFace_;

                        // local ip, ordinals, etc
                        const label currentGaussPointId =
                            dgInfo->currentGaussPointId_;

                        // pointer to face data
                        const scalar* c_areaVec = stk::mesh::field_data(
                            exposedAreaVecSTKFieldRef, currentFace);
                        const scalar* ncqDot = stk::mesh::field_data(
                            qDotSideSTKFieldRef, currentFace);

                        // area associated to ip
                        scalar c_amag = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar c_axj =
                                c_areaVec[currentGaussPointId * SPATIAL_DIM +
                                          j];
                            c_amag += c_axj * c_axj;
                        }
                        c_amag = std::sqrt(c_amag);

                        const scalar heat_flow = ncqDot[currentGaussPointId];
                        if (heat_flow < 0.0) // heat-in
                        {
                            interfaceData.in += heat_flow;
                            interfaceData.in_area += c_amag;
                        }
                        else // heat-out
                        {
                            interfaceData.out += heat_flow;
                            interfaceData.out_area += c_amag;
                        }
                        interfaceData.total_area += c_amag;
                    }
                }
            }
        }
        else
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(domain->index());

            interfaceDataVector.push_back(HeatBoundaryData());
            HeatBoundaryData& interfaceData = interfaceDataVector.back();

            // set primary traits
            std::strncpy(interfaceData.name,
                         interfaceSideInfoPtr->name().c_str(),
                         sizeof(interfaceData.name));
            interfaceData.name[sizeof(interfaceData.name) - 1] = '\0';

            auto type = interf->type();
            std::strncpy(interfaceData.type,
                         toString(type).c_str(),
                         sizeof(interfaceData.type));
            interfaceData.type[sizeof(interfaceData.type) - 1] = '\0';

            // consider only if qDot is defined on the boundary
            if (!this->qDotRef().sideFieldRef().definedOn(
                    interfaceSideInfoPtr->currentPartVec_))
            {
                continue;
            }

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

                    // consider only owned sides
                    if (bulkData.parallel_owner_rank(dgInfo->currentElement_) !=
                        messager::myProcNo())
                        continue;

                    // if gauss point is exposed (non-overlapping),
                    // then treat as a wall
                    if (dgInfo->gaussPointExposed_)
                    {
                        continue;
                    }

                    // extract current/opposing face/element
                    stk::mesh::Entity currentFace = dgInfo->currentFace_;

                    // local ip, ordinals, etc
                    const label currentGaussPointId =
                        dgInfo->currentGaussPointId_;

                    // pointer to face data
                    const scalar* c_areaVec = stk::mesh::field_data(
                        exposedAreaVecSTKFieldRef, currentFace);
                    const scalar* ncqDot =
                        stk::mesh::field_data(qDotSideSTKFieldRef, currentFace);

                    // area associated to ip
                    scalar c_amag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar c_axj =
                            c_areaVec[currentGaussPointId * SPATIAL_DIM + j];
                        c_amag += c_axj * c_axj;
                    }
                    c_amag = std::sqrt(c_amag);

                    const scalar heat_flow = ncqDot[currentGaussPointId];

                    if (heat_flow < 0.0) // heat-in
                    {
                        interfaceData.in += heat_flow;
                        interfaceData.in_area += c_amag;
                    }
                    else // heat-out
                    {
                        interfaceData.out += heat_flow;
                        interfaceData.out_area += c_amag;
                    }
                    interfaceData.total_area += c_amag;
                }
            }
        }
    }

    if (!noInterfaces)
    {
        std::vector<HeatBoundaryData> globalInterfaceDataVector(
            interfaceDataVector.size());

        MPI_Reduce(interfaceDataVector.data(),
                   globalInterfaceDataVector.data(),
                   interfaceDataVector.size(),
                   MPIHeatBoundaryData,
                   MPIHeatBoundaryData_SUM,
                   0,
                   messager::comm());

        heatInterfaceDataVector_.insert(heatInterfaceDataVector_.end(),
                                        globalInterfaceDataVector.begin(),
                                        globalInterfaceDataVector.end());
    }
}
#endif /* HAS_INTERFACE */

} /* namespace accel */
