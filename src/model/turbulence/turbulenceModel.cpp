// File : turbulenceModel.cpp
// Created : Mon Mar 25 2024 16:48:19 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "turbulenceModel.h"
#include "initialConditions.h"
#include "messager.h"
#include "simulation.h"

namespace accel
{

turbulenceModel::turbulenceModel(realm* realm) : model(realm)
{
    // create field instances
    mutRef();
    muEffRef();

    // wall-function
    yPlusRef();
    uPlusRef();
    yStarRef();
    uStarRef();
    uTauRef();
    duPlusdyPlusRef();

    // create heat transfer related fields in case required
    for (const std::shared_ptr<domain>& domain :
         realmPtr_->simulationRef().domainVector())
    {
        if (domain->heatTransfer_.option_ != heatTransferOption::none)
        {
            // create turb thermal field instances
            lambdatRef();
            lambdaEffRef();

            // thermal wall-function
            TPlusRef();

            break;
        }
    }
}

void turbulenceModel::clipMinDistToWall(const std::shared_ptr<domain> domain)
{
    STKScalarField* minDistanceToWallSTKFieldPtr = yMinRef().stkFieldPtr();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // extract fields required
    const STKScalarField* wallNormalDistanceSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::wall_normal_distance_ID);

    // Process for no-slip wall parts
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

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

        // mapping from ip to nodes for this
        // ordinal; face perspective (use with
        // face_node_relations)
        const label* faceIpNodeMap = meFC->ipNodeMap();

        const stk::mesh::Bucket::size_type length = sideBucket.size();

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // get face
            stk::mesh::Entity side = sideBucket[k];
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);

            label numSideNodes = bulkData.num_nodes(side);

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label localFaceNode = faceIpNodeMap[ip];

                // right is on the face
                stk::mesh::Entity node = sideNodeRels[localFaceNode];

                // assemble to nodal quantities
                scalar* minD =
                    stk::mesh::field_data(*minDistanceToWallSTKFieldPtr, node);
                (*minD) = std::max(*minD, wallNormalDistanceBip[ip]);
            }
        }
    }

    // sync in case of parallel
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData,
                                          {minDistanceToWallSTKFieldPtr});
    }
}

void turbulenceModel::updateEffectiveDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    // Get fields
    auto& muSTKFieldRef = muRef().stkFieldRef();
    auto& muEffSTKFieldRef = muEffRef().stkFieldRef();

    stk::mesh::MetaData& metaData = meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    // define some common selectors
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    auto& mutSTKFieldRef = mutRef().stkFieldRef();

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        const scalar* mu = stk::mesh::field_data(muSTKFieldRef, b);
        const scalar* mut = stk::mesh::field_data(mutSTKFieldRef, b);
        scalar* muEff = stk::mesh::field_data(muEffSTKFieldRef, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            muEff[k] = mu[k] + mut[k];
        }
    }
}

void turbulenceModel::updateTurbulentThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (domain->heatTransfer_.option_ == heatTransferOption::none)
    {
        return;
    }

    // Get fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& mutSTKFieldRef = mutRef().stkFieldRef();
    auto& lambdatSTKFieldRef = lambdatRef().stkFieldRef();

    stk::mesh::MetaData& metaData = meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    scalar Prt = domain->turbulence_.turbulentFluxClosureForHeatTransfer_
                     .turbulentPrandtlNumber_;

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

        const scalar* cp = stk::mesh::field_data(cpSTKFieldRef, b);
        const scalar* mut = stk::mesh::field_data(mutSTKFieldRef, b);
        scalar* lambdat = stk::mesh::field_data(lambdatSTKFieldRef, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            lambdat[k] = cp[k] * mut[k] / Prt;
        }
    }
}

void turbulenceModel::updateEffectiveThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (domain->heatTransfer_.option_ == heatTransferOption::none)
    {
        return;
    }

    // common
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    // Get effective thermal conductivity
    auto& lambdaEffSTKFieldRef = lambdaEffRef().stkFieldRef();

    // update from laminar and turbulent lambdas on all nodes
    {
        const auto& lambdaSTKFieldRef = lambdaRef().stkFieldRef();
        const auto& lambdatSTKFieldRef = lambdatRef().stkFieldRef();

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

            const scalar* lambdab = stk::mesh::field_data(lambdaSTKFieldRef, b);
            const scalar* lambdatb =
                stk::mesh::field_data(lambdatSTKFieldRef, b);
            scalar* lambdaEffb = stk::mesh::field_data(lambdaEffSTKFieldRef, b);

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                lambdaEffb[k] = lambdab[k] + lambdatb[k];
            }
        }
    }
}

const stk::mesh::PartVector
turbulenceModel::collectNoSlipWallParts_(const std::shared_ptr<domain> domain)
{
    stk::mesh::PartVector parts;

    // no-slip boundary walls
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        switch (domain->zonePtr()->boundaryRef(iBoundary).type())
        {
            case boundaryPhysicalType::wall:
                {
                    boundaryConditionType bcType =
                        URef()
                            .boundaryConditionRef(domain->index(), iBoundary)
                            .type();

                    switch (bcType)
                    {
                        case boundaryConditionType::noSlip:
                            {
                                for (const auto part :
                                     domain->zonePtr()
                                         ->boundaryRef(iBoundary)
                                         .parts())
                                {
                                    parts.push_back(part);
                                }
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

    return parts;
}

// Initialize

void turbulenceModel::initializeTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::initializeTurbulentKineticEnergy(domain);

    updateTurbulentKineticEnergySideFields_(domain);
}

void turbulenceModel::initializeTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::initializeTurbulentEddyFrequency(domain);

    updateTurbulentEddyFrequencySideFields_(domain);
}

void turbulenceModel::initializeTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::initializeTurbulentDissipationRate(domain);

    updateTurbulentDissipationRateSideFields_(domain);
}

// Update

void turbulenceModel::updateTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::updateTurbulentKineticEnergy(domain);

    updateTurbulentKineticEnergySideFields_(domain);
}

void turbulenceModel::updateTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::updateTurbulentEddyFrequency(domain);

    updateTurbulentEddyFrequencySideFields_(domain);
}

void turbulenceModel::updateTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::updateTurbulentDissipationRate(domain);

    updateTurbulentDissipationRateSideFields_(domain);
}

void turbulenceModel::updateTurbulentKineticEnergySideFields_(
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
            kRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::intensityAndLengthScale:
                            {
                                updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndLengthScale_(
                                    domain, boundary);
                            }
                            break;

                        case boundaryConditionType::
                            intensityAndEddyViscosityRatio:
                            {
                                updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
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

void turbulenceModel::updateTurbulentEddyFrequencySideFields_(
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
            omegaRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::intensityAndLengthScale:
                            {
                                updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndLengthScale_(
                                    domain, boundary);
                            }
                            break;

                        case boundaryConditionType::
                            intensityAndEddyViscosityRatio:
                            {
                                updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
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

void turbulenceModel::updateTurbulentDissipationRateSideFields_(
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
            epsilonRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::intensityAndLengthScale:
                            {
                                updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndLengthScale_(
                                    domain, boundary);
                            }
                            break;

                        case boundaryConditionType::
                            intensityAndEddyViscosityRatio:
                            {
                                updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
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

void turbulenceModel::
    updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc = kRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& I = bc.data<1>("I");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues = kRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (I.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar IValue = *I.value();

                // Get fields
                const STKScalarField* USTKFieldPtr = this->URef().stkFieldPtr();
                STKScalarField* nodeSideKSTKFieldPtr =
                    this->kRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

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
                    scalar* kb = stk::mesh::field_data(*nodeSideKSTKFieldPtr,
                                                       nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        scalar Umag = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            Umag += Ub[SPATIAL_DIM * iNode + i] *
                                    Ub[SPATIAL_DIM * iNode + i];
                        }
                        Umag = std::sqrt(Umag);

                        kb[iNode] =
                            3.0 / 2.0 *
                            std::pow(IValue * std::max(0.01, Umag), 2.0);
                    }
                }

                // interpolate to node-side field
                kRef().sideFieldRef().interpolate(kRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  this->kRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    kRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

void turbulenceModel::
    updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc =
        omegaRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& l = bc.data<1>("l");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues = omegaRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (l.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar lValue = *l.value();

                // Get fields
                const STKScalarField* kSTKFieldPtr = this->kRef().stkFieldPtr();
                STKScalarField* nodeSideOmegaSTKFieldPtr =
                    this->omegaRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);

                const scalar Cmu = this->Cmu();

                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* kb =
                        stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
                    scalar* omegab = stk::mesh::field_data(
                        *nodeSideOmegaSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        omegab[iNode] = std::sqrt(kb[iNode]) / (Cmu * lValue);
                    }
                }

                // interpolate to node-side field
                omegaRef().sideFieldRef().interpolate(
                    omegaRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->omegaRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    omegaRef().correctBoundaryNodes(domain->index(),
                                                    boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

void turbulenceModel::
    updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc =
        epsilonRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& l = bc.data<1>("l");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues =
        epsilonRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (l.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar lValue = *l.value();

                // Get fields
                const STKScalarField* kSTKFieldPtr = this->kRef().stkFieldPtr();
                STKScalarField* nodeSideEpsilonSTKFieldPtr =
                    this->epsilonRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

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
                    const scalar* kb =
                        stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
                    scalar* epsilonb = stk::mesh::field_data(
                        *nodeSideEpsilonSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        epsilonb[iNode] =
                            std::pow(kb[iNode], 3.0 / 2.0) / lValue;
                    }
                }

                // interpolate to node-side field
                epsilonRef().sideFieldRef().interpolate(
                    epsilonRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->epsilonRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    epsilonRef().correctBoundaryNodes(domain->index(),
                                                      boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

void turbulenceModel::
    updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc = kRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& I = bc.data<1>("I");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues = kRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (I.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar IValue = *I.value();

                // Get fields
                const STKScalarField* USTKFieldPtr = this->URef().stkFieldPtr();
                STKScalarField* nodeSideKSTKFieldPtr =
                    this->kRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

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
                    scalar* kb = stk::mesh::field_data(*nodeSideKSTKFieldPtr,
                                                       nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        scalar Umag = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            Umag += Ub[SPATIAL_DIM * iNode + i] *
                                    Ub[SPATIAL_DIM * iNode + i];
                        }
                        Umag = std::sqrt(Umag);

                        kb[iNode] =
                            3.0 / 2.0 *
                            std::pow(IValue * std::max(0.01, Umag), 2.0);
                    }
                }

                // interpolate to node-side field
                kRef().sideFieldRef().interpolate(kRef().nodeSideFieldRef(),
                                                  domain->index(),
                                                  boundary->index(),
                                                  this->kRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    kRef().correctBoundaryNodes(domain->index(),
                                                boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

void turbulenceModel::
    updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc =
        omegaRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& r = bc.data<1>("r");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues = omegaRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (r.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar rValue = *r.value();

                // Get fields
                const STKScalarField* kSTKFieldPtr = this->kRef().stkFieldPtr();
                const STKScalarField* muSTKFieldPtr =
                    this->muRef().stkFieldPtr();
                const STKScalarField* rhoSTKFieldPtr =
                    this->rhoRef().stkFieldPtr();
                STKScalarField* nodeSideOmegaSTKFieldPtr =
                    this->omegaRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

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
                    const scalar* kb =
                        stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
                    const scalar* mub =
                        stk::mesh::field_data(*muSTKFieldPtr, nodeBucket);
                    const scalar* rhob =
                        stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
                    scalar* omegab = stk::mesh::field_data(
                        *nodeSideOmegaSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        omegab[iNode] =
                            kb[iNode] / (mub[iNode] / rhob[iNode] * rValue);
                    }
                }

                // interpolate to node-side field
                omegaRef().sideFieldRef().interpolate(
                    omegaRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->omegaRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    omegaRef().correctBoundaryNodes(domain->index(),
                                                    boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

void turbulenceModel::
    updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    // get boundary condition data
    auto& bc =
        epsilonRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& r = bc.data<1>("r");

    // boolean to translate boundary values to field
    bool correctedBoundaryNodeValues =
        epsilonRef().correctedBoundaryNodeValues();

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (r.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                const scalar rValue = *r.value();

                // Get fields
                const STKScalarField* kSTKFieldPtr = this->kRef().stkFieldPtr();
                const STKScalarField* muSTKFieldPtr =
                    this->muRef().stkFieldPtr();
                const STKScalarField* rhoSTKFieldPtr =
                    this->rhoRef().stkFieldPtr();
                STKScalarField* nodeSideEpsilonSTKFieldPtr =
                    this->epsilonRef().nodeSideFieldRef().stkFieldPtr();

                // get interior parts the domain is defined on
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                // define some common selectors; select owned nodes
                stk::mesh::Selector selUniversalNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK,
                                         selUniversalNodes);

                const scalar Cmu = this->Cmu();

                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    // field chunks in bucket
                    const scalar* kb =
                        stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
                    const scalar* mub =
                        stk::mesh::field_data(*muSTKFieldPtr, nodeBucket);
                    const scalar* rhob =
                        stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
                    scalar* epsilonb = stk::mesh::field_data(
                        *nodeSideEpsilonSTKFieldPtr, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        epsilonb[iNode] = Cmu * std::pow(kb[iNode], 2.0) /
                                          (mub[iNode] / rhob[iNode] * rValue);
                    }
                }

                // interpolate to node-side field
                epsilonRef().sideFieldRef().interpolate(
                    epsilonRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->epsilonRef().isShifted());

                // translate to field
                if (correctedBoundaryNodeValues)
                {
                    epsilonRef().correctBoundaryNodes(domain->index(),
                                                      boundary->index());
                }
            }
            break;

        default:
            errorMsg("Not implemented");
    }
}

} /* namespace accel */
