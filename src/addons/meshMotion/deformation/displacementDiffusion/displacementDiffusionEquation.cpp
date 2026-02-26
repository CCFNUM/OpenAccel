// File       : displacementDiffusionEquation.cpp
// Created    : Tue Nov 26 2024
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "displacementDiffusionEquation.h"
#include "initialConditions.h"
#include "meshMotion.h"

namespace accel
{

displacementDiffusionEquation::displacementDiffusionEquation(realm* realm)
    : equation("Displacement Diffusion"), displacementDiffusionModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    // Collect equation names
    std::vector<std::string> eqNames;

#if SPATIAL_DIM == 2
    eqNames.push_back("D-x");
    eqNames.push_back("D-y");
#elif SPATIAL_DIM == 3
    eqNames.push_back("D-x");
    eqNames.push_back("D-y");
    eqNames.push_back("D-z");
#endif

    // set
    this->setEquationName(eqNames);

    // disable residual plot of the equation
    plot_res_ = false;
}

bool displacementDiffusionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void displacementDiffusionEquation::setup()
{
    // setup of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(setupDisplacement);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler
    assembler_->setup(&this->DRef(),
                      diffusion,
                      domainVector_,
                      // anonymous function to compute Gamma for displacement
                      // diffusion equation:
                      [this](const domain* domain, STKScalarField& Gamma)
    {
        const auto& mesh = this->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        // define selector for domain
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        const BucketVector& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

        using stk::mesh::field_data;

        if (domain->zonePtr()
                ->deformationRef()
                .displacementDiffusion()
                .meshStiffnessSpecification_ ==
            meshStiffnessSpecificationType::value)
        {
            scalar meshStiffnessValue = domain->zonePtr()
                                            ->deformationRef()
                                            .displacementDiffusion()
                                            .meshStiffnessValue_;
            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                scalar* Gammab = field_data(Gamma, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    Gammab[i] = meshStiffnessValue;
                }
            }
        }
        else if (domain->zonePtr()
                     ->deformationRef()
                     .displacementDiffusion()
                     .meshStiffnessSpecification_ ==
                 meshStiffnessSpecificationType::increaseNearSmallVolumes)
        {
            scalar modelExponent = domain->zonePtr()
                                       ->deformationRef()
                                       .displacementDiffusion()
                                       .modelExponent_;

            const auto& dualNodalVolumeSTKFieldRef =
                *metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                            mesh::dual_nodal_volume_ID);

            // calculate reference volume: mean control volume
            scalar volRef = 0.0;
            label volCount = 0;
            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                const scalar* volb =
                    field_data(dualNodalVolumeSTKFieldRef, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    volRef += volb[i];
                    volCount++;
                }
            }

            messager::sumReduce(volRef);
            messager::sumReduce(volCount);

            // average
            volRef /= static_cast<scalar>(volCount);

            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                scalar* Gammab = field_data(Gamma, nodeBucket);
                const scalar* volb =
                    field_data(dualNodalVolumeSTKFieldRef, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    scalar stiffness =
                        std::pow(volRef / volb[i], modelExponent);

                    // store and clip
                    Gammab[i] = std::max(std::min(stiffness, 1e15), 1e-15);
                }
            }
        }
        else if (domain->zonePtr()
                     ->deformationRef()
                     .displacementDiffusion()
                     .meshStiffnessSpecification_ ==
                 meshStiffnessSpecificationType::increaseNearBoundaries)
        {
            scalar modelExponent = domain->zonePtr()
                                       ->deformationRef()
                                       .displacementDiffusion()
                                       .modelExponent_;

            scalar referenceLengthScale = domain->zonePtr()
                                              ->deformationRef()
                                              .displacementDiffusion()
                                              .referenceLengthScale_;

            const auto& yMinSTKFieldRef = this->yMinRef().stkFieldRef();

            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                scalar* Gammab = field_data(Gamma, nodeBucket);
                const scalar* yMinb = field_data(yMinSTKFieldRef, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    scalar stiffness = std::pow(referenceLengthScale / yMinb[i],
                                                modelExponent);

                    // store and clip
                    Gammab[i] = std::max(std::min(stiffness, 1e15), 1e-15);
                }
            }
        }
        else if (domain->zonePtr()
                     ->deformationRef()
                     .displacementDiffusion()
                     .meshStiffnessSpecification_ ==
                 meshStiffnessSpecificationType::blendedDistanceAndSmallVolumes)
        {
            const auto& dict =
                domain->zonePtr()->deformationRef().displacementDiffusion();

            const scalar A = dict.blendedVolumeWeight_;
            const scalar B = dict.blendedDistanceWeight_;
            const scalar C_vol = dict.blendedVolumeExponent_;
            const scalar C_dis = dict.blendedDistanceExponent_;

            const auto& dualNodalVolumeSTKFieldRef =
                *metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                            mesh::dual_nodal_volume_ID);

            const auto& yMinSTKFieldRef = this->yMinRef().stkFieldRef();

            // compute reference values from local domain mesh:
            // V_ref = mean control volume
            // V_min = minimum control volume (for d_wall)
            scalar volRef = 0.0;
            scalar volMin = std::numeric_limits<scalar>::max();
            scalar volTotal = 0.0;
            label volCount = 0;

            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                const scalar* volb =
                    field_data(dualNodalVolumeSTKFieldRef, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    volTotal += volb[i];
                    volCount++;
                    volMin = std::min(volMin, volb[i]);
                }
            }

            messager::sumReduce(volTotal);
            messager::sumReduce(volCount);
            messager::minReduce(volMin);

            // V_ref = mean control volume in domain
            volRef = volTotal / static_cast<scalar>(volCount);

            // L_ref = 0.5 * (volume of domain)^(1/3)
            const scalar domainVolume = domain->zonePtr()->stats().volume_;
            const scalar L_ref = 0.5 * std::cbrt(domainVolume);

            // d_wall = 10.0 * (minimum control volume in domain)^(1/3)
            const scalar d_wall = 10.0 * std::cbrt(volMin);

            for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
            {
                const Bucket& nodeBucket = *nodeBuckets[ib];
                const Bucket::size_type nNodesPerBucket = nodeBucket.size();

                scalar* Gammab = field_data(Gamma, nodeBucket);
                const scalar* volb =
                    field_data(dualNodalVolumeSTKFieldRef, nodeBucket);
                const scalar* yMinb = field_data(yMinSTKFieldRef, nodeBucket);

                for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
                {
                    // Gamma = A*(V_ref/V)^C_vol
                    // + B*(L_ref/max(d, d_wall))^C_dis
                    const scalar volTerm =
                        A * std::pow(volRef / volb[i], C_vol);

                    const scalar effectiveDistance = std::max(yMinb[i], d_wall);
                    const scalar disTerm =
                        B * std::pow(L_ref / effectiveDistance, C_dis);

                    scalar stiffness = volTerm + disTerm;

                    // store and clip
                    Gammab[i] = std::max(std::min(stiffness, 1e15), 1e-15);
                }
            }
        }
    });

    // setup linear solver
    // FIXME: Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), fieldBroker::meshRef());

    equation::isCreated_ = true;
}

void displacementDiffusionEquation::initialize()
{
    FOREACH_DOMAIN(initializeDisplacement);

    // post initialization
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // update scale
    this->DRef().updateScale();

    equation::isInitialized_ = true;
}

void displacementDiffusionEquation::postInitialize()
{
}

void displacementDiffusionEquation::preSolve()
{
    FOREACH_DOMAIN(updateDisplacement);
}

void displacementDiffusionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    linearSystem::simulationRef().getProfiler().push("linear_system_assembly");

    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get());

    // fix system in all the stationary parts
    assembler_->fix(this->collectStationaryParts(), {}, ctx.get());

    // fix system on all dirichlet boundaries
    assembler_->fix(this->collectDirichletBoundaryParts_(), {}, ctx.get());

    linearSystem::simulationRef().getProfiler().pop();

    // solve linear system
    linearSystem::solve();

    // correction
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            this->DRef().stkFieldRef());

        // synchronize
        this->DRef().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // update gradient
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // update scale
    this->DRef().updateScale();
}

void displacementDiffusionEquation::postSolve()
{
}

stk::mesh::PartVector
displacementDiffusionEquation::collectDirichletBoundaryParts_()
{
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
#ifdef HAS_INTERFACE
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                for (auto part : interf->interfaceSideInfoPtr(domain->index())
                                     ->currentPartVec_)
                {
                    incPartVec.push_back(part);
                }
            }
        }
#endif /* HAS_INTERFACE */

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            switch (bcType)
            {
                case boundaryConditionType::specifiedValue:
                case boundaryConditionType::periodicDisplacement:
                case boundaryConditionType::rigidBodySolution:
                    {
                        const auto& boundaryRef =
                            domain->zonePtr()->boundaryRef(iBoundary);

                        for (auto part : boundaryRef.parts())
                        {
                            incPartVec.push_back(part);
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    return incPartVec;
}

void displacementDiffusionEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateDisplacementPrevTimeField);
}

void displacementDiffusionEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    const scalar D_scale_inv = 1.0 / (this->DRef().scale() + ::accel::SMALL);
    for (int i = 0; i < SPATIAL_DIM; i++)
    {
        residual_scales[i] = D_scale_inv;
    }
}

} /* namespace accel */
