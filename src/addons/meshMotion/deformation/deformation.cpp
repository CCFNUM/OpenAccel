// File       : deformation.cpp
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "deformation.h"
#include "displacementDiffusionEquation.h"
#include "meshMotion.h"
#include "realm.h"

namespace accel
{

deformation::deformation(meshMotion* meshMotionPtr)
    : meshMotionPtr_(meshMotionPtr)
{
    // register the domain for displacement diffusion equation
    for (auto domain :
         meshMotionPtr_->realmRef().simulationRef().domainVector())
    {
        // add domains on which the displacement diffusion will be solved for
        if (domain->hasEquation(equationID::displacementDiffusion))
        {
            if (displacementDiffusionEquation_ == nullptr)
            {
                displacementDiffusionEquation_ =
                    std::make_unique<displacementDiffusionEquation>(
                        &meshMotionPtr_->realmRef());
            }

            // register domain for equation
            displacementDiffusionEquation_->addDomain(domain);
        }
    }

#ifndef NDEBUG
    // Sanity check
    for (auto domain :
         meshMotionPtr_->realmRef().simulationRef().domainVector())
    {
        if (domain->zonePtr()->meshDeforming())
        {
            if (domain->zonePtr()->deformationRef().specification() ==
                meshDeformationSpecificationType::regionsOfMotionSpecified)
            {
                if (messager::master())
                {
                    std::cout << "updating total displacement field due to "
                                 "deformation "
                                 "in a region of motion specified "
                                 "in zone: "
                              << domain->zonePtr()->name() << std::endl;
                }
            }
            else if (domain->zonePtr()->deformationRef().specification() ==
                     meshDeformationSpecificationType::inherent)
            {
                if (messager::master())
                {
                    std::cout << "updating total displacement field due to "
                                 "deformation "
                                 "in a region of inherent motion "
                                 "in zone: "
                              << domain->zonePtr()->name() << std::endl;
                }
            }
            else
            {
                errorMsg("Must not reach here. deformation is instantiated but "
                         "no valid mesh deformation specification");
            }
        }
    }
#endif
}

void deformation::setup()
{
    if (displacementDiffusionEquation_)
    {
        displacementDiffusionEquation_->setup();
    }
}

void deformation::initialize()
{
    if (displacementDiffusionEquation_)
    {
        displacementDiffusionEquation_->initialize();
    }
}

void deformation::update()
{
    // solve for the displacement field due to wall motion
    if (displacementDiffusionEquation_)
    {
        label maxIterations =
            meshMotionPtr_->meshRef()
                .controlsRef()
                .solverRef()
                .solverControl_.advancedOptions_.equationControls_.meshMotion_
                .maxSmoothingIterations_;

        // re-use control iterations
        label& iter = meshMotionPtr_->meshRef().controlsRef().iter;

        // store the iter to reset afterwards
        label iterOrig = iter;

        // loop and solver displacement diffusion
        for (iter = 1; iter <= maxIterations; iter++)
        {
            displacementDiffusionEquation_->preSolve();
            displacementDiffusionEquation_->solve();
            displacementDiffusionEquation_->postSolve();

            // check if converged
            if (displacementDiffusionEquation_->isConverged())
            {
                break;
            }
        }

        // set back control iter not to screw up the other equations
        iter = iterOrig;
    }

    // update total displacement: inherent motion is already by now updated
    for (auto domain :
         meshMotionPtr_->realmRef().simulationRef().domainVector())
    {
        if (domain->zonePtr()->meshDeforming())
        {
            auto& mesh = meshMotionPtr_->meshRef();
            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
            stk::mesh::MetaData& metaData = mesh.metaDataRef();

            // Get fields
            auto& DtSTKFieldRef = meshMotionPtr_->DtRef().stkFieldRef();
            const auto& DtSTKFieldRefOld =
                meshMotionPtr_->DtRef().prevTimeRef().stkFieldRef();
            const auto& DSTKFieldRef = meshMotionPtr_->DRef().stkFieldRef();

            // define some common selectors; select all nodes
            stk::mesh::Selector selUniversalNodes =
                metaData.universal_part() &
                stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

            label relDisp = domain->zonePtr()
                                ->deformationRef()
                                .displacementRelativeToPreviousMesh();

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
                const scalar* Db =
                    stk::mesh::field_data(DSTKFieldRef, nodeBucket);
                const scalar* DtbOld =
                    stk::mesh::field_data(DtSTKFieldRefOld, nodeBucket);
                scalar* Dtb = stk::mesh::field_data(DtSTKFieldRef, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        Dtb[SPATIAL_DIM * iNode + i] =
                            relDisp * DtbOld[SPATIAL_DIM * iNode + i] +
                            Db[SPATIAL_DIM * iNode + i];
                    }
                }
            }
        }
    }
}

} /* namespace accel */
