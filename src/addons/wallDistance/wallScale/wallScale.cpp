// File       : wallScale.cpp
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "wallScale.h"
#include "realm.h"
#include "wallDistance.h"
#include "wallScaleDiffusionEquation.h"

namespace accel
{

wallScale::wallScale(wallDistance* wallDistancePtr)
    : wallDistancePtr_(wallDistancePtr)
{
    // register the domain for wall scale diffusion equation
    for (auto domain :
         wallDistancePtr_->realmRef().simulationRef().domainVector())
    {
        if (domain->isWallDistanceRequired())
        {
            if (wallScaleDiffusionEquation_ == nullptr)
            {
                wallScaleDiffusionEquation_ =
                    std::make_unique<wallScaleDiffusionEquation>(
                        &wallDistancePtr_->realmRef());
            }

            // register domain for equation
            wallScaleDiffusionEquation_->addDomain(domain);
        }
    }
}

void wallScale::setup()
{
    if (wallScaleDiffusionEquation_)
    {
        wallScaleDiffusionEquation_->setup();
    }
}

void wallScale::initialize()
{
    if (wallScaleDiffusionEquation_)
    {
        wallScaleDiffusionEquation_->initialize();
    }
}

void wallScale::update()
{
    // solve for the displacement field due to wall motion
    if (wallScaleDiffusionEquation_)
    {
        if (wallDistancePtr_->meshRef().controlsRef().isTransient())
        {
            label maxIterations = wallDistancePtr_->meshRef()
                                      .controlsRef()
                                      .solverRef()
                                      .solverControl_.basicSettings_
                                      .convergenceControl_.maxIterations_;

            // re-use control iterations
            label& iter = wallDistancePtr_->meshRef().controlsRef().iter;

            // store the iter to reset afterwards
            label iterOrig = iter;

            // loop and solver displacement diffusion
            for (iter = 1; iter <= maxIterations; iter++)
            {
                wallScaleDiffusionEquation_->preSolve();
                wallScaleDiffusionEquation_->solve();
                wallScaleDiffusionEquation_->postSolve();

                // check if converged
                if (wallScaleDiffusionEquation_->isConverged())
                {
                    break;
                }
            }

            // set back control iter not to screw up the other equations
            iter = iterOrig;
        }
        else
        {
            if (!wallScaleDiffusionEquation_->isConverged())
            {
                wallScaleDiffusionEquation_->preSolve();
                wallScaleDiffusionEquation_->solve();
                wallScaleDiffusionEquation_->postSolve();
            }
        }

        // update minimum distance field
        for (auto domain :
             wallDistancePtr_->realmRef().simulationRef().domainVector())
        {
            if (domain->isWallDistanceRequired())
            {
                auto& mesh = wallDistancePtr_->meshRef();
                stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
                stk::mesh::MetaData& metaData = mesh.metaDataRef();

                // required fields
                STKScalarField* yScaleSTKFieldPtr =
                    wallDistancePtr_->yScaleRef().stkFieldPtr();
                STKScalarField* gradYScaleSTKFieldPtr =
                    wallDistancePtr_->yScaleRef().gradRef().stkFieldPtr();

                STKScalarField* minDistanceToWallSTKFieldPtr =
                    wallDistancePtr_->yMinRef().stkFieldPtr();

                // define some common selectors
                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

                stk::mesh::BucketVector const& node_buckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         node_buckets.begin();
                     ib != node_buckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& b = **ib;
                    const stk::mesh::Bucket::size_type length = b.size();

                    scalar* yScaleb =
                        stk::mesh::field_data(*yScaleSTKFieldPtr, b);
                    scalar* yMinDistb =
                        stk::mesh::field_data(*minDistanceToWallSTKFieldPtr, b);
                    scalar* gradYScaleb =
                        stk::mesh::field_data(*gradYScaleSTKFieldPtr, b);

                    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
                    {
                        const scalar yScale = yScaleb[k];
                        const scalar* gradYScale =
                            &gradYScaleb[SPATIAL_DIM * k];

                        scalar dpdxsq = 0;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            dpdxsq += gradYScale[i] * gradYScale[i];
                        }

                        yMinDistb[k] =
                            fmax(-std::sqrt(dpdxsq) +
                                     std::sqrt(dpdxsq + 2.0 * yScale),
                                 0.0);
                    }
                }
            }
        }
    }
}

} // namespace accel
