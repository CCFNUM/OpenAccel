// File : nodeVectorField.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "nodeVectorField.h"
#include "boundary.h"
#include "domain.h"
#include "realm.h"
#include "simulation.h"
#include "zone.h"

namespace accel
{

nodeVectorField::nodeVectorField(realm* realmPtr,
                                 std::string name,
                                 unsigned numberOfStates,
                                 bool prevIter,
                                 bool highResolution,
                                 bool computeGradient,
                                 bool correctedBoundaryNodeValues)
    : nodeField<SPATIAL_DIM, SPATIAL_DIM * SPATIAL_DIM>(
          realmPtr->meshPtr(),
          name,
          numberOfStates,
          prevIter,
          highResolution,
          computeGradient,
          correctedBoundaryNodeValues),
      realmPtr_(realmPtr)
{
}

nodeVectorField::nodeVectorField(realm* realmPtr, STKScalarField* stkField_ptr)
    : nodeField<SPATIAL_DIM, SPATIAL_DIM * SPATIAL_DIM>(realmPtr->meshPtr(),
                                                        stkField_ptr),
      realmPtr_(realmPtr)
{
}

// Methods

void nodeVectorField::correctGradientField_(label iZone)
{
    // ensure zone is active
    assert(this->isZoneSet(iZone));

    const auto& bulkData = this->meshRef().bulkDataRef();
    const auto& metaData = this->meshRef().metaDataRef();

    for (label iBoundary = 0;
         iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
         iBoundary++)
    {
        const auto& boundaryRef =
            this->meshRef().zoneRef(iZone).boundaryRef(iBoundary);
        boundaryPhysicalType physicalType = boundaryRef.type();

        switch (physicalType)
        {
            case boundaryPhysicalType::symmetry:
                {
                    // get fields
                    auto* gradPhiSTKFieldPtr = this->gradRef().stkFieldPtr();

                    const auto& assembledSymmSTKFieldRef =
                        *metaData.get_field<scalar>(
                            stk::topology::NODE_RANK,
                            mesh::assembled_symm_area_ID);

                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundaryRef.parts());
                    const auto& sideNodeBuckets = bulkData.get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);

                    // Velocity correction using full tensor projection
                    std::vector<scalar> normalGrad(SPATIAL_DIM);
                    std::vector<scalar> gradNormalVec(SPATIAL_DIM);

                    for (const stk::mesh::Bucket* bucket : sideNodeBuckets)
                    {
                        const stk::mesh::Bucket& sideNodeBucket = *bucket;
                        const auto nSideNodesPerBucket = sideNodeBucket.size();

                        for (size_t iNode = 0; iNode < nSideNodesPerBucket;
                             ++iNode)
                        {
                            stk::mesh::Entity node = sideNodeBucket[iNode];
                            scalar* grad = stk::mesh::field_data(
                                *gradPhiSTKFieldPtr, node);
                            const scalar* aarea = stk::mesh::field_data(
                                assembledSymmSTKFieldRef, node);

                            scalar asq = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj = aarea[j];
                                asq += axj * axj;
                            }
                            const scalar amag = std::sqrt(asq);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                normalGrad[i] = 0;
                                gradNormalVec[i] = 0;
                            }

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar nij = aarea[j] / amag;
                                    normalGrad[i] +=
                                        grad[i * SPATIAL_DIM + j] * nij;
                                    gradNormalVec[j] +=
                                        grad[i * SPATIAL_DIM + j] * aarea[i] /
                                        amag;
                                }
                            }

                            scalar normalGradNormalVec = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                normalGradNormalVec +=
                                    gradNormalVec[j] * aarea[j] / amag;
                            }

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar nij = aarea[j] / amag;
                                    const scalar nji = aarea[i] / amag;
                                    grad[i * SPATIAL_DIM + j] -=
                                        normalGrad[i] * nij;
                                    grad[i * SPATIAL_DIM + j] -=
                                        nji * gradNormalVec[j];
                                    grad[i * SPATIAL_DIM + j] +=
                                        2.0 * nji * nij * normalGradNormalVec;
                                }
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

// Access

simulation* nodeVectorField::simulationPtr()
{
    return realmRef().simulationPtr();
}

const simulation* nodeVectorField::simulationPtr() const
{
    return realmRef().simulationPtr();
}

simulation& nodeVectorField::simulationRef()
{
    return realmRef().simulationRef();
}

const simulation& nodeVectorField::simulationRef() const
{
    return realmRef().simulationRef();
}

realm* nodeVectorField::realmPtr()
{
    return realmPtr_;
}

const realm* nodeVectorField::realmPtr() const
{
    return realmPtr_;
}

realm& nodeVectorField::realmRef()
{
    return *realmPtr_;
}

const realm& nodeVectorField::realmRef() const
{
    return *realmPtr_;
}

} // namespace accel
