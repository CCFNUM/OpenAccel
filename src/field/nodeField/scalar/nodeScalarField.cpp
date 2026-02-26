// File       : nodeScalarField.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "nodeScalarField.h"
#include "boundary.h"
#include "domain.h"
#include "messager.h"
#include "realm.h"
#include "simulation.h"
#include "zone.h"

namespace accel
{

nodeScalarField::nodeScalarField(realm* realmPtr,
                                 std::string name,
                                 unsigned numberOfStates,
                                 bool prevIter,
                                 bool highResolution,
                                 bool computeGradient,
                                 bool correctedBoundaryNodeValues)
    : nodeField<1, SPATIAL_DIM>(realmPtr->meshPtr(),
                                name,
                                numberOfStates,
                                prevIter,
                                highResolution,
                                computeGradient,
                                correctedBoundaryNodeValues),
      realmPtr_(realmPtr)
{
}

nodeScalarField::nodeScalarField(realm* realmPtr, STKScalarField* stkField_ptr)
    : nodeField<1, SPATIAL_DIM>(realmPtr->meshPtr(), stkField_ptr),
      realmPtr_(realmPtr)
{
}

// Methods

void nodeScalarField::correctGradientField_(label iZone)
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
                        *metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            mesh::assembled_symm_area_ID);

                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundaryRef.parts());
                    const auto& sideNodeBuckets = bulkData.get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);

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

                            scalar dot = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                dot += grad[j] * aarea[j] / amag;
                            }

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                grad[j] -= dot * aarea[j] / amag;
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

simulation* nodeScalarField::simulationPtr()
{
    return realmRef().simulationPtr();
}

const simulation* nodeScalarField::simulationPtr() const
{
    return realmRef().simulationPtr();
}

simulation& nodeScalarField::simulationRef()
{
    return realmRef().simulationRef();
}

const simulation& nodeScalarField::simulationRef() const
{
    return realmRef().simulationRef();
}

realm* nodeScalarField::realmPtr()
{
    return realmPtr_;
}

const realm* nodeScalarField::realmPtr() const
{
    return realmPtr_;
}

realm& nodeScalarField::realmRef()
{
    return *realmPtr_;
}

const realm& nodeScalarField::realmRef() const
{
    return *realmPtr_;
}

} // namespace accel
