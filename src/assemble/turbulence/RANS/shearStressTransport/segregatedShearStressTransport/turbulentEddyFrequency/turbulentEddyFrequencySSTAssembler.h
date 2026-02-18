// File : turbulentEddyFrequencySSTAssembler.h
// Created : Thu Feb 22 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Assembler for the turbulent eddy frequency equation in SST model
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTEDDYFREQUENCYSSTASSEMBLER_H
#define TURBULENTEDDYFREQUENCYSSTASSEMBLER_H

#include "shearStressTransportModel.h"
#include "turbulentEddyFrequencyAssembler.h"

namespace accel
{

class turbulentEddyFrequencySSTAssembler
    : public turbulentEddyFrequencyAssembler
{
private:
    shearStressTransportModel* model_;

public:
    turbulentEddyFrequencySSTAssembler(shearStressTransportModel* model)
        : turbulentEddyFrequencyAssembler(model), model_(model)
    {
    }

protected:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;
    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;
    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    // Boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    void postAssemble_(const domain* domain, Context* ctx) override
    {
        turbulentEddyFrequencyAssembler::postAssemble_(domain, ctx);

        stk::mesh::BulkData& bulkData = model_->meshRef().bulkDataRef();
        stk::mesh::MetaData& metaData = model_->meshRef().metaDataRef();

        auto& A = ctx->getAMatrix();
        auto& b = ctx->getBVector();

        // collect no-slip walls

        stk::mesh::PartVector parts;

#ifdef HAS_INTERFACE
        // fluid-solid interface side
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                // get interface side that is sitting in this domain
                const auto* interfaceSideInfoPtr =
                    interf->interfaceSideInfoPtr(domain->index());

                for (const auto part : interfaceSideInfoPtr->currentPartVec_)
                {
                    parts.push_back(part);
                }
            }
        }
#endif /* HAS_INTERFACE */

        // no-slip boundary walls
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            switch (domain->zonePtr()->boundaryRef(iBoundary).type())
            {
                case boundaryPhysicalType::wall:
                    {
                        boundaryConditionType bcType =
                            model_->URef()
                                .boundaryConditionRef(domain->index(),
                                                      iBoundary)
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

        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(parts);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const stk::mesh::Bucket::size_type length = nodeBucket.size();

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                stk::mesh::EntityId id = bulkData.local_id(nodeBucket[k]);

                auto rowVals = A.rowVals(id);

                // save diag
                scalar d = *A.diag(id);

                // zero row
                for (label i = 0; i < rowVals.size(); i++)
                {
                    rowVals[i] = 0.0;
                }

                // restore diag
                *(A.diag(id)) = d;

                // zero rhs
                b[id] = 0.0;
            }
        }
    }
};

} /* namespace accel */

#endif // TURBULENTEDDYFREQUENCYSSTASSEMBLER_H
