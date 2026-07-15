// File       : turbulentDissipationRateAssembler.h
// Created    : Thu Feb 22 2025 13:38:51 (+0100)
// Author     : Achraf Nagihi
// Description: Assembler for the turbulent dissipation rate equation in
//              k-epsilon model
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTDISSIPATIONRATEASSEMBLER_H
#define TURBULENTDISSIPATIONRATEASSEMBLER_H

#include "kEpsilonModel.h"
#include "phiAssembler.h"

namespace accel
{

class turbulentDissipationRateAssembler : public phiAssembler<1>
{
private:
    kEpsilonModel* model_;

public:
    using Base = phiAssembler<1>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

public:
    turbulentDissipationRateAssembler(kEpsilonModel* model)
        : Base(model), model_(model)
    {
    }

    void postAssemble(const domain* domain, Context* ctx) override
    {
        Base::postAssemble(domain, ctx);
        assembleBoundaryRelaxation_(domain, ctx, 0.75);

        stk::mesh::BulkData& bulkData = model_->meshRef().bulkDataRef();
        stk::mesh::MetaData& metaData = model_->meshRef().metaDataRef();

        auto& A = ctx->getAMatrix();
        auto& b = ctx->getBVector();

        // collect no-slip walls

        stk::mesh::PartVector parts;

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
                const int64_t id =
                    A.getGraph()->localToRow(bulkData.local_id(nodeBucket[k]));
                if (id < 0) // node not part of this (subset) system
                    continue;

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
};

} /* namespace accel */

#endif // TURBULENTDISSIPATIONRATEASSEMBLER_H
