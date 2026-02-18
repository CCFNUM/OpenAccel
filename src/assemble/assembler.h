// File : assembler.h
// Created : Thu Feb 22 2024 10:10:15 (+0100)
// Author : Fabian Wermelinger
// Description: Abstract base class for equation assembly
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "controls.h"
#include "fieldBroker.h"
#include "linearSolverContext.h"
#include "types.h"

namespace accel
{

class domain;

template <size_t N>
class assembler
{
public:
    static constexpr int BLOCKSIZE = N;
    using Context = ::linearSolver::Context<N>;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;

    assembler() = delete;

    assembler(fieldBroker* field_broker) : field_broker_(field_broker)
    {
    }

    virtual ~assembler()
    {
    }

    assembler(const assembler& c) = default;
    assembler& operator=(const assembler& c) = default;
    assembler(assembler&& c) = default;
    assembler& operator=(assembler&& c) = default;

    // public API
    virtual void assemble(const domain* domain, Context* ctx)
    {
        preAssemble_(domain, ctx);
        assemble_(domain, ctx);
        postAssemble_(domain, ctx);
    }

    // overloaded for parts input
    virtual void zero(stk::mesh::PartVector incPartVec,
                      stk::mesh::PartVector excPartVec,
                      Context* ctx,
                      std::vector<label> ignoreDofs = {})
    {
        using Bucket = stk::mesh::Bucket;
        using BucketVec = stk::mesh::BucketVector;

        // select all locally owned nodes for this domain
        const auto& mesh = field_broker_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        stk::mesh::Selector selOwnedNodes = metaData.locally_owned_part() &
                                            stk::mesh::selectUnion(incPartVec) &
                                            !stk::mesh::selectUnion(excPartVec);

        const BucketVec& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        // create target dof's first
        std::vector<label> dofs;
        {
            std::unordered_set<label> ignore(ignoreDofs.begin(),
                                             ignoreDofs.end());
            dofs.reserve(BLOCKSIZE - ignoreDofs.size());

            for (int i = 0; i < BLOCKSIZE; ++i)
            {
                if (ignore.find(i) == ignore.end())
                {
                    dofs.push_back(i);
                }
            }
        }

        // loop over local nodes and override lhs and rhs
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& bucket = *nodeBuckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            for (Bucket::size_type i_ent = 0; i_ent < n_entities; ++i_ent)
            {
                const stk::mesh::Entity entity = bucket[i_ent];
                const auto lid = bulkData.local_id(entity);

                auto rowVals = A.rowVals(lid);
                auto rowCols = A.rowCols(lid);

                // zero the whole row
                for (label icol = 0; icol < static_cast<label>(rowCols.size());
                     icol++)
                {
                    for (auto i : dofs)
                    {
                        for (label j = 0; j < BLOCKSIZE; j++)
                        {
                            rowVals[BLOCKSIZE * BLOCKSIZE * icol +
                                    BLOCKSIZE * i + j] = 0.0;
                        }
                    }
                }

                // zero the rhs
                scalar* rhs = &b[BLOCKSIZE * lid];
                for (auto i : dofs)
                {
                    rhs[i] = 0.0;
                }
            }
        }
    }

    virtual void fix(stk::mesh::PartVector incPartVec,
                     stk::mesh::PartVector excPartVec,
                     Context* ctx,
                     std::vector<label> ignoreDofs = {},
                     bool normalise = false)
    {
        using Bucket = stk::mesh::Bucket;
        using BucketVec = stk::mesh::BucketVector;

        // select all locally owned nodes for this domain
        const auto& mesh = field_broker_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        stk::mesh::Selector selOwnedNodes = metaData.locally_owned_part() &
                                            stk::mesh::selectUnion(incPartVec) &
                                            !stk::mesh::selectUnion(excPartVec);

        const BucketVec& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        // create target dof's first
        std::vector<label> dofs;
        {
            std::unordered_set<label> ignore(ignoreDofs.begin(),
                                             ignoreDofs.end());
            dofs.reserve(BLOCKSIZE - ignoreDofs.size());

            for (int i = 0; i < BLOCKSIZE; ++i)
            {
                if (ignore.find(i) == ignore.end())
                {
                    dofs.push_back(i);
                }
            }
        }

        // loop over local nodes and override lhs and rhs
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& bucket = *nodeBuckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            for (Bucket::size_type i_ent = 0; i_ent < n_entities; ++i_ent)
            {
                const stk::mesh::Entity entity = bucket[i_ent];
                const auto lid = bulkData.local_id(entity);

                auto rowVals = A.rowVals(lid);
                auto rowCols = A.rowCols(lid);

                // zero the whole row
                for (auto icol = 0; icol < rowCols.size(); icol++)
                {
                    for (auto i : dofs)
                    {
                        for (auto j = 0; j < BLOCKSIZE; j++)
                        {
                            rowVals[BLOCKSIZE * BLOCKSIZE * icol +
                                    BLOCKSIZE * i + j] = 0.0;
                        }
                    }
                }

                // set diagonal of diagonal block to 1
                auto* diag = A.diag(lid);
                for (auto i : dofs)
                {
                    diag[BLOCKSIZE * i + i] = 1.0;
                }

                // zero the rhs
                auto* rhs = &b[BLOCKSIZE * lid];
                for (auto i : dofs)
                {
                    rhs[i] = 0.0;
                }
            }
        }

        if (normalise && nodeBuckets.size() > 0)
        {
            ctx->getCoefficients().normalize();
        }
    }

    virtual void
    assembleDiagonalDominance(const domain* domain, Context* ctx, label dof = 0)
    {
        using Bucket = stk::mesh::Bucket;
        using BucketVec = stk::mesh::BucketVector;

        // select all locally owned nodes for this domain
        const auto& mesh = field_broker_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        const zone* zonePtr = domain->zonePtr();

        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(zonePtr->interiorParts());

        const BucketVec& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        const auto& diagOffset = A.diagOffsetRef();

        // loop over local nodes and override lhs and rhs
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& bucket = *nodeBuckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            for (Bucket::size_type i = 0; i < n_entities; ++i)
            {
                const stk::mesh::Entity entity = bucket[i];
                const auto lid = bulkData.local_id(entity);

                auto rowVals = A.rowVals(lid);
                auto rowCols = A.rowCols(lid);

                // get reference to diagonal of equation
                scalar& dia = A.dofDiag(lid, dof);

                // sum-up all off-diagonal entries of coefficients
                scalar sum = 0.0;
                for (label icol = 0; icol < rowCols.size(); icol++)
                {
                    if (diagOffset[lid] == icol)
                        continue;

                    scalar& offDiag = rowVals[BLOCKSIZE * BLOCKSIZE * icol +
                                              BLOCKSIZE * dof + dof];
                    if (offDiag > 0.0)
                    {
                        offDiag = 0.0;
                    }
                    sum -= offDiag;
                }

                // set diagonal
                dia = std::max(dia, sum);
            }
        }
    }

protected:
    fieldBroker* field_broker_;

    // work arrays used for for `applyCoeff_`
    std::vector<typename Context::Index> sortPermutation_;
    std::vector<typename Context::Index> matrixColumnIds_;

    // polymorphic protected interface
    virtual void preAssemble_(const domain*, Context*)
    {
    }

    virtual void postAssemble_(const domain*, Context*)
    {
    }

    virtual void assemble_(const domain*, Context*) = 0;

    // helper methods
    void applyCoeff_(
        Matrix& A,
        Vector& b,
        const std::vector<stk::mesh::Entity>& sym_meshobj,
        std::vector<label>& scratchIds,
        std::vector<scalar>& /* scratchVals */, // FIXME: Unused parameter
        const std::vector<scalar>& rhs,
        const std::vector<scalar>& lhs,
        bool deductUnfound = true);
};

template <size_t N>
void assembler<N>::applyCoeff_(
    Matrix& A,
    Vector& b,
    const std::vector<stk::mesh::Entity>& connectedNodes,
    std::vector<label>& scratchIds,
    std::vector<scalar>& /* scratchVals */, // FIXME: Unused parameter
    const std::vector<scalar>& rhs,
    const std::vector<scalar>& lhs,
    bool deductUnfound)
{
    const auto& mesh = field_broker_->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const size_t nConnectedNodes = connectedNodes.size();
    const unsigned numRows = nConnectedNodes;
    const label nOwnedRows = mesh.nNodes();

    STK_ThrowAssert(BLOCKSIZE * numRows <= static_cast<unsigned>(rhs.size()));
    STK_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <=
                    static_cast<unsigned>(lhs.size()));

    scratchIds.resize(numRows);
    matrixColumnIds_.resize(numRows);
    sortPermutation_.resize(numRows);

    if (A.getGraph()->isLocalColumnOrder())
    {
        for (size_t iNode = 0; iNode < connectedNodes.size(); iNode++)
        {
            stk::mesh::Entity node = connectedNodes[iNode];
            matrixColumnIds_[iNode] = bulkData.local_id(node);
            scratchIds[iNode] = bulkData.local_id(node);
        }
    }
    else
    {
        for (size_t iNode = 0; iNode < connectedNodes.size(); iNode++)
        {
            stk::mesh::Entity node = connectedNodes[iNode];
            matrixColumnIds_[iNode] = bulkData.global_id(node);
            scratchIds[iNode] = bulkData.local_id(node);
        }
    }

    // Sort
    std::iota(sortPermutation_.begin(), sortPermutation_.end(), 0);
    std::sort(sortPermutation_.begin(),
              sortPermutation_.end(),
              [this](size_t i1, size_t i2)
    { return this->matrixColumnIds_[i1] < this->matrixColumnIds_[i2]; });

    for (unsigned r = 0; r < numRows; r++)
    {
        const label cur_perm_index = sortPermutation_[r];
        const label rowID = scratchIds[cur_perm_index];

        const scalar* const lhsCurrent =
            &lhs[BLOCKSIZE * BLOCKSIZE * cur_perm_index * numRows];
        const scalar* const rhsCurrent = &rhs[BLOCKSIZE * cur_perm_index];

#ifndef NDEBUG
        for (label k = 0; k < BLOCKSIZE; k++)
        {
            STK_ThrowAssertMsg(std::isfinite(rhsCurrent[k]), "Invalid rhs");
        }
#endif /* NDEBUG */

        // assemble only owned rows
        if (rowID < nOwnedRows)
        {
            auto rowVals = A.rowVals(rowID);
            const auto rowCols = A.rowCols(rowID);

            const label length = rowCols.size();
            const label numCols = nConnectedNodes;

            label offset = 0;
            label offsetOld = 0;
            for (label j = 0; j < numCols; j++)
            {
                const label permIndex = sortPermutation_[j];
                const label curColumnIndex = matrixColumnIds_[permIndex];

                while (rowCols[offset] != curColumnIndex && offset < length)
                {
                    ++offset;
                }

                if (offset < length)
                {
                    for (label m = 0; m < BLOCKSIZE; m++)
                    {
                        for (label n = 0; n < BLOCKSIZE; n++)
                        {
                            const label column_idx =
                                BLOCKSIZE * BLOCKSIZE * offset + BLOCKSIZE * m +
                                n;
                            const label local_idx = BLOCKSIZE * numRows * m +
                                                    BLOCKSIZE * permIndex + n;
#ifndef NDEBUG
                            STK_ThrowAssertMsg(
                                std::isfinite(lhsCurrent[local_idx]),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            rowVals[column_idx] += lhsCurrent[local_idx];
                        }
                    }

                    // store this offset
                    offsetOld = offset;
                }
                else if (deductUnfound)
                {
                    // reduced stencil (remove reduced nodes contribution to
                    // diagonal)
                    scalar* diag = A.diag(rowID);
                    for (label m = 0; m < BLOCKSIZE; m++)
                    {
                        for (label n = 0; n < BLOCKSIZE; n++)
                        {
                            const label local_idx = BLOCKSIZE * numRows * m +
                                                    BLOCKSIZE * permIndex + n;
#ifndef NDEBUG
                            STK_ThrowAssertMsg(
                                std::isfinite(lhsCurrent[local_idx]),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            diag[BLOCKSIZE * m + n] += lhsCurrent[local_idx];
                        }
                    }

                    // reset offset
                    offset = offsetOld;
                }
            }
            for (label k = 0; k < BLOCKSIZE; k++)
            {
                b[BLOCKSIZE * rowID + k] += rhsCurrent[k];
            }
        }
    }
}

} /* namespace accel */

#endif // ASSEMBLER_H
