// File       : assembler.h
// Created    : Thu Feb 22 2024 10:10:15 (+0100)
// Author     : Fabian Wermelinger
// Description: Abstract base class for equation assembly
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "coefficients.h"
#include "fieldBroker.h"
#include "linearSolverContext.h"
#include "scaling.h"
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
    virtual void preAssemble(const domain*, Context*)
    {
    }

    virtual void postAssemble(const domain*, Context*)
    {
    }

    virtual void assemble(const domain*, Context*) = 0;

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
                const int64_t row =
                    A.getGraph()->localToRow(bulkData.local_id(entity));
                if (row < 0) // node not part of this (subset) graph
                    continue;

                auto rowVals = A.rowVals(row);
                auto rowCols = A.rowCols(row);

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
                scalar* rhs = &b[BLOCKSIZE * row];
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
                const int64_t row =
                    A.getGraph()->localToRow(bulkData.local_id(entity));
                if (row < 0) // inactive node absent from this (subset) graph
                    continue;

                auto rowVals = A.rowVals(row);
                auto rowCols = A.rowCols(row);

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
                auto* diag = A.diag(row);
                for (auto i : dofs)
                {
                    diag[BLOCKSIZE * i + i] = 1.0;
                }

                // zero the rhs
                auto* rhs = &b[BLOCKSIZE * row];
                for (auto i : dofs)
                {
                    rhs[i] = 0.0;
                }
            }
        }

        if (normalise && nodeBuckets.size() > 0)
        {
            linearSolver::normalize(ctx->getCoefficients().getAMatrix(),
                                    ctx->getCoefficients().getBVector());
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
                const int64_t lid =
                    A.getGraph()->localToRow(bulkData.local_id(entity));
                if (lid < 0) // node not part of this (subset) system
                    continue;

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

    // Cross-rank fold+constraint for split conformal/periodic pairs over dof
    // sub-range [r0,r1); T is the NxN column transform (rotation/identity).
    void applyCrossRankFold(
        Context& ctx,
        const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
            pairs,
        const std::array<scalar, N * N>& T,
        int r0,
        int r1,
        const std::function<void(stk::mesh::Entity, scalar*)>& gatherPhiSub);
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

    const auto* graph = A.getGraph();
    const size_t nConnectedNodes = connectedNodes.size();
    const unsigned numRows = nConnectedNodes;
    // rows are graph rows (mesh local-id for full graph, remapped for subset)
    const label nOwnedRows = static_cast<label>(graph->nOwnedNodes());

    STK_ThrowAssert(BLOCKSIZE * numRows <= static_cast<unsigned>(rhs.size()));
    STK_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <=
                    static_cast<unsigned>(lhs.size()));

    scratchIds.resize(numRows);
    matrixColumnIds_.resize(numRows);
    sortPermutation_.resize(numRows);

    if (graph->isLocalColumnOrder())
    {
        for (size_t iNode = 0; iNode < connectedNodes.size(); iNode++)
        {
            stk::mesh::Entity node = connectedNodes[iNode];
            const label row =
                static_cast<label>(graph->localToRow(bulkData.local_id(node)));
            matrixColumnIds_[iNode] = row;
            scratchIds[iNode] = row;
        }
    }
    else
    {
        for (size_t iNode = 0; iNode < connectedNodes.size(); iNode++)
        {
            stk::mesh::Entity node = connectedNodes[iNode];
            // global column: subset id for a subset graph, else STK global id
            matrixColumnIds_[iNode] =
                graph->hasRowRemap()
                    ? static_cast<label>(
                          graph->localToGlobalRow(bulkData.local_id(node)))
                    : static_cast<label>(bulkData.global_id(node));
            scratchIds[iNode] =
                static_cast<label>(graph->localToRow(bulkData.local_id(node)));
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

        // assemble only owned rows (rowID < 0 => node not in this subset graph)
        if (rowID >= 0 && rowID < nOwnedRows)
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

                // column outside this subset graph (curColumnIndex < 0): skip
                if (curColumnIndex < 0)
                    continue;

                while (offset < length && rowCols[offset] != curColumnIndex)
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

template <size_t N>
void assembler<N>::applyCrossRankFold(
    Context& ctx,
    const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>& pairs,
    const std::array<scalar, N * N>& T,
    int r0,
    int r1,
    const std::function<void(stk::mesh::Entity, scalar*)>& gatherPhiSub)
{
    if (!messager::parallel())
        return;

    constexpr label B = static_cast<label>(N);
    constexpr label BB = B * B;
    const label M = r1 - r0; // constrained sub-range size
    const label myProc = messager::myProcNo();

    const stk::mesh::BulkData& bulkData =
        field_broker_->meshRef().bulkDataRef();

    auto& A = ctx.getAMatrix();
    auto& bvec = ctx.getBVector();
    const auto& diagOffsets = A.diagOffsetRef();
    const auto* graph = A.getGraph();
    using GIdx = std::decay_t<decltype(A.localToGlobal(0))>;

    // Slave global row id -> slave node (phi2 reads) and slave column row ->
    // master-entity remap.  The column map is keyed on the graph column index
    // so it covers shared/ghost slave nodes as well as owned rows: a split row
    // can contain interface columns owned by a third rank.
    std::unordered_map<GIdx, stk::mesh::Entity> slaveGToNode2;
    std::unordered_map<int64_t, stk::mesh::EntityId> slaveColToMasterEnt;
    slaveGToNode2.reserve(pairs.size());
    slaveColToMasterEnt.reserve(pairs.size());
    for (const auto& p : pairs)
    {
        const int64_t row2 = graph->localToRow(bulkData.local_id(p.second));
        if (row2 < 0) // pair not in this (subset) graph
            continue;
        slaveColToMasterEnt[row2] = bulkData.identifier(p.first);
        if (bulkData.parallel_owner_rank(p.second) == myProc)
        {
            // The slave row is owned here, so this conversion is valid.
            const GIdx g2 = A.localToGlobal(row2);
            slaveGToNode2[g2] = p.second;
        }
    }

    // ---- Round 1: slave owners send row 2 (global-id columns + NxN blocks +
    // rhs) to the master's owner, which folds T*row2*T^T into row 1 ----
    stk::CommSparse comm1(messager::comm());
    auto pack1 = [&]()
    {
        for (const auto& p : pairs)
        {
            const label o1 = bulkData.parallel_owner_rank(p.first);
            const label o2 = bulkData.parallel_owner_rank(p.second);
            if (o1 == o2 || o2 != myProc)
                continue;
            const label lid2 = static_cast<label>(
                graph->localToRow(bulkData.local_id(p.second)));
            if (lid2 < 0) // pair not in this (subset) graph
                continue;
            auto vals2 = A.rowVals(lid2);
            auto cl2 = graph->rowLocalIndices(lid2);
            auto cg2 = graph->rowGlobalIndices(lid2);
            stk::CommBuffer& buf = comm1.send_buffer(o1);
            buf.pack<stk::mesh::EntityId>(bulkData.identifier(p.first));
            buf.pack<GIdx>(A.localToGlobal(lid2));
            for (label i = 0; i < B; ++i)
                buf.pack<scalar>(bvec[lid2 * B + i]);
            buf.pack<unsigned>(static_cast<unsigned>(cg2.size()));
            for (label c = 0; c < static_cast<label>(cg2.size()); ++c)
            {
                buf.pack<GIdx>(cg2[c]);
                // cl2/cg2 are aligned: resolve the master partner straight
                // from the local column index of this row.
                const auto mit = slaveColToMasterEnt.find(cl2[c]);
                buf.pack<stk::mesh::EntityId>(
                    mit != slaveColToMasterEnt.end() ? mit->second : 0);
                for (label m = 0; m < BB; ++m)
                    buf.pack<scalar>(vals2[c * BB + m]);
            }
        }
    };
    stk::pack_and_communicate(comm1, pack1);

    struct Reply
    {
        GIdx slaveG;
        GIdx masterG;
        std::array<scalar, N * N> diag1;
        std::array<scalar, N> phi1; // master phi sub-range ([0,M))
        label slaveProc;
    };

    std::vector<Reply> replies;
    stk::unpack_communications(comm1,
                               [&](label fromProc)
    {
        stk::CommBuffer& buf = comm1.recv_buffer(fromProc);
        // received column block and its column-transformed copy, row-major BxB;
        // both are fully overwritten per column, so neither needs clearing
        std::array<scalar, N * N> blk;
        std::array<scalar, N * N> bt;
        while (buf.remaining())
        {
            stk::mesh::EntityId masterEnt;
            GIdx slaveG;
            buf.unpack(masterEnt);
            buf.unpack(slaveG);
            std::array<scalar, N> b2;
            for (label i = 0; i < B; ++i)
                buf.unpack(b2[i]);
            unsigned nc;
            buf.unpack(nc);

            const stk::mesh::Entity n1 =
                bulkData.get_entity(stk::topology::NODE_RANK, masterEnt);
            const label lid1 =
                static_cast<label>(graph->localToRow(bulkData.local_id(n1)));
            assert(lid1 >= 0); // master of a subset pair is in the subset
            const GIdx masterG = A.localToGlobal(lid1);
            auto vals1 = A.rowVals(lid1);
            auto cl1 = graph->rowLocalIndices(lid1);
            auto cg1 = graph->rowGlobalIndices(lid1);
            // Row stencils are short (~30-80); a linear scan over contiguous
            // indices beats building a node-based map per received row.
            const auto findPos = [](const auto& idx, const int64_t v) -> label
            {
                for (label c = 0; c < static_cast<label>(idx.size()); ++c)
                    if (static_cast<int64_t>(idx[c]) == v)
                        return c;
                return -1;
            };

            for (unsigned k = 0; k < nc; ++k)
            {
                GIdx cg;
                stk::mesh::EntityId mEnt;
                buf.unpack(cg);
                buf.unpack(mEnt);
                for (label m = 0; m < BB; ++m)
                    buf.unpack(blk[m]);
                // full remap: interface-slave column -> its master column
                label pos = -1;
                bool remappedColumn = false;
                if (mEnt != 0)
                {
                    const stk::mesh::Entity nm =
                        bulkData.get_entity(stk::topology::NODE_RANK, mEnt);
                    if (bulkData.is_valid(nm))
                    {
                        const int64_t rm =
                            graph->localToRow(bulkData.local_id(nm));
                        if (rm >= 0)
                        {
                            pos = findPos(cl1, rm);
                            remappedColumn = pos >= 0;
                        }
                    }
                }
                if (pos < 0) // no remap: keep the original slave column
                    pos = findPos(cg1, static_cast<int64_t>(cg));
                if (pos < 0)
                {
                    // Pair-row graph equalization can leave structural zeros
                    // that are absent from the receiving row. A nonzero
                    // original slave column must always be present; when its
                    // remapped master column is absent we deliberately retain
                    // the slave column and let the explicit constraint enforce
                    // equivalence.
                    const bool nonzero =
                        std::any_of(blk.begin(), blk.end(), [](const scalar v) {
                        return v != 0.0;
                    });
                    if (nonzero)
                        errorMsg("conformal row-merge: nonzero slave column is "
                                 "missing from the master row stencil");
                    continue;
                }
                // Always transform the slave equation rows.  Transform a
                // column as well only when that slave-interface unknown was
                // remapped to its master partner.  Unmapped interior columns
                // remain global-Cartesian unknowns.
                const label off = pos * BB;
                if (remappedColumn)
                {
                    // bt = blk * T^T first, so the fold stays O(B^3)
                    for (label kk = 0; kk < B; ++kk)
                        for (label j = 0; j < B; ++j)
                        {
                            scalar s = 0.0;
                            for (label l = 0; l < B; ++l)
                                s += blk[kk * B + l] * T[j * B + l];
                            bt[kk * B + j] = s;
                        }
                }
                // row1[i in [r0,r1), all j] += (T src)[i,j]
                const scalar* src = remappedColumn ? bt.data() : blk.data();
                for (label i = r0; i < r1; ++i)
                    for (label j = 0; j < B; ++j)
                    {
                        scalar s = 0.0;
                        for (label kk = 0; kk < B; ++kk)
                            s += T[i * B + kk] * src[kk * B + j];
                        vals1[off + i * B + j] += s;
                    }
            }
            // rhs1[i in [r0,r1)] += (T b2)[i]
            for (label i = r0; i < r1; ++i)
            {
                scalar s = 0.0;
                for (label j = 0; j < B; ++j)
                    s += T[i * B + j] * b2[j];
                bvec[lid1 * B + i] += s;
            }

            Reply r;
            r.slaveG = slaveG;
            r.masterG = masterG;
            r.slaveProc = fromProc;
            const label d1 = diagOffsets[lid1] * BB;
            for (label m = 0; m < BB; ++m)
                r.diag1[m] = vals1[d1 + m];
            std::array<scalar, N> phiTmp{};
            gatherPhiSub(n1, phiTmp.data());
            for (label i = 0; i < M; ++i)
                r.phi1[i] = phiTmp[i];
            replies.push_back(r);
        }
    });

    // ---- Round 2: master owners return the folded diagonal block + master
    // phi; slave owners write the scaled constraint row over [r0,r1) ----
    stk::CommSparse comm2(messager::comm());
    auto pack2 = [&]()
    {
        for (const auto& r : replies)
        {
            stk::CommBuffer& buf = comm2.send_buffer(r.slaveProc);
            buf.pack<GIdx>(r.slaveG);
            buf.pack<GIdx>(r.masterG);
            for (label m = 0; m < BB; ++m)
                buf.pack<scalar>(r.diag1[m]);
            for (label i = 0; i < M; ++i)
                buf.pack<scalar>(r.phi1[i]);
        }
    };
    stk::pack_and_communicate(comm2, pack2);

    stk::unpack_communications(comm2,
                               [&](label fromProc)
    {
        stk::CommBuffer& buf = comm2.recv_buffer(fromProc);
        while (buf.remaining())
        {
            GIdx slaveG, masterG;
            buf.unpack(slaveG);
            buf.unpack(masterG);
            std::array<scalar, N * N> D1;
            for (label m = 0; m < BB; ++m)
                buf.unpack(D1[m]);
            std::array<scalar, N> phi1{};
            for (label i = 0; i < M; ++i)
                buf.unpack(phi1[i]);

            const label lid2 = A.globalToLocal(slaveG);
            const stk::mesh::Entity n2 = slaveGToNode2[slaveG];
            auto vals2 = A.rowVals(lid2);
            auto cg2 = graph->rowGlobalIndices(lid2);

            // S = Tsub^T D1sub Tsub (MxM, stored top-left); Tsub = T[r0..,r0..]
            std::array<scalar, N * N> S{};
            for (label a = 0; a < M; ++a)
                for (label bb = 0; bb < M; ++bb)
                {
                    scalar s = 0.0;
                    for (label m = 0; m < M; ++m)
                        for (label n = 0; n < M; ++n)
                            s += T[(r0 + m) * B + (r0 + a)] *
                                 D1[(r0 + m) * B + (r0 + n)] *
                                 T[(r0 + n) * B + (r0 + bb)];
                    S[a * M + bb] = s;
                }

            // zero rows [r0,r1) of row 2 (all columns) and rhs sub-range
            for (label c = 0; c < static_cast<label>(cg2.size()); ++c)
                for (label i = r0; i < r1; ++i)
                    for (label j = 0; j < B; ++j)
                        vals2[c * BB + i * B + j] = 0.0;
            for (label i = r0; i < r1; ++i)
                bvec[lid2 * B + i] = 0.0;

            // diagonal sub-block = S
            const label dOff2 = diagOffsets[lid2] * BB;
            for (label a = 0; a < M; ++a)
                for (label bb = 0; bb < M; ++bb)
                    vals2[dOff2 + (r0 + a) * B + (r0 + bb)] = S[a * M + bb];

            // off-diagonal at master column = -(S Tsub^T)
            label posN1 = -1;
            for (label c = 0; c < static_cast<label>(cg2.size()); ++c)
                if (cg2[c] == masterG)
                {
                    posN1 = c;
                    break;
                }
            if (posN1 >= 0)
            {
                const label o = posN1 * BB;
                for (label a = 0; a < M; ++a)
                    for (label bb = 0; bb < M; ++bb)
                    {
                        scalar s = 0.0;
                        for (label kk = 0; kk < M; ++kk)
                            s += S[a * M + kk] * T[(r0 + bb) * B + (r0 + kk)];
                        vals2[o + (r0 + a) * B + (r0 + bb)] = -s;
                    }
            }

            // rhs = S * (-(phi2 - Tsub^T phi1))
            std::array<scalar, N> phi2{};
            gatherPhiSub(n2, phi2.data());
            std::array<scalar, N> res{};
            for (label a = 0; a < M; ++a)
            {
                scalar tp = 0.0;
                for (label m = 0; m < M; ++m)
                    tp += T[(r0 + m) * B + (r0 + a)] * phi1[m];
                res[a] = -(phi2[a] - tp);
            }
            for (label a = 0; a < M; ++a)
            {
                scalar s = 0.0;
                for (label bb = 0; bb < M; ++bb)
                    s += S[a * M + bb] * res[bb];
                bvec[lid2 * B + (r0 + a)] = s;
            }
        }
    });
}

} /* namespace accel */

#endif // ASSEMBLER_H
