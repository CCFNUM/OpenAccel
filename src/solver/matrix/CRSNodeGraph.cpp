// File       : CRSNodeGraph.cpp
// Created    : Mon Sep 30 2024 10:34:24 (+0200)
// Author     : Fabian Wermelinger
// Description: Node graph generator base class implementation details
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "CRSNodeGraph.h"
#include <cstdint>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <unordered_set>

namespace linearSolver
{

CRSNodeGraph::CRSNodeGraph(const MPI_Comm comm, const GraphLayout layout)
    : comm_(comm), layout_(layout)
{
    resetGraph_();
    setup_();
}

void CRSNodeGraph::buildGraph()
{
    buildGraph_();

    assert(row_ptr_.size() > 1);
    assert(static_cast<Index>(primary_indices_.size()) == row_ptr_.back());
    n_owned_nodes_ = static_cast<Index>(row_ptr_.size()) - 1;
    is_built_ = true;

    computeAuxiliaryData_();
    computePackInfos_();

#ifndef NDEBUG
    // sanity checks:
    assert(is_built_);
    assert(n_owned_nodes_ > 0);
    assert(n_ghost_nodes_ != ~0);
    assert(global_row_offset_ != ~0);
    assert(global_number_nodes_ != ~0);
    assert(global_number_indices_ != ~0ull);

    // a.)
    Index sum_nnz = 0;
    assert(n_owned_nodes_ == static_cast<Index>(row_nnz_owned_.size()));
    assert(n_owned_nodes_ == static_cast<Index>(row_nnz_ghost_.size()));
    for (Index i = 0; i < n_owned_nodes_; i++)
    {
        const Index row_sum_nnz = nnzOwned(i) + nnzGhost(i);
        sum_nnz += row_sum_nnz;
        assert(row_sum_nnz == row_ptr_[i + 1] - row_ptr_[i]);
    }
    assert(sum_nnz == this->nIndices());
    // b.)
    Index sum_received = 0;
    for (const PackInfo& p : pack_infos_)
    {
        sum_received += static_cast<Index>(p.recv_idx.size());
    }
    assert(sum_received <= n_ghost_nodes_);
    if (comm_ != MPI_COMM_NULL)
    {
        for (int r = 0; r < size_; r++)
        {
            if (r == rank_ && sum_received < n_ghost_nodes_)
            {
                std::cout << "WARNING (rank=" << rank_
                          << "): receiving fewer ghost nodes than number of "
                             "ghosts in graph (receiving="
                          << sum_received << "; n_ghosts=" << n_ghost_nodes_
                          << ")" << std::endl;
            }
            MPI_Barrier(comm_);
        }
    }
    // c.)
    assert(n_owned_nodes_ == static_cast<Index>(diagonal_row_offset_.size()));
    for (const Index i : diagonal_row_offset_)
    {
        assert(i >= 0);
    }
#endif /* NDEBUG */
}

std::span<const CRSNodeGraph::Index> CRSNodeGraph::localIndices() const
{
    if (this->isLocalColumnOrder())
    {
        assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
        return std::span<const Index>(primary_indices_); // sorted
    }
    else
    {
        assert(static_cast<Index>(secondary_indices_.size()) ==
               this->nIndices());
        return std::span<const Index>(secondary_indices_); // not sorted
    }
}

std::span<const CRSNodeGraph::Index> CRSNodeGraph::globalIndices() const
{
    if (this->isGlobalColumnOrder())
    {
        assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
        return std::span<const Index>(primary_indices_); // sorted
    }
    else
    {
        assert(static_cast<Index>(secondary_indices_.size()) ==
               this->nIndices());
        return std::span<const Index>(secondary_indices_); // not sorted
    }
}

std::span<const CRSNodeGraph::Index>
CRSNodeGraph::rowLocalIndices(const Index i_row) const
{
    assert(0 <= i_row);
    assert(i_row < this->nOwnedNodes());
    if (this->isLocalColumnOrder())
    {
        assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
        return std::span<const Index>(primary_indices_) // sorted
            .subspan(row_ptr_[i_row], row_ptr_[i_row + 1] - row_ptr_[i_row]);
    }
    else
    {
        assert(static_cast<Index>(secondary_indices_.size()) ==
               this->nIndices());
        return std::span<const Index>(secondary_indices_) // not sorted
            .subspan(row_ptr_[i_row], row_ptr_[i_row + 1] - row_ptr_[i_row]);
    }
}

std::span<const CRSNodeGraph::Index>
CRSNodeGraph::rowGlobalIndices(const Index i_row) const
{
    assert(0 <= i_row);
    assert(i_row < this->nOwnedNodes());
    if (this->isGlobalColumnOrder())
    {
        assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
        return std::span<const Index>(primary_indices_) // sorted
            .subspan(row_ptr_[i_row], row_ptr_[i_row + 1] - row_ptr_[i_row]);
    }
    else
    {
        assert(static_cast<Index>(secondary_indices_.size()) ==
               this->nIndices());
        return std::span<const Index>(secondary_indices_) // not sorted
            .subspan(row_ptr_[i_row], row_ptr_[i_row + 1] - row_ptr_[i_row]);
    }
}

void CRSNodeGraph::getMemoryFootprint(MemoryFootprint& data) const
{
    assert(is_built_);
    std::memset(&data, 0, sizeof(MemoryFootprint));

    if (comm_ != MPI_COMM_NULL)
    {
        // integer array data
        unsigned long long int n_elements = row_ptr_.size();
        n_elements += primary_indices_.size();
        n_elements += secondary_indices_.size();
        n_elements += row_nnz_owned_.size();
        n_elements += row_nnz_ghost_.size();
        n_elements += diagonal_row_offset_.size();
        for (const PackInfo& p : pack_infos_)
        {
            n_elements += p.send_idx.size();
            n_elements += p.recv_idx.size();
        }
        MPI_Reduce(&n_elements,
                   &data.sum_byte,
                   1,
                   MPI_UNSIGNED_LONG_LONG,
                   MPI_SUM,
                   0,
                   comm_);
        data.sum_byte *= sizeof(Index);
    }
}

void CRSNodeGraph::serialize(std::ofstream& out) const
{
    uint64_t v64; // always write 64-bit
    const char* p64 = reinterpret_cast<char*>(&v64);

    // header
    v64 = sizeof(Index); // size of index type
    out.write(p64, sizeof(uint64_t));
    v64 = n_owned_nodes_;
    out.write(p64, sizeof(uint64_t));
    v64 = row_ptr_.size();
    out.write(p64, sizeof(uint64_t));
    v64 = primary_indices_.size();
    out.write(p64, sizeof(uint64_t));

    // data
    for (const Index i : row_ptr_)
    {
        v64 = i;
        out.write(p64, sizeof(uint64_t));
    }
    for (const Index i : primary_indices_)
    {
        v64 = i;
        out.write(p64, sizeof(uint64_t));
    }
    for (const Index i : secondary_indices_)
    {
        v64 = i;
        out.write(p64, sizeof(uint64_t));
    }
}

void CRSNodeGraph::deserialize(std::ifstream& in)
{
    uint64_t v64; // always read 64-bit
    char* p64 = reinterpret_cast<char*>(&v64);
    in.read(p64, sizeof(uint64_t));
    assert(sizeof(Index) == v64);
    in.read(p64, sizeof(uint64_t));
    n_owned_nodes_ = static_cast<Index>(v64);
    in.read(p64, sizeof(uint64_t));
    row_ptr_.resize(v64, 0);
    in.read(p64, sizeof(uint64_t));
    primary_indices_.resize(v64, 0);
    secondary_indices_.resize(v64, 0);

    for (size_t i = 0; i < row_ptr_.size(); i++)
    {
        in.read(p64, sizeof(uint64_t));
        row_ptr_[i] = static_cast<Index>(v64);
    }
    for (size_t i = 0; i < primary_indices_.size(); i++)
    {
        in.read(p64, sizeof(uint64_t));
        primary_indices_[i] = static_cast<Index>(v64);
    }
    for (size_t i = 0; i < secondary_indices_.size(); i++)
    {
        in.read(p64, sizeof(uint64_t));
        secondary_indices_[i] = static_cast<Index>(v64);
    }
}

typename CRSNodeGraph::Index CRSNodeGraph::filterGhostsForOwnerRank_(
    int& owner_rank,
    std::vector<typename CRSNodeGraph::Index>::const_iterator ghost_start,
    const std::vector<typename CRSNodeGraph::Index>::const_iterator ghost_end,
    const int myrank,
    const int size,
    const std::vector<typename CRSNodeGraph::Index>& n_local_nodes) const
{
    const int r0 = (*ghost_start < this->global_row_offset_) ? 0 : myrank + 1;
    const int r1 = (*ghost_start < this->global_row_offset_) ? myrank : size;
    Index j_global =
        (0 == r0) ? 0 : this->global_row_offset_ + n_local_nodes[myrank];
    for (int r = r0; r < r1; r++)
    {
        const Index k_global = j_global + n_local_nodes[r];
        if (j_global <= *ghost_start && *ghost_start < k_global)
        {
            Index n_ghosts = 0;
            while (ghost_start != ghost_end && *ghost_start < k_global)
            {
                ++ghost_start;
                ++n_ghosts;
            }
            owner_rank = r;
            return n_ghosts;
        }
        j_global = k_global;
    }

    // if you end up here something is wrong...
    assert(
        false &&
        "CRSNodeGraph::filterGhostsForOwnerRank_: could not find owner rank");
    return 0;
}

void CRSNodeGraph::computePackInfos_()
{
    if (1 == size_ || comm_ == MPI_COMM_NULL)
    {
        this->pack_infos_.clear();
        this->n_ghost_nodes_ = 0;
        return;
    }

    std::vector<Index> rank_n_nodes(size_);
    MPI_Allgather(&n_owned_nodes_,
                  1,
                  MPIDataType<Index>::type(),
                  rank_n_nodes.data(),
                  1,
                  MPIDataType<Index>::type(),
                  comm_);

    // 1.) determine locally required ghosts
    std::vector<Index> global_ghost;
    std::vector<Index> local_ghost;
    std::unordered_set<Index> ghost_query;
    const auto local_idx = this->localIndices();
    const auto global_idx = this->globalIndices();
    for (Index i = 0; i < n_owned_nodes_; i++)
    {
        for (Index j = row_ptr_[i]; j < row_ptr_[i + 1]; j++)
        {
            const Index j_local = local_idx[j];
            if (j_local >= n_owned_nodes_) // a ghost
            {
                const Index j_global = global_idx[j];
                const auto status = ghost_query.insert(j_global);
                if (status.second) // push this ghost ID once
                {
                    // global index relationship is known, check it is satisfied
                    // a.) total bound
                    assert(j_global >= 0);
                    assert(j_global < global_number_nodes_);
                    // b.) since this is a ghost, it better not be one of my
                    // nodes
                    assert(j_global < global_row_offset_ ||
                           global_row_offset_ + n_owned_nodes_ <= j_global);
                    global_ghost.push_back(j_global);
                    local_ghost.push_back(j_local);
                }
            }
        }
    }
    assert(global_ghost.size() == local_ghost.size());
    // NOTE: [2024-10-08] the fewest (most compact) possible
    // number of ghost nodes is given by `global_ghost.size()`.  The member
    // `this->n_ghost_nodes_` could be set to this value at this point.  It may
    // not always be possible to enforce this most compact representation
    // depending on the buildGraph_ implementation that was executed prior to
    // this method.  To use the most compact number of ghosts the following
    // switch can be enabled
    //
    // determine_n_ghosts_ = true (default is false)
    //
    // The switch is disabled by default for which it is assumed that
    // this->n_ghost_nodes_ has been set during buildGraph_ execution.  If the
    // previously set number of ghosts is smaller than the most compact
    // described above, the assertion below will be raised.
    //
    // Example where number of ghosts must be set externally is the mesh node
    // graph.  The algorithm implemented for the mesh graph involves manual
    // ghosting due to interfaces for example.  Since local and global ID's that
    // are determined in this algorithm depend on the selected nodes, it is not
    // possible to use the most compact number of ghost nodes as it is possible
    // that there are "holes" in the ghost section that would lead to indexing
    // errors in the linear system structure later on. ("Holes" may exist due to
    // nodes that are selected and used for indices enumeration in the mesh node
    // graph algorithm but are not actually ghosts.)
    if (this->determine_n_ghosts_)
    {
        this->n_ghost_nodes_ = static_cast<Index>(global_ghost.size());
    }
    assert(static_cast<Index>(global_ghost.size()) <= this->n_ghost_nodes_);

    // sort ghosts based on global indices
    std::vector<Index> permute(global_ghost.size());
    std::iota(permute.begin(), permute.end(), 0);
    std::sort(permute.begin(),
              permute.end(),
              [&global_ghost](const size_t i, const size_t j)
    { return global_ghost[i] < global_ghost[j]; });
    CRSNodeGraph::permuteInPlace_(global_ghost.data(), permute);
    CRSNodeGraph::permuteInPlace_(local_ghost.data(), permute);

    // 2.) assign recv indices and prepare request buffer
    std::vector<PackInfo> infos(size_);
    std::vector<std::vector<Index>> request_ghosts(size_);
    auto l_ghost = local_ghost.cbegin();
    auto g_ghost = global_ghost.cbegin();
    while (g_ghost != global_ghost.end())
    {
        // filter number of ghosts owned by owner_rank that owns these n
        // consecutive ghost(s) starting at current g_ghost iterator position
        Index owner_rank;
        const Index n_ghosts = filterGhostsForOwnerRank_(owner_rank,
                                                         g_ghost,
                                                         global_ghost.end(),
                                                         rank_,
                                                         size_,
                                                         rank_n_nodes);

        for (Index i = 0; i < n_ghosts; i++)
        {
            assert(g_ghost != global_ghost.end());
            request_ghosts[owner_rank].push_back(*g_ghost++);
            infos[owner_rank].recv_idx.push_back(*l_ghost++);
        }
    }

    // 3.) inform other ranks what messages I expect and exchange indices
    std::vector<Index> n_expected(size_, 0);
    std::vector<Index> n_send_count(size_, 0);
    for (int i = 0; i < size_; i++)
    {
        n_expected[i] = static_cast<Index>(infos[i].recv_idx.size());
    }
    // TODO: [2024-10-06] better way to do this?
    // probing?
    MPI_Alltoall(n_expected.data(),
                 1,
                 MPIDataType<Index>::type(),
                 n_send_count.data(),
                 1,
                 MPIDataType<Index>::type(),
                 comm_);
    for (int i = 0; i < size_; i++)
    {
        infos[i].send_idx.resize(n_send_count[i]);
    }

    std::vector<MPI_Request> requests;
    for (int i = 0; i < size_; i++)
    {
        if (infos[i].send_idx.size() > 0)
        {
            MPI_Request req;
            MPI_Irecv(infos[i].send_idx.data(),
                      infos[i].send_idx.size(),
                      MPIDataType<Index>::type(),
                      i,
                      i,
                      comm_,
                      &req);
            requests.push_back(req);
        }
    }

    for (int i = 0; i < size_; i++)
    {
        if (request_ghosts[i].size() > 0)
        {
            MPI_Send(request_ghosts[i].data(),
                     request_ghosts[i].size(),
                     MPIDataType<Index>::type(),
                     i,
                     rank_,
                     comm_);
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // 4.) convert send indices to local and sanitize pack infos
    this->pack_infos_.clear();
    for (int i = 0; i < size_; i++)
    {
        PackInfo& info = infos[i];
        info.remote_rank = i;
        if (info.send_idx.size() > 0)
        {
            for (size_t k = 0; k < info.send_idx.size(); k++)
            {
                info.send_idx[k] = this->globalToLocalIndex(info.send_idx[k]);
            }
        }
        if (info.recv_idx.size() > 0 || info.send_idx.size() > 0)
        {
            this->pack_infos_.push_back(info);
        }
    }
}

void CRSNodeGraph::computeDiagonalIndices_()
{
    assert(is_built_);
    if (!diagonal_row_offset_.empty())
    {
        assert(n_owned_nodes_ ==
               static_cast<Index>(diagonal_row_offset_.size()));
        return;
    }
    diagonal_row_offset_.resize(n_owned_nodes_, -1);
    const auto local_idx = this->localIndices();
    for (Index i = 0; i < n_owned_nodes_; i++)
    {
        for (Index J = row_ptr_[i]; J < row_ptr_[i + 1]; J++)
        {
            if (i == local_idx[J])
            {
                diagonal_row_offset_[i] = J - row_ptr_[i];
                break;
            }
        }
    }
}

void CRSNodeGraph::computeAuxiliaryData_()
{
    assert(is_built_);

    computeGlobalRowOffset_();
    computeGlobalGraphSize_();
    // following may depend on global data
    computeRowNonZeros_();
    computeDiagonalIndices_();
}

void CRSNodeGraph::computeRowNonZeros_()
{
    assert(is_built_);
    if (!row_nnz_owned_.empty() || !row_nnz_ghost_.empty())
    {
        assert(n_owned_nodes_ == static_cast<Index>(row_nnz_owned_.size()));
        assert(row_nnz_owned_.size() == row_nnz_ghost_.size());
        return;
    }
    row_nnz_owned_.resize(n_owned_nodes_, 0);
    row_nnz_ghost_.resize(n_owned_nodes_, 0);
    const auto local_idx = this->localIndices();
    for (Index i = 0; i < n_owned_nodes_; i++)
    {
        for (Index J = row_ptr_[i]; J < row_ptr_[i + 1]; J++)
        {
            if (local_idx[J] < n_owned_nodes_)
            {
                ++row_nnz_owned_[i];
            }
            else
            {
                ++row_nnz_ghost_[i];
            }
        }
    }
}

void CRSNodeGraph::computeGlobalRowOffset_()
{
    assert(comm_ != MPI_COMM_NULL);
    if (global_row_offset_ == ~0)
    {
        assert(n_owned_nodes_ > 0);
        const unsigned long long n_rows = n_owned_nodes_;
        unsigned long long my_offset = 0;
        MPI_Exscan(
            &n_rows, &my_offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm_);
        global_row_offset_ = static_cast<Index>(my_offset);
    }
}

void CRSNodeGraph::computeGlobalGraphSize_()
{
    assert(is_built_);
    assert(comm_ != MPI_COMM_NULL);
    if (global_number_nodes_ == ~0 || global_number_indices_ == ~0ull)
    {
        const unsigned long long local[2] = {
            static_cast<unsigned long long>(this->nOwnedNodes()),
            static_cast<unsigned long long>(this->nIndices())};
        unsigned long long global[2] = {0};
        MPI_Allreduce(local, global, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm_);
        global_number_nodes_ = static_cast<Index>(global[0]);
        global_number_indices_ = global[1];
    }
}

void CRSNodeGraph::sortPrimaryIndices_()
{
    assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
    const bool has_secondary =
        (static_cast<Index>(secondary_indices_.size()) == this->nIndices())
            ? true
            : false;

#ifndef NDEBUG
    if (has_secondary)
    {
        assert(static_cast<Index>(secondary_indices_.size()) ==
               this->nIndices());
    }
#endif /* NDEBUG */

    std::vector<Index> buffer;
    std::vector<Index> permute;
    for (Index i_row = 0; i_row < n_owned_nodes_; ++i_row)
    {
        std::span<Index> primary_idx =
            std::span<Index>(primary_indices_)
                .subspan(row_ptr_[i_row],
                         row_ptr_[i_row + 1] - row_ptr_[i_row]);

        permute.resize(primary_idx.size());
        std::iota(permute.begin(), permute.end(), 0);
        std::sort(permute.begin(),
                  permute.end(),
                  [primary_idx](const size_t i, const size_t j)
        { return primary_idx[i] < primary_idx[j]; });
        CRSNodeGraph::permuteCopy_(primary_idx.data(), buffer, permute);

        if (has_secondary)
        {
            std::span<Index> secondary_idx =
                std::span<Index>(secondary_indices_)
                    .subspan(row_ptr_[i_row],
                             row_ptr_[i_row + 1] - row_ptr_[i_row]);
            CRSNodeGraph::permuteCopy_(secondary_idx.data(), buffer, permute);
        }
    }
}

void CRSNodeGraph::setup_()
{
    rank_ = -1;
    size_ = -1;
    if (comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);
    }
}

void CRSNodeGraph::resetGraph_()
{
    is_built_ = false;
    n_owned_nodes_ = 0;
    n_ghost_nodes_ = ~0;

    determine_n_ghosts_ = false;
    global_row_offset_ = ~0;
    global_number_nodes_ = ~0;
    global_number_indices_ = ~0ull;

    row_ptr_.clear();
    primary_indices_.clear();
    secondary_indices_.clear();

    row_nnz_owned_.clear();
    row_nnz_ghost_.clear();
    diagonal_row_offset_.clear();

    pack_infos_.clear();
}

} /* namespace linearSolver */
