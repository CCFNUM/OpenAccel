// File       : UtopiaFSIGraph.hpp
// Created    : Wed Oct 16 2024 10:50:46 (+0200)
// Author     : Fabian Wermelinger
// Description: Utopia FSI node graph implementation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// no include guard: must only be included in UtopiaFSIGraph.h

#include <cassert>
#include <set>

namespace accel
{

namespace fsi
{

template <size_t N>
void UtopiaFSIGraph<N>::buildGraph(const UMatrix& umat)
{
    this->resetGraph_();

    n_ghosts_ = 0;
    o_local_block_map_.clear();
    o_local_to_global_ = nullptr;
    if (1 == this->commSize())
    {
        d_petsc_ = umat.raw_type();
    }
    else
    {
        MatMPIAIJGetSeqAIJ(
            umat.raw_type(), &d_petsc_, &o_petsc_, &o_local_to_global_);
        MatGetGhosts(umat.raw_type(), &n_ghosts_, nullptr);
    }
    CRSNodeGraph::buildGraph();

    d_petsc_ = nullptr;
    o_petsc_ = nullptr;
}

template <size_t N>
void UtopiaFSIGraph<N>::buildGraph_()
{
    assert(d_petsc_);

    PetscBool ok;
    PetscInt d_rows;
    PetscInt o_rows;
    const PetscInt* d_row_ptr = nullptr;
    const PetscInt* d_indices = nullptr;
    const PetscInt* o_row_ptr = nullptr;
    const PetscInt* o_indices = nullptr;

    MatGetRowIJ(d_petsc_,
                0,
                PETSC_FALSE,
                PETSC_FALSE,
                &d_rows,
                &d_row_ptr,
                &d_indices,
                &ok);
    assert(ok);
    assert(d_rows % BLOCKSIZE == 0);
    assert(d_row_ptr);
    assert(d_indices);

    std::vector<Index> o_l2g_block_id; // local to global block idx map
    if (o_local_to_global_)
    {
        assert(o_petsc_);
        MatGetRowIJ(o_petsc_,
                    0,
                    PETSC_FALSE,
                    PETSC_FALSE,
                    &o_rows,
                    &o_row_ptr,
                    &o_indices,
                    &ok);
        assert(ok);
        assert(o_rows % BLOCKSIZE == 0);
        assert(o_row_ptr);
        assert(o_indices);
        assert(d_rows == o_rows);

        assert(n_ghosts_ > 0);
        o_local_block_map_.resize(n_ghosts_);
        const Index local_shift = o_local_to_global_[0] / BLOCKSIZE;
        Index prev_block = 0;
        Index next_block_id = 0;
        o_l2g_block_id.push_back(local_shift);
        for (PetscInt i = 0; i < n_ghosts_; i++)
        {
            const Index global_id = o_local_to_global_[i] / BLOCKSIZE;
            assert(global_id >= local_shift); // must be sorted
            const Index local_id = global_id - local_shift;
            if (prev_block != local_id)
            {
                ++next_block_id;
                prev_block = local_id;
                o_l2g_block_id.push_back(global_id);
            }
            o_local_block_map_[i] = next_block_id;
        }
    }

    // scan flat rows and determined block nnz
    this->n_owned_nodes_ = d_rows / BLOCKSIZE;
    this->computeGlobalRowOffset_();

    row_ptr_.resize(this->n_owned_nodes_ + 1);
    row_ptr_[0] = 0;
    assert(this->getLayout() &
           linearSolver::GraphLayout::ColumnIndexOrder__Local);
    for (Index i = 0; i < this->n_owned_nodes_; i++)
    {
        std::set<Index> d_local_idx;
        std::set<Index> o_local_idx;

        // scan diagonal part
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            const PetscInt start = d_row_ptr[BLOCKSIZE * i + k];
            const PetscInt end = d_row_ptr[BLOCKSIZE * i + k + 1];
            for (PetscInt j = start; j < end; j++)
            {
                d_local_idx.insert(d_indices[j] / BLOCKSIZE);
            }
        }

        // scan off-diagonal part
        if (o_local_to_global_)
        {
            for (Index k = 0; k < BLOCKSIZE; k++)
            {
                const PetscInt start = o_row_ptr[BLOCKSIZE * i + k];
                const PetscInt end = o_row_ptr[BLOCKSIZE * i + k + 1];
                for (PetscInt j = start; j < end; j++)
                {
                    o_local_idx.insert(o_local_block_map_[o_indices[j]]);
                }
            }
        }

        std::vector<Index> local_idx;  // sorted
        std::vector<Index> global_idx; // unsorted
        const Index row_block_nnz = d_local_idx.size() + o_local_idx.size();
        local_idx.reserve(row_block_nnz);
        global_idx.reserve(row_block_nnz);
        for (const Index j : d_local_idx)
        {
            local_idx.push_back(j);
            global_idx.push_back(j + this->global_row_offset_);
        }
        for (const Index j : o_local_idx)
        {
            local_idx.push_back(j + this->n_owned_nodes_);
            assert(j < static_cast<Index>(o_l2g_block_id.size()));
            global_idx.push_back(o_l2g_block_id[j]);
        }
        assert(static_cast<Index>(local_idx.size()) == row_block_nnz);
        assert(static_cast<Index>(global_idx.size()) == row_block_nnz);

        row_ptr_[i + 1] = row_ptr_[i] + row_block_nnz;
        primary_indices_.insert(
            primary_indices_.end(), local_idx.begin(), local_idx.end());
        secondary_indices_.insert(
            secondary_indices_.end(), global_idx.begin(), global_idx.end());
    }

    this->determine_n_ghosts_ = true;

    // clang-format off
    MatRestoreRowIJ(d_petsc_, 0, PETSC_FALSE, PETSC_FALSE, &d_rows, &d_row_ptr, &d_indices, &ok);
    assert(ok);
    if (o_local_to_global_) {
        MatRestoreRowIJ(o_petsc_, 0, PETSC_FALSE, PETSC_FALSE, &o_rows, &o_row_ptr, &o_indices, &ok);
        assert(ok);
    }
    // clang-format on
}

} /* namespace fsi */
} /* namespace accel */
