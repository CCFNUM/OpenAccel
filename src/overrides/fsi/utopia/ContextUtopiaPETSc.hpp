// File       : ContextUtopiaPETSc.hpp
// Created    : Fri Oct 18 2024 17:48:35 (+0200)
// Author     : Fabian Wermelinger
// Description: PETSc related function implementation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// no include guard: must only be included in ContextUtopia.h

#include "residual.h"

#include <cassert>
#include <functional>
#ifdef HAS_PETSC
#include <petscmat.h>
#include <petscvec.h>
#endif /* HAS_PETSC */

namespace accel
{
namespace fsi
{

template <size_t N>
void ContextUtopia<N>::convertPETScToNative_(const Mat& mat_seq,
                                             const bool diagonal_part)
{
    using Matrix = typename linearSolver::coefficients<BLOCKSIZE>::Matrix;
    using MatrixElement = typename Matrix::DataType;

    auto& coeff = this->getCoefficients();
    const Index n_owned = coeff.nCoefficients();
    const UtopiaFSIGraph<BLOCKSIZE>* fsi_graph = this->fsi_graph_.get();
    assert(fsi_graph->isBuilt());
    assert(fsi_graph == coeff.getGraph());
    assert(fsi_graph != this->fluid_graph_);

#ifndef DEBUG_GRAPH
    assert(n_owned == this->fluid_graph_->nOwnedNodes());
#endif // DEBUG_GRAPH
    PetscBool petsc_ok = PETSC_FALSE;
    PetscInt petsc_rows = 0;

    PetscScalar* petsc_src = nullptr;
    const PetscInt* petsc_row_ptr = nullptr;
    const PetscInt* petsc_indices = nullptr;

    MatGetRowIJ(mat_seq,
                0,
                PETSC_FALSE,
                PETSC_FALSE,
                &petsc_rows,
                &petsc_row_ptr,
                &petsc_indices,
                &petsc_ok);
    assert(petsc_ok);
    assert(petsc_rows / BLOCKSIZE == n_owned);
    assert(petsc_row_ptr);
    assert(petsc_indices);

    MatSeqAIJGetArray(mat_seq, &petsc_src);

    auto dst = coeff.getAMatrix().valuesRef();
    const auto& row_ptr = fsi_graph->offsets();
    const auto& indices = fsi_graph->indices();
    assert(coeff.getAMatrix().getMemoryLayout() &
           linearSolver::GraphLayout::MemoryLayout__BlockRowMajor);

    std::function<Index(const Index)> block_map;
    std::function<Index(const Index)> block_idx;
    const auto& local_block_map = fsi_graph->getLocalBlockMap();
    if (diagonal_part)
    {
        // owned part
        block_map = [](const Index petsc_idx) { return petsc_idx / BLOCKSIZE; };
        block_idx = [](const Index petsc_idx)
        { return petsc_idx - (petsc_idx / BLOCKSIZE) * BLOCKSIZE; };
    }
    else
    {
        // ghost part
        const auto local_to_global = fsi_graph->getPETScGhostLocalToGlobalMap();
        block_map = [n_owned, &local_block_map](const Index petsc_idx)
        {
            assert(petsc_idx < static_cast<Index>(local_block_map.size()));
            return local_block_map[petsc_idx] + n_owned;
        };
        block_idx = [local_to_global](const Index petsc_idx)
        {
            assert(petsc_idx < static_cast<Index>(local_to_global.size()));
            const Index global_idx = local_to_global[petsc_idx];
            return global_idx - (global_idx / BLOCKSIZE) * BLOCKSIZE;
        };
    }

    for (Index i = 0; i < n_owned; i++)
    {
        PetscInt petsc_j[BLOCKSIZE];
        const PetscInt* petsc_block_row_ptr = petsc_row_ptr + BLOCKSIZE * i;
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            petsc_j[k] = petsc_block_row_ptr[k];
        }

        for (Index j_idx = row_ptr[i]; j_idx < row_ptr[i + 1]; j_idx++)
        {
            const Index j_block = indices[j_idx];
            MatrixElement* block_base = &dst[BLOCKSIZE * BLOCKSIZE * j_idx];
            for (Index k = 0; k < BLOCKSIZE; k++)
            {
                MatrixElement* block_row = block_base + BLOCKSIZE * k;
                while (petsc_j[k] < petsc_block_row_ptr[k + 1] &&
                       block_map(petsc_indices[petsc_j[k]]) == j_block)
                {
                    const Index j = block_idx(petsc_indices[petsc_j[k]]);
                    assert(j >= 0);
                    assert(j < BLOCKSIZE);
                    block_row[j] = petsc_src[petsc_j[k]++];
                }
            }
        }
#ifndef NDEBUG
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            assert(petsc_j[k] == petsc_block_row_ptr[k + 1]);
        }
#endif /* NDEBUG */
    }

    MatSeqAIJRestoreArray(mat_seq, &petsc_src);
    petsc_ok = PETSC_FALSE;
    MatRestoreRowIJ(mat_seq,
                    0,
                    PETSC_FALSE,
                    PETSC_FALSE,
                    &petsc_rows,
                    &petsc_row_ptr,
                    &petsc_indices,
                    &petsc_ok);
    assert(petsc_ok);
}

template <size_t N>
typename ContextUtopia<N>::Array
ContextUtopia<N>::PETScRMSError_(const PetscScalar* R)
{
    const auto& coeff = ctx_wrap_->getCoefficients();
    return linearSolver::residual::RMSError<BLOCKSIZE>(
        coeff.getCommunicator(),
        coeff.nCoefficients(),
        coeff.nGlobalCoefficients(),
        R);
}

template <size_t N>
typename ContextUtopia<N>::Array
ContextUtopia<N>::PETScRMSErrorDiagNormalized_(const PetscScalar* R)
{
    const auto& coeff = ctx_wrap_->getCoefficients();

    const Index n_local = coeff.nCoefficients();
    const Index n_global = coeff.nGlobalCoefficients();
    Vec diag;
    MatCreateVecs(hessian_fluid_->raw_type(), nullptr, &diag);
    MatGetDiagonal(hessian_fluid_->raw_type(), diag);

    const PetscScalar* D;
    VecGetArrayRead(diag, &D);
    Array rms = {0}; // root-mean-square of residual vector
    for (Index i = 0; i < n_local; i++)
    {
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            assert(D[BLOCKSIZE * i + k] > static_cast<PetscScalar>(0.0));
            const double r = R[BLOCKSIZE * i + k] / D[BLOCKSIZE * i + k];
            rms[k] += r * r;
        }
    }
    VecRestoreArrayRead(diag, &D);
    VecDestroy(&diag);

    return linearSolver::residual::RMSReduceAndScale<BLOCKSIZE>(
        coeff.getCommunicator(), rms, 1.0 / n_global);
}

template <size_t N>
void ContextUtopia<N>::PETScApplyStrongBC_()
{
    if (!use_strong_velocity_bc_)
    {
        return;
    }

    assert(fluid_graph_);

    std::vector<PetscInt> global_idx;
    const size_t n_nodes = this->zero_dirichlet_vel_idx_.size();
    global_idx.reserve(n_nodes * N); // N == BLOCKSIZE

    auto& rhs = this->getCoefficients().getBVector();
    for (const auto l_idx : zero_dirichlet_vel_idx_)
    {
        const Index g_idx = fluid_graph_->localToGlobalIndex(l_idx);
        for (int i = 0; i < SPATIAL_DIM; i++) // momentum components only
        {
            assert(static_cast<size_t>(l_idx * N + i) < rhs.size());
            rhs[l_idx * N + i] = static_cast<scalar>(0.0); // zero RHS
            global_idx.push_back(g_idx * N + i);
        }
    }

    // zero matrix rows in PETSc
    PetscBool val = PETSC_TRUE;
    MatGetOption( // get current value
        hessian_fluid_->raw_type(),
        MAT_KEEP_NONZERO_PATTERN,
        &val);
    MatSetOption( // freeze non-zero pattern
        hessian_fluid_->raw_type(),
        MAT_KEEP_NONZERO_PATTERN,
        PETSC_TRUE);
    MatZeroRows(hessian_fluid_->raw_type(),
                global_idx.size(),
                global_idx.data(),
                1.0, // set diagonal to 1.0
                nullptr,
                nullptr);
    MatSetOption( // restore value
        hessian_fluid_->raw_type(),
        MAT_KEEP_NONZERO_PATTERN,
        val);
}

} /* namespace fsi */
} /* namespace accel */
