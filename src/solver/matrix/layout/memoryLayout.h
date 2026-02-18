// File       : memoryLayout.h
// Created    : Tue Oct 01 2024 10:27:55 (+0200)
// Author     : Fabian Wermelinger
// Description: CRS matrix memory layout transformations
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MEMORYLAYOUT_H
#define MEMORYLAYOUT_H

#include "CRSMatrix.h"
#include "CRSNodeGraph.h"
#include <cassert>
#include <span>

namespace linearSolver
{
namespace matrixLayout
{

template <typename Matrix>
void blockRowToRowMajor(const typename Matrix::Index block_row_idx,
                        const Matrix& A,
                        std::span<const typename Matrix::Index> block_col_idx,
                        std::vector<typename Matrix::Index>& row_nnz,
                        std::vector<typename Matrix::Index>& row_idx,
                        std::vector<typename Matrix::Index>& col_idx,
                        std::vector<typename Matrix::DataType>& values)
{
    using Index = typename Matrix::Index;
    using DataType = typename Matrix::DataType;
    constexpr Index BLOCKSIZE = Matrix::BLOCKSIZE;

    // TODO: [2024-10-01] if BLOCKSIZE == 1 can do direct
    // copy!

    assert(A.getMemoryLayout() & GraphLayout::MemoryLayout__BlockRowMajor);

    const auto flat_block_row_values = A.rowVals(block_row_idx);
    assert(block_col_idx.size() ==
           flat_block_row_values.size() / (BLOCKSIZE * BLOCKSIZE));

    row_nnz.resize(BLOCKSIZE);
    row_idx.resize(BLOCKSIZE);
    col_idx.resize(flat_block_row_values.size());
    values.resize(flat_block_row_values.size());

    const Index idx_offset = A.localToGlobal(0) * BLOCKSIZE;
    const Index n_col_blocks = block_col_idx.size();
    const Index block_nnz = n_col_blocks * BLOCKSIZE;
    for (Index k = 0; k < BLOCKSIZE; k++)
    {
        row_nnz[k] = block_nnz;
        row_idx[k] = idx_offset + block_row_idx * BLOCKSIZE + k;
        Index* idx = &col_idx[k * block_nnz];
        DataType* dst = &values[k * block_nnz];
        for (Index j = 0; j < n_col_blocks; j++)
        {
            const DataType* src =
                &flat_block_row_values[(j * BLOCKSIZE + k) * BLOCKSIZE];
            for (Index l = 0; l < BLOCKSIZE; l++)
            {
                *idx++ = block_col_idx[j] * BLOCKSIZE + l;
                *dst++ = src[l];
            }
        }
    }
}

} /* namespace matrixLayout */
} /* namespace linearSolver */

#endif /* MEMORYLAYOUT_H */
