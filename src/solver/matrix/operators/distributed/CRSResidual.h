// File       : CRSResidual.h
// Created    : Thu Apr 18 2024 13:09:04 (+0200)
// Author     : Fabian Wermelinger
// Description: Distributed residual computation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CRSRESIDUAL_H
#define CRSRESIDUAL_H

#include "CRSMatrix.h"
#include "CRSNodeGraph.h"
#include "CRSNodeGraphSynchronizer.h"
#include "blockMatrixOperators.h"
#include <algorithm>
#include <cassert>

#define BLOCKING 1

namespace linearSolver
{

template <size_t N>
void CRSResidual(const CRSMatrix<N>& A,
                 const typename CRSMatrix<N>::Vector& b,
                 typename CRSMatrix<N>::Vector& x,
                 typename CRSMatrix<N>::Vector& r)
{
    using Index = typename CRSMatrix<N>::Index;
    using Real = typename CRSMatrix<N>::DataType;
    static constexpr Index BLOCKSIZE = N;

    assert(static_cast<const void*>(x.data()) !=
           static_cast<const void*>(r.data()));
    assert(BLOCKSIZE * A.nRows() <= static_cast<Index>(b.size()));
    assert(BLOCKSIZE * A.nRows() <= static_cast<Index>(r.size()));
    // TODO: [2024-04-18] this could be relaxed if recv_buf
    // is used directly in algorithm
    assert(static_cast<Index>(x.size()) ==
           BLOCKSIZE * A.getGraph()->nAllNodes());

    const CRSNodeGraph* graph = A.getGraph();
    CRSNodeGraphSynchronizer<N, Real> sync(graph);
    sync.async(x);

    // r = b - Ax
#if BLOCKING
    sync.waitAll(x);
    for (Index i = 0; i < A.nRows(); i++)
    {
        const auto local_idx = graph->rowLocalIndices(i);
        const auto coeffs = A.rowVals(i);

        Real* r_blk = &r[BLOCKSIZE * i];
        const Real* b_blk = &b[BLOCKSIZE * i];
        BlockMatrix::vectorAssign<BLOCKSIZE>(b_blk, r_blk);
        for (Index j = 0; j < static_cast<Index>(local_idx.size()); j++)
        {
            const Real* a_blk = &coeffs[BLOCKSIZE * BLOCKSIZE * j];
            const Real* x_blk = &x[BLOCKSIZE * local_idx[j]];
            BlockMatrix::matrixVectorSub<BLOCKSIZE>(a_blk, x_blk, r_blk);
        }
    }
#else
    // WIP: [2024-04-18] asynchronous
#endif /* BLOCKING */
}

} /* namespace linearSolver */

#undef BLOCKING

#endif /* CRSRESIDUAL_H */
