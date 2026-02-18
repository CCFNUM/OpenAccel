// File       : CRSMatrixVector.h
// Created    : Wed Apr 17 2024 16:03:52 (+0200)
// Author     : Fabian Wermelinger
// Description: Distributed Matrix-Vector product
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CRSMATRIXVECTOR_H
#define CRSMATRIXVECTOR_H

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
void CRSMatrixVector(const CRSMatrix<N>& A,
                     typename CRSMatrix<N>::Vector& x,
                     typename CRSMatrix<N>::Vector& y)
{
    using Index = typename CRSMatrix<N>::Index;
    using Real = typename CRSMatrix<N>::DataType;
    static constexpr Index BLOCKSIZE = N;

    assert(static_cast<const void*>(x.data()) !=
           static_cast<const void*>(y.data()));
    assert(BLOCKSIZE * A.nRows() <= static_cast<Index>(y.size()));
    // TODO: [2024-04-18] this could be relaxed if recv_buf
    // is used directly in algorithm
    assert(static_cast<Index>(x.size()) ==
           BLOCKSIZE * A.getGraph()->nAllNodes());

    const CRSNodeGraph* graph = A.getGraph();
    CRSNodeGraphSynchronizer<N, Real> sync(graph);
    sync.async(x);
    std::fill(y.begin(), y.end(), static_cast<Real>(0));

#if BLOCKING
    sync.waitAll(x);
    for (Index i = 0; i < A.nRows(); i++)
    {
        const auto local_idx = graph->rowLocalIndices(i);
        const auto coeffs = A.rowVals(i);
        Real* y_blk = &y[BLOCKSIZE * i];
        for (Index j = 0; j < static_cast<Index>(local_idx.size()); j++)
        {
            const Real* a_blk = &coeffs[BLOCKSIZE * BLOCKSIZE * j];
            const Real* x_blk = &x[BLOCKSIZE * local_idx[j]];
            BlockMatrix::matrixVectorAdd<BLOCKSIZE>(a_blk, x_blk, y_blk);
        }
    }
#else
    // WIP: [2024-04-18] asynchronous
    // for (auto i : container) {
    //     // owned coeffs
    // }

    // std::vector<const PackData *> avail;
    // do {
    //     sync.waitSome(avail);
    //     for (const PackData *pd : avail) {
    //         const PackInfo *info = sync.unpack(*pd, x);
    //         // update ghost
    //     }
    // } while (avail.size() > 0);
#endif /* BLOCKING */
}

} /* namespace linearSolver */

#undef BLOCKING

#endif /* CRSMATRIXVECTOR_H */
