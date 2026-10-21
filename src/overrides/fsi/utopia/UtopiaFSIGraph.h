// File       : UtopiaFSIGraph.h
// Created    : Wed Oct 16 2024 10:21:42 (+0200)
// Author     : Fabian Wermelinger
// Description: Utopia FSI node graph
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef UTOPIAFSIGRAPH_H
#define UTOPIAFSIGRAPH_H

#include "CRSNodeGraph.h"
#include "utopia_fe_base.hpp"

#include <span>
#include <vector>
#ifdef HAS_PETSC
#include <petscmat.h>
#endif /* HAS_PETSC */

namespace accel
{

namespace fsi
{

template <size_t N>
class UtopiaFSIGraph : public linearSolver::CRSNodeGraph
{
    using Index = typename linearSolver::CRSNodeGraph::Index;
    using UMatrix = typename utopia::USparseMatrix;
    static constexpr Index BLOCKSIZE = N;

public:
    UtopiaFSIGraph() = delete;

    UtopiaFSIGraph(const MPI_Comm comm, const linearSolver::GraphLayout layout)
        : CRSNodeGraph(comm, layout), n_ghosts_(0), o_local_to_global_(nullptr)
    {
    }

    void buildGraph(const UMatrix& m);

    std::span<const PetscInt> getPETScGhostLocalToGlobalMap() const
    {
        assert(o_local_to_global_);
        assert(n_ghosts_ > 0);
        return std::span<const PetscInt>(o_local_to_global_, n_ghosts_);
    }

    Index getPETScNumberOfGhosts() const
    {
        return n_ghosts_;
    }

    const std::vector<Index>& getLocalBlockMap() const
    {
        return o_local_block_map_;
    }

protected:
    using linearSolver::CRSNodeGraph::primary_indices_;
    using linearSolver::CRSNodeGraph::row_ptr_;
    using linearSolver::CRSNodeGraph::secondary_indices_;

    void buildGraph_() override;

private:
    Mat d_petsc_; // diagonal matrix part
    Mat o_petsc_; // off-diagonal matrix part

    // PETSc ghost management
    PetscInt n_ghosts_;
    const PetscInt* o_local_to_global_;
    std::vector<Index> o_local_block_map_;
};

} /* namespace fsi */
} /* namespace accel */

#include "UtopiaFSIGraph.hpp"

#endif // UTOPIAFSIGRAPH_H
