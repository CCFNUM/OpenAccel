// File       : CRSConnectivity.h
// Created    : Wed Oct 02 2024 11:37:50 (+0200)
// Author     : Fabian Wermelinger
// Description: Graph connectivity interface for matrices
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CRSCONNECTIVITY_H
#define CRSCONNECTIVITY_H

#include "CRSNodeGraph.h"

#include <cassert>
#include <mpi.h>
#include <span>
#include <vector>

namespace linearSolver
{

class CRSConnectivity
{
public:
    using Index = CRSNodeGraph::Index;

    CRSConnectivity() = delete;

    CRSConnectivity(const CRSNodeGraph* graph) : graph_(graph)
    {
    }

    virtual ~CRSConnectivity()
    {
        graph_ = nullptr;
    }

    const CRSNodeGraph* getGraph() const
    {
        return graph_;
    }

    MPI_Comm getCommunicator() const
    {
        return graph_->getCommunicator();
    }

    int commRank() const
    {
        return graph_->commRank();
    }

    int commSize() const
    {
        return graph_->commSize();
    }

    inline Index globalRowOffset() const
    {
        return graph_->globalRowOffset();
    };

    inline Index nRows() const
    {
        return graph_->nRows();
    };

    inline Index nGlobalRows() const
    {
        return graph_->nGlobalRows();
    };

    inline Index nnzBlocks() const
    {
        return graph_->nIndices();
    }

    inline unsigned long long nnzGlobalBlocks() const
    {
        return graph_->nGlobalIndices();
    }

    inline const Index* offsetsPtr() const
    {
        return graph_->offsets().data();
    }

    inline const std::vector<Index>& offsetsRef() const
    {
        return graph_->offsets();
    }

    inline const Index* indicesPtr() const
    {
        return graph_->indices().data();
    }

    inline const std::vector<Index>& indicesRef() const
    {
        return graph_->indices();
    }

    inline const Index* diagOffsetPtr() const
    {
        return graph_->diagonalIndicesOffset().data();
    }

    inline const std::vector<Index>& diagOffsetRef() const
    {
        return graph_->diagonalIndicesOffset();
    }

    inline std::span<const Index> rowCols(Index iRow) const
    {
        assert(0 <= iRow);
        assert(iRow < this->nRows());
        const auto& row_ptr = this->offsetsRef();
        return std::span<const Index>(this->indicesRef())
            .subspan(row_ptr[iRow], row_ptr[iRow + 1] - row_ptr[iRow]);
    }

    inline Index localToGlobal(Index localID) const
    {
        return graph_->localToGlobalIndex(localID);
    };

    inline Index globalToLocal(Index globalID) const
    {
        return graph_->globalToLocalIndex(globalID);
    };

protected:
    const CRSNodeGraph* graph_; // never owned by this class
};

} // namespace linearSolver

#endif /* CRSCONNECTIVITY_H */
