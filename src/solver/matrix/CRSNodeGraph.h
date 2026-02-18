// File       : CRSNodeGraph.h
// Created    : Sat Apr 13 2024 23:13:14 (+0200)
// Author     : Mhamad Mahdi Alloush, Fabian Wermelinger
// Description: Node graph generator base class
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CRSNODEGRAPH_H
#define CRSNODEGRAPH_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <span>
#include <vector>

#ifdef GRAPH_INDEX_64BIT
typedef int64_t TGraphIndex;
#else
typedef int32_t TGraphIndex;
#endif /* GRAPH_INDEX_64BIT */

namespace linearSolver
{

// MPI type handles
template <typename T>
struct MPIDataType
{
    static inline MPI_Datatype type()
    {
        return MPI_DATATYPE_NULL;
    }
};

template <>
struct MPIDataType<double>
{
    static inline MPI_Datatype type()
    {
        return MPI_DOUBLE;
    }
};

template <>
struct MPIDataType<float>
{
    static inline MPI_Datatype type()
    {
        return MPI_FLOAT;
    }
};

template <>
struct MPIDataType<unsigned int>
{
    static inline MPI_Datatype type()
    {
        return MPI_UNSIGNED;
    }
};

template <>
struct MPIDataType<int32_t>
{
    static inline MPI_Datatype type()
    {
        return MPI_INT32_T;
    }
};

template <>
struct MPIDataType<int64_t>
{
    static inline MPI_Datatype type()
    {
        return MPI_INT64_T;
    }
};

struct MemoryFootprint
{
    static constexpr char UNIT[] = "MB";
    static constexpr unsigned int BYTE_DIVIDE = (1 << 20); // MB
    unsigned long long int sum_byte;
};

class GraphLayout
{
public:
    enum : int
    {
        // identifiers must be powers of two
        //
        // graph layout
        ColumnIndexOrder__Local = 0x001,
        ColumnIndexOrder__Global = 0x002,
        Stencil__Reduced = 0x004,
        // memory layout
        MemoryLayout__ConvertInplace = 0x008,
        MemoryLayout__BlockRowMajor = 0x010,
        MemoryLayout__RowMajor = 0x020,
        // convenience
        MemoryLayout__PETScAIJ = MemoryLayout__RowMajor,
        MemoryLayout__HYPRE_IJ = MemoryLayout__RowMajor,
        MemoryLayout__AllBits =
            MemoryLayout__BlockRowMajor | MemoryLayout__RowMajor,
        ColumnIndexOrder__AllBits =
            ColumnIndexOrder__Local | ColumnIndexOrder__Global,
        UNDEFINED = 0x000,
        Stencil__Full = UNDEFINED
    };

    GraphLayout() : value_(UNDEFINED)
    {
    }

    GraphLayout(const int c) : value_(c)
    {
    }

    operator int()
    {
        return value_;
    }

    operator int() const
    {
        return value_;
    }

private:
    int value_;
};

class CRSNodeGraph
{
public:
    // Main graph index type.  This type is used to compute differences and
    // should be a signed integral.
    using Index = TGraphIndex;

    struct PackInfo
    {
        int remote_rank;
        std::vector<Index> send_idx;
        std::vector<Index> recv_idx;
    };

    CRSNodeGraph() = delete;

    explicit CRSNodeGraph(const MPI_Comm comm, const GraphLayout layout);

    virtual ~CRSNodeGraph() = default;

    void buildGraph();

    inline MPI_Comm getCommunicator() const
    {
        return comm_;
    }

    inline int commRank() const
    {
        return rank_;
    }

    inline int commSize() const
    {
        return size_;
    }

    inline Index nOwnedNodes() const
    {
        assert(is_built_);
        assert(n_owned_nodes_ > 0);
        return n_owned_nodes_;
    }

    inline Index nGhostNodes() const
    {
        assert(is_built_);
        return n_ghost_nodes_;
    }

    inline Index nAllNodes() const
    {
        assert(is_built_);
        return nOwnedNodes() + nGhostNodes();
    }

    inline Index nGlobalNodes() const
    {
        assert(is_built_);
        return global_number_nodes_;
    }

    inline Index globalRowOffset() const
    {
        assert(is_built_);
        return global_row_offset_;
    }

    inline Index nRows() const
    {
        return nOwnedNodes();
    }

    inline Index nGlobalRows() const
    {
        return nGlobalNodes();
    }

    inline Index nIndices() const
    {
        return row_ptr_.back();
    }

    inline unsigned long long nGlobalIndices() const
    {
        assert(is_built_);
        return global_number_indices_;
    }

    inline const std::vector<Index>& offsets() const
    {
        assert(is_built_);
        assert(row_ptr_.size() > 1);
        return row_ptr_;
    }

    inline const std::vector<Index>& indices() const
    {
        assert(is_built_);
        assert(static_cast<Index>(primary_indices_.size()) == this->nIndices());
        return primary_indices_;
    }

    inline const std::vector<Index>& diagonalIndicesOffset() const
    {
        assert(is_built_);
        assert(static_cast<Index>(diagonal_row_offset_.size()) ==
               this->nRows());
        return diagonal_row_offset_;
    }

    inline const std::vector<PackInfo>& getPackInfos() const
    {
        assert(is_built_);
        return pack_infos_;
    }

    inline GraphLayout getLayout() const
    {
        return layout_;
    }

    inline bool isLocalColumnOrder() const
    {
        return (layout_ & GraphLayout::ColumnIndexOrder__Local);
    }

    inline bool isGlobalColumnOrder() const
    {
        return (layout_ & GraphLayout::ColumnIndexOrder__Global);
    }

    inline bool isBuilt() const
    {
        return is_built_;
    }

    inline Index nnzOwned(const Index i_row) const
    {
        assert(is_built_);
        assert(0 <= i_row);
        assert(i_row < this->nRows());
        return row_nnz_owned_[i_row];
    }

    inline std::span<const Index> nnzOwned() const
    {
        assert(is_built_);
        return std::span<const Index>(row_nnz_owned_);
    }

    inline Index nnzGhost(const Index i_row) const
    {
        assert(is_built_);
        assert(0 <= i_row);
        assert(i_row < this->nRows());
        return row_nnz_ghost_[i_row];
    }

    inline std::span<const Index> nnzGhost() const
    {
        assert(is_built_);
        return std::span<const Index>(row_nnz_ghost_);
    }

    inline Index localToGlobalIndex(const Index j_local) const
    {
        // only defined if `j_local` is owned
        // clang-format off
        assert(0 <= j_local &&
               "CRSNodeGraph::localToGlobalIndex: `j_local` must be owned");
        assert(j_local <= this->nOwnedNodes() &&
               "CRSNodeGraph::localToGlobalIndex: `j_local` must be owned (inclusive)"); // inclusive
        // clang-format on
        return j_local + this->global_row_offset_;
    }

    inline Index globalToLocalIndex(const Index j_global) const
    {
        // only defined if `j_global` is owned
        // clang-format off
        assert(this->global_row_offset_ <= j_global &&
               "CRSNodeGraph::localToGlobalIndex: `j_global` must be owned");
        assert(j_global <= this->global_row_offset_ + this->nOwnedNodes() &&
               "CRSNodeGraph::localToGlobalIndex: `j_global` must be owned (inclusive)"); // inclusive
        // clang-format on
        return j_global - this->global_row_offset_;
    }

    std::span<const Index> localIndices() const;
    std::span<const Index> globalIndices() const;

    std::span<const Index> rowLocalIndices(const Index i_row) const;
    std::span<const Index> rowGlobalIndices(const Index i_row) const;

    void getMemoryFootprint(MemoryFootprint& data) const;
    void serialize(std::ofstream& out) const;
    void deserialize(std::ifstream& in);

protected:
    MPI_Comm comm_;
    int rank_;
    int size_;

    Index n_owned_nodes_; // locally owned nodes
    Index n_ghost_nodes_; // nodes in graph not locally owned by process
    bool determine_n_ghosts_;
    bool is_built_;

    // CRS data structures
    std::vector<Index> row_ptr_;           // row offsets
    std::vector<Index> primary_indices_;   // main column index order (sorted)
    std::vector<Index> secondary_indices_; // complement column index order

    // per row nnz
    std::vector<Index> row_nnz_owned_; // nnz per row (owned)
    std::vector<Index> row_nnz_ghost_; // nnz per row (total ghosts)

    // diagonal
    std::vector<Index> diagonal_row_offset_; // index into *_indices_

    // MPI
    Index global_row_offset_;
    Index global_number_nodes_;
    unsigned long long global_number_indices_;
    std::vector<PackInfo> pack_infos_;

    virtual void buildGraph_() = 0;
    virtual void computePackInfos_();
    virtual void computeDiagonalIndices_();

    void computeAuxiliaryData_();

    void computeGlobalRowOffset_();
    void computeGlobalGraphSize_();
    void computeRowNonZeros_();

    void sortPrimaryIndices_();

    static inline void permuteInPlace_(Index* indices,
                                       const std::vector<Index>& permute)
    {
        std::vector<bool> done(permute.size());
        for (Index i = 0; i < static_cast<Index>(permute.size()); ++i)
        {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            Index prev_j = i;
            Index j = permute[i];
            while (i != j)
            {
                std::swap(indices[prev_j], indices[j]);
                done[j] = true;
                prev_j = j;
                j = permute[j];
            }
        }
    }

    static inline void permuteCopy_(Index* indices,
                                    std::vector<Index>& buffer,
                                    const std::vector<Index>& permute)
    {
        buffer.resize(permute.size());
        for (Index i = 0; i < static_cast<Index>(permute.size()); i++)
        {
            buffer[i] = indices[permute[i]];
        }
        std::copy(buffer.begin(), buffer.end(), indices);
    }

    void setup_();
    void resetGraph_();

private:
    const GraphLayout layout_;

    Index filterGhostsForOwnerRank_(
        int& owner_rank,
        std::vector<Index>::const_iterator ghost_start, // must be sorted
        const std::vector<Index>::const_iterator ghost_end,
        const int myrank, // of caller
        const int size,
        const std::vector<Index>& n_local_nodes) const;
};

} /* namespace linearSolver */

#endif /* CRSNODEGRAPH_H */
