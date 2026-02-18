// File       : CRSMatrix.hpp
// Created    : Mon Sep 30 2024 16:29:54 (+0200)
// Author     : Fabian Wermelinger
// Description: CRS matrix implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace linearSolver
{

template <size_t N>
template <typename T>
std::vector<char> CRSMatrix<N>::gatherToGlobal_(const T* values,
                                                const size_t n_values,
                                                const bool exscan) const
{
    std::vector<T> val(n_values);
    std::copy(values, values + n_values, val.begin());

    const int size = this->commSize();
    const int rank = this->commRank();
    const MPI_Comm comm = this->getCommunicator();
    if (exscan)
    {
        const int last = val.back();
        int offset = 0;
        MPI_Exscan(&last, &offset, 1, MPI_INT, MPI_SUM, comm);
        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] += offset;
        }
        if (rank < size - 1)
        {
            val.pop_back();
        }
    }
    const int local_n = val.size() * sizeof(T);
    std::vector<int> bytes(size, 0);
    MPI_Gather(&local_n, 1, MPI_INT, bytes.data(), 1, MPI_INT, 0, comm);

    // if (0 == rank) {
    //     std::cout << "Size = " << sizeof(T) << std::endl;
    // }

    std::vector<int> displ;
    int offset = 0;
    for (const int b : bytes)
    {
        displ.push_back(offset);
        offset += b;
    }

    std::vector<char> ret(offset);
    MPI_Gatherv(val.data(),
                local_n,
                MPI_BYTE,
                ret.data(),
                bytes.data(),
                displ.data(),
                MPI_BYTE,
                0,
                comm);
    return ret;
}

template <size_t N>
void CRSMatrix<N>::dumpBinary_(const std::string fname,
                               const std::vector<char>& v) const
{
    if (0 == this->commRank())
    {
        std::ofstream out(fname, std::ios::binary);
        out.write(v.data(), v.size());
        out.close();
    }
}

template <size_t N>
std::ostream& CRSMatrix<N>::stream_(std::ostream& os,
                                    Index max_rows,
                                    Index max_cols,
                                    const Index width,
                                    const Index precision) const
{
    std::ios og_fmt(nullptr);
    og_fmt.copyfmt(os);
    os << std::scientific << std::setprecision(precision);

    static constexpr Index BLOCKSIZE = N;
    const MPI_Comm comm = this->getCommunicator();
    max_rows = (max_rows == 0) ? this->nRows() : max_rows;
    max_cols = (max_cols == 0) ? this->nRows() : max_cols;
    if (this->commSize() > 1)
    {
        for (Index irank = 0; irank < this->commSize(); irank++)
        {
            if (irank == this->commRank())
            {
                os << "RANK: " << irank << '\n';

                os << "  " << std::setw(width) << " ";
                for (Index j = 0; j < max_cols * BLOCKSIZE; j++)
                {
                    os << std::setw(width) << j << " ";
                }
                os << "\n\n";
                for (Index i = 0; i < max_rows * BLOCKSIZE; i++)
                {
                    os << std::setw(width) << i << "  ";
                    for (Index j = 0; j < max_cols * BLOCKSIZE; j++)
                    {
                        os << std::setw(width) << (*this)(i, j) << " ";
                    }
                    os << '\n';
                }
            }
            MPI_Barrier(comm);
        }
    }
    else
    {
        os << "  " << std::setw(width) << " ";
        for (Index j = 0; j < max_cols * BLOCKSIZE; j++)
        {
            os << std::setw(width) << j << " ";
        }
        os << "\n\n";
        for (Index i = 0; i < max_rows * BLOCKSIZE; i++)
        {
            os << std::setw(width) << i << "  ";
            for (Index j = 0; j < max_cols * BLOCKSIZE; j++)
            {
                os << std::setw(width) << (*this)(i, j) << " ";
            }
            os << '\n';
        }
    }
    os.copyfmt(og_fmt);
    return os;
}

template <size_t N>
void CRSMatrix<N>::dump(Index maxRow,
                        Index maxCol,
                        Index width,
                        Index precision) const
{
    this->stream_(std::cout, maxRow, maxCol, width, precision);
}

template <size_t N>
void CRSMatrix<N>::dumpRow(Index rowID, Index width, Index precision) const
{
    const MPI_Comm comm = this->getCommunicator();
    if (this->commSize() > 1)
    {
        throw std::runtime_error("not implemented for parallel");
    }
    else
    {
        auto rowVals = this->rowVals(rowID);
        auto rowCols = this->rowCols(rowID);

        for (Index j = 0; j < rowCols.size(); j++)
        {
            std::cout << "col id " << rowCols[j] << ":\n\n";

            for (Index ii = 0; ii < N; ii++)
            {
                for (Index jj = 0; jj < N; jj++)
                {
                    std::cout << std::scientific << std::setprecision(precision)
                              << std::setw(width)
                              << rowVals[j * N * N + ii * N + jj] << "\t";
                }

                std::cout << "\n";
            }

            std::cout << "\n\n\n";
        }
    }
}

template <size_t N>
void CRSMatrix<N>::writeMatrix_(const char* name) const
{
    // convert layout from blocked to row based
    std::vector<Index> offsets_row(BLOCKSIZE * nRows() + 1);
    std::vector<Index> indices_row(this->nnz());
    std::vector<DataType> values_row(this->nnz());

    Index* dst_row = &offsets_row[0];
    Index* dst_col = &indices_row[0];
    DataType* dst_val = &values_row[0];
    Index col_count = 0;
    for (Index i = 0; i < nRows(); i++)
    {
        const auto block_col_idx = rowCols(i);
        const auto flat_block_row_values = rowVals(i);
        const Index n_col_blocks = block_col_idx.size();
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            *dst_row++ = col_count;
            for (Index j = 0; j < n_col_blocks; j++)
            {
                const DataType* src =
                    &flat_block_row_values[(j * BLOCKSIZE + k) * BLOCKSIZE];
                for (Index l = 0; l < BLOCKSIZE; l++)
                {
                    *dst_col++ = block_col_idx[j] * BLOCKSIZE + l;
                    *dst_val++ = src[l];
                    ++col_count;
                }
            }
        }
    }
    assert(col_count == static_cast<Index>(this->nnz()));
    *dst_row = col_count;

    const std::string basename(name);

#if 1
    this->dumpBinary_(
        basename + "_rows.bin",
        this->gatherToGlobal_(offsets_row.data(), offsets_row.size(), true));
    this->dumpBinary_(
        basename + "_cols.bin",
        this->gatherToGlobal_(indices_row.data(), indices_row.size()));
    this->dumpBinary_(
        basename + "_vals.bin",
        this->gatherToGlobal_(values_row.data(), values_row.size()));
#else
    // non-flattened block-based indexing (only relevant for testing purposes)
    const auto& offsets = this->offsetsRef();
    const auto& indices = this->indicesRef();
    const auto& values = this->valuesRef();
    this->dumpBinary_(
        basename + "_rows.bin",
        this->gatherToGlobal_(offsets.data(), offsets.size(), true));
    this->dumpBinary_(basename + "_cols.bin",
                      this->gatherToGlobal_(indices.data(), indices.size()));
    this->dumpBinary_(basename + "_vals.bin",
                      this->gatherToGlobal_(values.data(), values.size()));
#endif
}

template <size_t N>
typename CRSMatrix<N>::DataType CRSMatrix<N>::operator()(Index i, Index j) const
{
    Index I = i / N; // block row index
    Index J = j / N; // block col index

    auto rowVals = this->rowVals(I);
    auto rowCols = this->rowCols(I);

    Index length = rowCols.size();
    Index offset = 0;
    while (J != rowCols[offset] && offset < length)
    {
        offset++;
    }

    if (offset >= length)
    {
        return 0;
    }

    Index blockOffset = (i - I * N) * N + (j - J * N);

    return rowVals[N * N * offset + blockOffset];
};

template <size_t N>
typename CRSMatrix<N>::Index CRSMatrix<N>::bandwidth() const
{
    assert(this->commSize() == 1 || this->graph_->isGlobalColumnOrder());

    // NOTE [faw 2025-12-15]: serial implementation

    Index beta_max = 0;
    for (Index i = 0; i < this->nRows(); i++)
    {
        // the following search would not be necessary if rowCols(i) is sorted
        // (we do it anyway to be independent of that property)
        Index j_min = std::numeric_limits<Index>::max();
        for (const Index j : this->rowCols(i))
        {
            j_min = j < j_min ? j : j_min;
        }
        const Index beta = i - j_min;
        beta_max = beta > beta_max ? beta : beta_max;
    }
    return beta_max;
}

template <size_t N>
typename CRSMatrix<N>::Index CRSMatrix<N>::profile() const
{
    assert(this->commSize() == 1 || this->graph_->isGlobalColumnOrder());

    // NOTE [faw 2025-12-15]: serial implementation

    Index envelope = 0;
    for (Index i = 0; i < this->nRows(); i++)
    {
        // the following search would not be necessary if rowCols(i) is sorted
        // (we do it anyway to be independent of that property)
        Index j_min = std::numeric_limits<Index>::max();
        for (const Index j : this->rowCols(i))
        {
            j_min = j < j_min ? j : j_min;
        }
        envelope += i - j_min;
    }
    return envelope;
}

template <size_t N>
typename CRSMatrix<N>::DataType norm__frobenius(const CRSMatrix<N>* A)
{
    using TReal = CRSMatrix<N>::DataType;
    TReal norm = 0;
    for (const TReal a_ij : A->valuesRef())
    {
        norm += a_ij * a_ij;
    }

    TReal gnorm = 0;
    MPI_Allreduce(&norm,
                  &gnorm,
                  1,
                  ::linearSolver::MPIDataType<TReal>::type(),
                  MPI_SUM,
                  A->getCommunicator());

    return std::sqrt(gnorm);
}

template <size_t N>
typename CRSMatrix<N>::DataType CRSMatrix<N>::norm(const MatrixNorm type) const
{
    DataType norm = 0;
    if (MatrixNorm::Frobenius == type)
    {
        norm = norm__frobenius(this);
    }
    else
    {
        throw std::runtime_error("Requested matrix norm is not implemented");
    }
    return norm;
}

} // namespace linearSolver
