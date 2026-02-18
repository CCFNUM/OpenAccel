// File       : CRSNodeGraphSynchronizer.h
// Created    : Mon Apr 15 2024 16:41:14 (+0200)
// Author     : Fabian Wermelinger
// Description: MPI synchronizer for CRSNodeGraph data structure
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CRSNODEGRAPHSYNCHRONIZER_H
#define CRSNODEGRAPHSYNCHRONIZER_H

#include "CRSNodeGraph.h"
#include "coefficients.h"
#include <cassert>
#include <mpi.h>
#include <stdexcept>
#include <type_traits>
#include <vector>

#define MPISYNCH_DEBUG 0

#if MPISYNCH_DEBUG
#include <iostream>
#endif /* MPISYNCH_DEBUG */

namespace linearSolver
{

template <size_t N, typename T>
class CRSNodeGraphSynchronizer
{
public:
    using Index = typename CRSNodeGraph::Index;
    using PackInfo = typename CRSNodeGraph::PackInfo;

    static constexpr Index BLOCKSIZE = N;

    struct PackData
    {
        static constexpr Index BLOCKSIZE = N;
        const PackInfo* info;
        std::vector<T> send_buf;
        std::vector<T> recv_buf;
    };

    CRSNodeGraphSynchronizer() = delete;

    CRSNodeGraphSynchronizer(const CRSNodeGraph* graph)
        : graph_(graph), recv_complete_(0)
    {
        this->resizePackData();
    }

    CRSNodeGraphSynchronizer(const linearSolver::coefficients<N>& coeffs)
        : graph_(coeffs.getGraph()), recv_complete_(0)
    {
        this->resizePackData();
    }

    virtual ~CRSNodeGraphSynchronizer()
    {
        this->waitSend();
    }

    void resizePackData()
    {
        pack_data_.clear();
        const std::vector<PackInfo>& infos = graph_->getPackInfos();
        for (const PackInfo& info : infos)
        {
            assert(info.send_idx.size() > 0 || info.recv_idx.size() > 0);
            pack_data_.push_back(PackData());
            pack_data_.back().info = &info;
        }
        pending_recv_.reserve(pack_data_.size());
        pending_send_.reserve(pack_data_.size());
        recv_packs_.reserve(pack_data_.size());
    }

    const PackInfo* unpack(const PackData& data, std::vector<T>& dst) const
    {
        const PackInfo* info = data.info;

        assert(data.recv_buf.size() > 0);
        assert(BLOCKSIZE * info->recv_idx.size() <= data.recv_buf.size());

        const T* recv_src = data.recv_buf.data();
        for (const Index i : info->recv_idx)
        {
            assert(BLOCKSIZE * (i + 1) - 1 < static_cast<Index>(dst.size()));
            T* recv_dst = &dst[BLOCKSIZE * i];
            for (int k = 0; k < BLOCKSIZE; k++)
            {
                recv_dst[k] = *recv_src++;
            }
        }
        return info;
    }

    void async(const std::vector<T>& src)
    {
        if (pending_recv_.size() > 0)
        {
            throw std::runtime_error(
                "CRSNodeGraphSynchronizer: `async` method called while "
                "previous receives are pending.");
        }

        if (pending_send_.size() > 0)
        {
            MPI_Waitall(pending_send_.size(),
                        pending_send_.data(),
                        MPI_STATUSES_IGNORE);
        }
        send_clear_();
        recv_clear_();

        if (0 == pack_data_.size())
        {
            return;
        }

        // pack data
        for (PackData& pd : pack_data_)
        {
            pack_(pd, src);
        }

        // post receives
        const int rank = graph_->commRank();
        const MPI_Comm comm = graph_->getCommunicator();
        for (PackData& pd : pack_data_)
        {
            const PackInfo* info = pd.info;
            if (pd.recv_buf.size() > 0)
            {
                MPI_Request req;
                MPI_Irecv(pd.recv_buf.data(),
                          pd.recv_buf.size(),
                          MPIDataType<T>::type(),
                          info->remote_rank,
                          info->remote_rank,
                          comm,
                          &req);
                pending_recv_.push_back(req);
                recv_packs_.push_back(&pd);
            }
        }

        // post sends
        for (const PackData& pd : pack_data_)
        {
            const PackInfo* info = pd.info;
            if (pd.send_buf.size() > 0)
            {
                MPI_Request req;
                MPI_Isend(pd.send_buf.data(),
                          pd.send_buf.size(),
                          MPIDataType<T>::type(),
                          info->remote_rank,
                          rank,
                          comm,
                          &req);
                pending_send_.push_back(req);
            }
        }

#if MPISYNCH_DEBUG
        std::vector<const PackData*> avail;
        this->waitAll(avail);

        int size;
        MPI_Comm_size(comm, &size);
        for (int r = 0; r < size; r++)
        {
            for (const PackData* pd : avail)
            {
                if (r == rank)
                {
                    std::cout << "RANK " << rank
                              << ": remote_rank=" << pd->info->remote_rank
                              << '\n';
                    std::cout << "\tn nodes:\t" << graph_->nOwnedNodes()
                              << '\n';
                    std::cout << "\tn ghosts:\t" << graph_->nGhostNodes()
                              << '\n';
                    std::cout << "\tn active:\t" << graph_->nAllNodes() << '\n';
                    std::cout << "\tsend_idx (" << pd->info->send_idx.size()
                              << "):";
                    for (const Index i : pd->info->send_idx)
                    {
                        std::cout << '\t' << i;
                    }
                    std::cout << "\n\tsend_buf (" << pd->send_buf.size()
                              << "):";
                    for (const T v : pd->send_buf)
                    {
                        std::cout << '\t' << v;
                    }
                    std::cout << "\n\trecv_idx (" << pd->info->recv_idx.size()
                              << "):";
                    for (const Index i : pd->info->recv_idx)
                    {
                        std::cout << '\t' << i;
                    }
                    std::cout << "\n\trecv_buf (" << pd->recv_buf.size()
                              << "):";
                    for (const T v : pd->recv_buf)
                    {
                        std::cout << '\t' << v;
                    }
                    std::cout << '\n';
                }
            }
            MPI_Barrier(comm);
        }
#endif /* MPISYNCH_DEBUG */
    }

    void getSendOnlyPacks(std::vector<const PackData*>& avail)
    {
        avail.clear();
        if (pending_send_.size() > 0)
        {
            for (const PackData& pd : pack_data_)
            {
                if (pd.send_buf.size() > 0 && pd.recv_buf.size() == 0)
                {
                    avail.push_back(&pd);
                }
            }
        }
    }

    void waitSome(std::vector<const PackData*>& avail)
    {
        avail.clear();
        if (recv_complete_ == pending_recv_.size())
        {
            recv_clear_();
            return;
        }

        available_.resize(pending_recv_.size());
        int n_avail = 0;
        MPI_Waitsome(pending_recv_.size(),
                     pending_recv_.data(),
                     &n_avail,
                     available_.data(),
                     MPI_STATUSES_IGNORE);
        recv_complete_ += n_avail;
        avail.resize(n_avail);
        for (int i = 0; i < n_avail; i++)
        {
            avail[i] = recv_packs_[available_[i]];
        }
    }

    void waitAll(std::vector<const PackData*>& avail)
    {
        MPI_Waitall(
            pending_recv_.size(), pending_recv_.data(), MPI_STATUSES_IGNORE);
        avail.resize(pending_recv_.size());
        assert(pending_recv_.size() == recv_packs_.size());
        for (int i = 0; i < static_cast<int>(recv_packs_.size()); i++)
        {
            avail[i] = recv_packs_[i];
        }
        recv_clear_();
    }

    void waitAll(std::vector<T>& dst)
    {
        std::vector<const PackData*> avail;
        this->waitAll(avail);
        for (const PackData* pd : avail)
        {
            this->unpack(*pd, dst);
        }
    }

    void waitSend()
    {
        if (pending_send_.size() > 0)
        {
            MPI_Waitall(pending_send_.size(),
                        pending_send_.data(),
                        MPI_STATUSES_IGNORE);
        }
        send_clear_();
    }

private:
    const CRSNodeGraph* graph_;
    std::vector<PackData> pack_data_;
    std::vector<MPI_Request> pending_recv_;
    std::vector<MPI_Request> pending_send_;

    size_t recv_complete_;
    std::vector<int> available_;
    std::vector<const PackData*> recv_packs_;

    void recv_clear_()
    {
        recv_complete_ = 0;
        pending_recv_.clear();
        recv_packs_.clear();
    }

    void send_clear_()
    {
        pending_send_.clear();
    }

    void pack_(PackData& data, const std::vector<T>& src)
    {
        const PackInfo* info = data.info;

        data.recv_buf.resize(BLOCKSIZE * info->recv_idx.size());

        data.send_buf.resize(BLOCKSIZE * info->send_idx.size());
        T* send_dst = data.send_buf.data();
        for (const Index i : info->send_idx)
        {
            assert(BLOCKSIZE * (i + 1) - 1 < BLOCKSIZE * graph_->nOwnedNodes());
            const T* send_src = &src[BLOCKSIZE * i];
            for (int k = 0; k < BLOCKSIZE; k++)
            {
                *send_dst++ = send_src[k];
            }
        }
    }
};

} /* namespace linearSolver */

#undef MPISYNCH_DEBUG

#endif /* CRSNODEGRAPHSYNCHRONIZER_H */
