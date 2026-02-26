// File       : messager.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "messager.h"
#include <iostream>

namespace accel
{

MPI_Comm messager::comm_ = MPI_COMM_WORLD;

MPI_Comm messager::comm(MPI_Comm comm)
{
    return comm;
}

bool messager::parallel(MPI_Comm comm)
{
    return nProcs(comm) > 1;
}

bool messager::master(MPI_Comm comm)
{
    return myProcNo(comm) == 0;
}

int messager::myProcNo(MPI_Comm comm)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

int messager::nProcs(MPI_Comm comm)
{
    int sz = 1;
    MPI_Comm_size(comm, &sz);
    return sz;
}

void messager::barrier(MPI_Comm comm)
{
    MPI_Barrier(comm);
}

void messager::exitAll(MPI_Comm comm)
{
    MPI_Barrier(comm);
    exit(1);
}

void messager::markProcReached(MPI_Comm comm)
{
    std::cout << myProcNo() << " reached here" << std::endl;
}

void messager::markProcReachedAndExit(MPI_Comm comm)
{
    markProcReached(comm);
    exitAll(comm);
}

void messager::print(std::string val, MPI_Comm comm)
{
    std::cout << "[" << myProcNo(comm) << "]\t" << val << std::endl;
}

void messager::print(int val, MPI_Comm comm)
{
    std::cout << "[" << myProcNo(comm) << "]\t" << val << std::endl;
}

void messager::print(double val, MPI_Comm comm)
{
    std::cout << "[" << myProcNo(comm) << "]\t" << val << std::endl;
}

void messager::procPrint(int val, int proc, MPI_Comm comm)
{
    if (myProcNo(comm) == proc)
    {
        std::cout << "[" << proc << "]\t" << val << std::endl;
    }
}

void messager::procPrint(double val, int proc, MPI_Comm comm)
{
    if (myProcNo(comm) == proc)
    {
        std::cout << "[" << proc << "]\t" << val << std::endl;
    }
}

// Methods

void messager::reduce(int& val, MPI_Op op, MPI_Comm comm)
{
    int localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_INT, op, comm);
}

void messager::reduce(double& val, MPI_Op op, MPI_Comm comm)
{
    double localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_DOUBLE, op, comm);
}

void messager::maxReduce(int& val, MPI_Comm comm)
{
    int localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_INT, MPI_MAX, comm);
}

void messager::maxReduce(double& val, MPI_Comm comm)
{
    double localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_DOUBLE, MPI_MAX, comm);
}

void messager::minReduce(int& val, MPI_Comm comm)
{
    int localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_INT, MPI_MIN, comm);
}

void messager::minReduce(double& val, MPI_Comm comm)
{
    double localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_DOUBLE, MPI_MIN, comm);
}

void messager::sumReduce(int& val, MPI_Comm comm)
{
    int localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_INT, MPI_SUM, comm);
}

void messager::sumReduce(double& val, MPI_Comm comm)
{
    double localVal(val);
    MPI_Allreduce(&localVal, &val, 1, MPI_DOUBLE, MPI_SUM, comm);
}

void messager::sumReduce(std::vector<int>& arr, MPI_Comm comm)
{
    std::vector<int> localArr(arr);
    MPI_Allreduce(
        localArr.data(), arr.data(), arr.size(), MPI_INT, MPI_SUM, comm);
}

void messager::sumReduce(std::vector<double>& arr, MPI_Comm comm)
{
    std::vector<double> localArr(arr);
    MPI_Allreduce(
        localArr.data(), arr.data(), arr.size(), MPI_DOUBLE, MPI_SUM, comm);
}

void messager::broadcast(int& val, int root, MPI_Comm comm)
{
    MPI_Bcast(&val, 1, MPI_INT, root, comm);
}

void messager::broadcast(double& val, int root, MPI_Comm comm)
{
    MPI_Bcast(&val, 1, MPI_DOUBLE, root, comm);
}

void messager::broadcast(std::vector<int>& valArray, int root, MPI_Comm comm)
{
    int size;
    if (myProcNo() == root)
    {
        size = valArray.size();
    }
    broadcast(size, root, comm);
    if (myProcNo() != root)
    {
        valArray.resize(size);
    }
    MPI_Bcast(valArray.data(), valArray.size(), MPI_INT, root, comm);
}

void messager::broadcast(std::vector<double>& valArray, int root, MPI_Comm comm)
{
    int size;
    if (myProcNo() == root)
    {
        size = valArray.size();
    }
    broadcast(size, root, comm);
    if (myProcNo() != root)
    {
        valArray.resize(size);
    }
    MPI_Bcast(valArray.data(), valArray.size(), MPI_DOUBLE, root, comm);
}

void messager::broadcast(int* data, int count, int root, MPI_Comm comm)
{
    broadcast(count, root, comm);
    MPI_Bcast(data, count, MPI_INT, root, MPI_COMM_WORLD);
}

void messager::broadcast(double* data, int count, int root, MPI_Comm comm)
{
    broadcast(count, root, comm);
    MPI_Bcast(data, count, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

std::vector<int> messager::gatherValues(const int& localValue, MPI_Comm comm)
{
    std::vector<int> allValues;

    if (nProcs() > 1)
    {
        if (master(comm))
        {
            allValues.resize(nProcs());
        }

        MPI_Request request;
        if (MPI_Igather(&localValue,
                        1,
                        MPI_INT,
                        allValues.data(),
                        1,
                        MPI_INT,
                        0,
                        comm,
                        &request))
        {
            throw std::runtime_error("MPI_Igather failed.");
        }

        if (master(comm))
        {
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        // non-parallel: return own value
        allValues.resize(1);
        allValues[0] = localValue;
    }

    return allValues;
}

std::vector<double> messager::gatherValues(const double& localValue,
                                           MPI_Comm comm)
{
    std::vector<double> allValues;

    if (nProcs() > 1)
    {
        if (master(comm))
        {
            allValues.resize(nProcs());
        }

        MPI_Request request;
        if (MPI_Igather(&localValue,
                        1,
                        MPI_DOUBLE,
                        allValues.data(),
                        1,
                        MPI_DOUBLE,
                        0,
                        comm,
                        &request))
        {
            throw std::runtime_error("MPI_Igather failed.");
        }
    }
    else
    {
        // non-parallel: return own value
        allValues.resize(1);
        allValues[0] = localValue;
    }

    return allValues;
}

void messager::gatherVector(const std::vector<double>& vec,
                            std::vector<double>& vec_v)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather the sizes of each process's vector
    std::vector<int> sizes(size);
    int local_size = vec.size();
    MPI_Allgather(
        &local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate displacements for MPI_Allgatherv
    std::vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        displs[i] = displs[i - 1] + sizes[i - 1];
    }

    MPI_Allgatherv(vec.data(),
                   vec.size(),
                   MPI_DOUBLE,
                   vec_v.data(),
                   sizes.data(),
                   displs.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);
}

void messager::gatherVector(const std::vector<int>& vec,
                            std::vector<int>& vec_v)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather the sizes of each process's vector
    std::vector<int> sizes(size);
    int local_size = vec.size();
    MPI_Allgather(
        &local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate displacements for MPI_Allgatherv
    std::vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        displs[i] = displs[i - 1] + sizes[i - 1];
    }

    MPI_Allgatherv(vec.data(),
                   vec.size(),
                   MPI_INT,
                   vec_v.data(),
                   sizes.data(),
                   displs.data(),
                   MPI_INT,
                   MPI_COMM_WORLD);
}

void messager::gatherVector(const std::vector<uint64_t>& vec,
                            std::vector<uint64_t>& vec_v)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather the sizes of each process's vector
    std::vector<int> sizes(size);
    int local_size = vec.size();
    MPI_Allgather(
        &local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate displacements for MPI_Allgatherv
    std::vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        displs[i] = displs[i - 1] + sizes[i - 1];
    }

    MPI_Allgatherv(vec.data(),
                   vec.size(),
                   MPI_UINT64_T,
                   vec_v.data(),
                   sizes.data(),
                   displs.data(),
                   MPI_UINT64_T,
                   MPI_COMM_WORLD);
}

} // namespace accel
