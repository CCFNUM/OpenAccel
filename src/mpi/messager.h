// File : messager.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: MPI communication wrapper for reductions, broadcasts, and
// gathers
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef MESSAGER_H
#define MESSAGER_H

// Built-in libraries
#include <cstdint>
#include <mpi.h>
#include <string>
#include <vector>

namespace accel
{

class messager
{
private:
    static MPI_Comm comm_;

public:
    static MPI_Comm comm(MPI_Comm comm = comm_);

    static bool parallel(MPI_Comm comm = comm_);

    static bool master(MPI_Comm comm = comm_);

    static int myProcNo(MPI_Comm comm = comm_);

    static int nProcs(MPI_Comm comm = comm_);

    static void barrier(MPI_Comm comm = comm_);

    static void exitAll(MPI_Comm comm = comm_);

    static void markProcReached(MPI_Comm comm = comm_);

    static void markProcReachedAndExit(MPI_Comm comm = comm_);

    static void print(std::string val, MPI_Comm comm = comm_);

    static void print(int val, MPI_Comm comm = comm_);

    static void print(double val, MPI_Comm comm = comm_);

    static void procPrint(int val, int proc, MPI_Comm comm = comm_);

    static void procPrint(double val, int proc, MPI_Comm comm = comm_);

    // Methods

    static void reduce(int& val, MPI_Op op, MPI_Comm comm = comm_);

    static void reduce(double& val, MPI_Op op, MPI_Comm comm = comm_);

    static void maxReduce(int& val, MPI_Comm comm = comm_);

    static void maxReduce(double& val, MPI_Comm comm = comm_);

    static void minReduce(int& val, MPI_Comm comm = comm_);

    static void minReduce(double& val, MPI_Comm comm = comm_);

    static void sumReduce(int& val, MPI_Comm comm = comm_);

    static void sumReduce(double& val, MPI_Comm comm = comm_);

    static void sumReduce(std::vector<int>& arr, MPI_Comm comm = comm_);

    static void sumReduce(std::vector<double>& arr, MPI_Comm comm = comm_);

    static void broadcast(int& val, int root, MPI_Comm comm = comm_);

    static void broadcast(double& val, int root, MPI_Comm comm = comm_);

    static void
    broadcast(std::vector<int>& valArray, int root, MPI_Comm comm = comm_);

    static void
    broadcast(std::vector<double>& valArray, int root, MPI_Comm comm = comm_);

    static void
    broadcast(int* data, int count, int root, MPI_Comm comm = comm_);

    static void
    broadcast(double* data, int count, int root, MPI_Comm comm = comm_);

    static std::vector<int> gatherValues(const int& localValue,
                                         MPI_Comm comm = comm_);

    static std::vector<double> gatherValues(const double& localValue,
                                            MPI_Comm comm = comm_);

    static void gatherVector(const std::vector<double>& vec,
                             std::vector<double>& vec_v);

    static void gatherVector(const std::vector<int>& vec,
                             std::vector<int>& vec_v);

    static void gatherVector(const std::vector<uint64_t>& vec,
                             std::vector<uint64_t>& vec_v);
};

} // namespace accel

#endif // MESSAGER_H
