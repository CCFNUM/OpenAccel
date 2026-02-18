// File       : AMG.h
// Created    : Fri Mar 14 2025 17:47:33 (+0100)
// Author     : Fabian Wermelinger
// Description: AMG instance factory
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef AMG_H
#define AMG_H

#include "CRSNodeGraph.h"
#include <yaml-cpp/yaml.h>

#if defined _WIN32 || defined __CYGWIN__
#define LIBAMG_EXPORT __declspec(dllexport)
#else
#if __GNUC__ >= 4
#define LIBAMG_EXPORT __attribute__((visibility("default")))
#else
#define LIBAMG_EXPORT
#endif
#endif

namespace linearSolver
{

LIBAMG_EXPORT void* getAMGSolverInstance(const size_t blocksize,
                                         const YAML::Node& node,
                                         const YAML::Node& solver_lookup,
                                         GraphLayout& layout);
LIBAMG_EXPORT void* getGMRESSolverInstance(const size_t blocksize,
                                           const YAML::Node& node,
                                           const YAML::Node& solver_lookup,
                                           GraphLayout& layout);
LIBAMG_EXPORT void* getDirectSolverInstance(const size_t blocksize,
                                            const YAML::Node& node,
                                            GraphLayout& layout);

} /* namespace linearSolver */

#endif /* AMG_H */
