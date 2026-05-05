// File       : nodeGraph.h
// Created    : Wed Jan 22 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: CRS node connectivity graph for linear solver sparsity patterns
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef NODEGRAPH_H
#define NODEGRAPH_H

// code
#include "CRSNodeGraph.h"

namespace accel
{

class mesh;

class nodeGraph : public ::linearSolver::CRSNodeGraph
{
public:
    nodeGraph(const MPI_Comm comm,
              mesh* meshPtr,
              const ::linearSolver::GraphLayout layout);

    // Access

    mesh& meshRef();

    const mesh& meshRef() const;

    // Methods

    void rebuildGraph();

private:
    mesh* meshPtr_;

    void buildGraph_() override;
};

// Out-of-line definitions

std::ostream& operator<<(std::ostream& os, const nodeGraph& graph);

} // namespace accel

#endif // nodeGraph
