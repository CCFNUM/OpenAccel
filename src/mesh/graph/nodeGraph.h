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
class zone;

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

    // Restrict this graph to a subset of mesh zones
    void setSubsetZones(const std::vector<const zone*>& zones)
    {
        subsetZones_ = zones;
    }

private:
    mesh* meshPtr_;

    // empty => full-mesh graph; non-empty => subset graph over these zones
    std::vector<const zone*> subsetZones_;

    void buildGraph_() override;
    void buildSubsetGraph_();
};

// Out-of-line definitions

std::ostream& operator<<(std::ostream& os, const nodeGraph& graph);

} // namespace accel

#endif // nodeGraph
