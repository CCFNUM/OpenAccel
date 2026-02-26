// File       : nodeGraph.cpp
// Created    : Wed Jan 22 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "nodeGraph.h"
#ifdef HAS_INTERFACE
#include "dgInfo.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#endif /* HAS_INTERFACE */
#include "mesh.h"
#include "messager.h"

using namespace linearSolver;

namespace accel
{

nodeGraph::nodeGraph(const MPI_Comm comm,
                     mesh* meshPtr,
                     const GraphLayout layout)
    : CRSNodeGraph(comm, layout), meshPtr_(meshPtr)
{
}

void nodeGraph::buildGraph_()
{
    stk::mesh::MetaData& metaData = meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = meshPtr_->bulkDataRef();

    n_owned_nodes_ = meshPtr_->nNodes();
    n_ghost_nodes_ = meshPtr_->nShadowNodes();
    std::vector<std::set<stk::mesh::Entity>> crsRowStencil(n_owned_nodes_);

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK,
                             meshPtr_->universalInteriorPartsSelector());
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        for (size_t iElement = 0; iElement < nElementsPerBucket; ++iElement)
        {
            // get elem
            const stk::mesh::Entity& elem = elementBucket[iElement];
            stk::mesh::Entity const* elemNodeRels =
                elementBucket.begin_nodes(iElement);

            // In case of a reduced stencil, only the neighbouring nodes that
            // share an ip with the central node are considered. This will
            // significantly cut down the required matrix storage
            if (this->getLayout() & GraphLayout::Stencil__Reduced)
            {
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = elemNodeRels[il];
                    stk::mesh::Entity nodeR = elemNodeRels[ir];

                    ulabel idL = bulkData.local_id(nodeL);
                    ulabel idR = bulkData.local_id(nodeR);

                    if (idL < static_cast<ulabel>(
                                  n_owned_nodes_)) // implies this is a local
                    // node in case of parallel
                    {
                        crsRowStencil[idL].insert(nodeL);
                        crsRowStencil[idL].insert(nodeR);
                    }

                    if (idR < static_cast<ulabel>(
                                  n_owned_nodes_)) // implies this is a local
                    // node in case of parallel
                    {
                        crsRowStencil[idR].insert(nodeL);
                        crsRowStencil[idR].insert(nodeR);
                    }
                }
            }
            else
            {
                label numNodes = bulkData.num_nodes(elem);

                // sanity check on num nodes
                STK_ThrowAssert(numNodes == nodesPerElement);

                for (label ni1 = 0; ni1 < numNodes; ++ni1)
                {
                    ulabel id1 = bulkData.local_id(elemNodeRels[ni1]);

                    if (id1 < static_cast<ulabel>(
                                  n_owned_nodes_)) // implies this is a local
                    // node in case of parallel
                    {
                        for (label ni2 = 0; ni2 < numNodes; ++ni2)
                        {
                            crsRowStencil[id1].insert(elemNodeRels[ni2]);
                        }
                    }
                }
            }
        }
    }

#ifdef HAS_INTERFACE
    // Additional connections from non-conformal boundaries (fully-implicit
    // matrix)
    if (meshPtr_->hasInterfaces())
    {
        for (label iInterface = 0; iInterface < meshPtr_->nInterfaces();
             iInterface++)
        {
            const auto& interf = meshPtr_->interfaceRef(iInterface);

            if (interf.isConformalTreatment())
            {
                // loop over matching node pairs
                const auto& nodePairs = interf.matchingNodePairVector();

                for (const auto& nodePair : nodePairs)
                {
                    // data for stencil 1 (of node 1)
                    const auto& node1 = nodePair.first;
                    const label& lid1 = bulkData.local_id(node1);

                    // data for stencil 2 (of node 2)
                    const auto& node2 = nodePair.second;
                    const label& lid2 = bulkData.local_id(node2);

                    // add stencil of node 2 to that of node 1
                    for (const auto& node : crsRowStencil[lid2])
                    {
                        crsRowStencil[lid1].insert(node);
                    }

                    // add stencil of node 1 to that of node 2
                    for (const auto& node : crsRowStencil[lid1])
                    {
                        crsRowStencil[lid2].insert(node);
                    }

                    // diagnostics
                    assert(crsRowStencil[lid1] == crsRowStencil[lid2]);
                }
            }
            else
            {
                // Note: In case of a reduced stencil, the node on the current
                // side will always be implicitly connected to all nodes on the
                // opposing side

                // Master
                {
                    const auto& masterInterface = interf.masterInfoRef();

                    // extract vector of dgInfo
                    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
                        masterInterface.dgInfoVec_;

                    for (label iSide = 0;
                         iSide < static_cast<label>(dgInfoVec.size());
                         iSide++)
                    {
                        const std::vector<dgInfo*>& faceDgInfoVec =
                            dgInfoVec[iSide];

                        // now loop over all the DgInfo objects on this
                        // particular exposed face
                        for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
                        {
                            dgInfo* dgInfo = faceDgInfoVec[k];

                            if (dgInfo->gaussPointExposed_)
                                continue;

                            // extract current/opposing face/element
                            stk::mesh::Entity currentFace =
                                dgInfo->currentFace_;
                            stk::mesh::Entity opposingElement =
                                dgInfo->opposingElement_;

                            // master element
                            MasterElement* meFCCurrent = dgInfo->meFCCurrent_;

                            // local ip, ordinals, etc
                            const label currentGaussPointId =
                                dgInfo->currentGaussPointId_;

                            // mapping from ip to nodes for this ordinal
                            const label* faceIpNodeMap =
                                meFCCurrent->ipNodeMap();
                            stk::mesh::Entity const* current_face_node_rels =
                                bulkData.begin_nodes(currentFace);
                            const label nn = faceIpNodeMap[currentGaussPointId];
                            stk::mesh::Entity node = current_face_node_rels[nn];

                            stk::mesh::EntityId nearestID =
                                bulkData.local_id(node);

                            if (nearestID <
                                static_cast<ulabel>(
                                    n_owned_nodes_)) // implies this is a local
                                                     // node in case of parallel
                            {
                                // gather opposing face data
                                stk::mesh::Entity const*
                                    opposing_elem_node_rels =
                                        bulkData.begin_nodes(opposingElement);
                                const label opposing_num_elem_nodes =
                                    bulkData.num_nodes(opposingElement);
                                for (label ni = 0; ni < opposing_num_elem_nodes;
                                     ++ni)
                                {
                                    crsRowStencil[nearestID].insert(
                                        opposing_elem_node_rels[ni]);
                                }
                            }
                        }
                    }
                }

                // Slave
                {
                    const auto& slaveInterface = interf.slaveInfoRef();

                    // extract vector of dgInfo
                    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
                        slaveInterface.dgInfoVec_;

                    for (label iSide = 0;
                         iSide < static_cast<label>(dgInfoVec.size());
                         iSide++)
                    {
                        const std::vector<dgInfo*>& faceDgInfoVec =
                            dgInfoVec[iSide];

                        // now loop over all the DgInfo objects on this
                        // particular exposed face
                        for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
                        {
                            dgInfo* dgInfo = faceDgInfoVec[k];

                            if (dgInfo->gaussPointExposed_)
                                continue;

                            // extract current/opposing face/element
                            stk::mesh::Entity currentFace =
                                dgInfo->currentFace_;
                            stk::mesh::Entity opposingElement =
                                dgInfo->opposingElement_;

                            // master element
                            MasterElement* meFCCurrent = dgInfo->meFCCurrent_;

                            // local ip, ordinals, etc
                            const label currentGaussPointId =
                                dgInfo->currentGaussPointId_;

                            // mapping from ip to nodes for this ordinal
                            const label* faceIpNodeMap =
                                meFCCurrent->ipNodeMap();
                            stk::mesh::Entity const* current_face_node_rels =
                                bulkData.begin_nodes(currentFace);
                            const label nn = faceIpNodeMap[currentGaussPointId];
                            stk::mesh::Entity node = current_face_node_rels[nn];

                            stk::mesh::EntityId nearestID =
                                bulkData.local_id(node);

                            if (nearestID <
                                static_cast<ulabel>(
                                    n_owned_nodes_)) // implies this is a local
                                                     // node in case of parallel
                            {
                                // gather opposing face data
                                stk::mesh::Entity const*
                                    opposing_elem_node_rels =
                                        bulkData.begin_nodes(opposingElement);
                                const label opposing_num_elem_nodes =
                                    bulkData.num_nodes(opposingElement);
                                for (label ni = 0; ni < opposing_num_elem_nodes;
                                     ++ni)
                                {
                                    crsRowStencil[nearestID].insert(
                                        opposing_elem_node_rels[ni]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif /* HAS_INTERFACE */

    // create CRS structure
    row_ptr_.resize(n_owned_nodes_ + 1);

    // Get number of non-zero nodes and populate CRS arrays
    row_ptr_[0] = 0;
    size_t nonZeroNodes = 0;
    ulabel elementPerRow = 0;
    for (Index i = 0; i < n_owned_nodes_; ++i)
    {
        auto& comp = crsRowStencil[i];
        for (auto it = comp.begin(); it != comp.end(); ++it)
        {
            nonZeroNodes++;
            elementPerRow++;
        }
        row_ptr_[i + 1] = elementPerRow;
    }

    primary_indices_.resize(nonZeroNodes);
    secondary_indices_.resize(nonZeroNodes);

#ifndef NDEBUG
    const Index __n_active_nodes =
        meshPtr_->nNodes() + meshPtr_->nShadowNodes();

    Index __n_global_nodes = 0;
    MPI_Allreduce(&n_owned_nodes_,
                  &__n_global_nodes,
                  1,
                  MPIDataType<Index>::type(),
                  MPI_SUM,
                  messager::comm());
#endif /* NDEBUG */

    Index k = 0;
    if (this->isLocalColumnOrder())
    {
        for (Index i = 0; i < n_owned_nodes_; ++i)
        {
            const auto& comp = crsRowStencil[i];
            for (auto it = comp.begin(); it != comp.end(); ++it)
            {
                const stk::mesh::Entity node = (*it);
                primary_indices_[k] = bulkData.local_id(node);
                secondary_indices_[k] = bulkData.global_id(node);
                assert(0 <= primary_indices_[k] &&
                       primary_indices_[k] < __n_active_nodes);
                assert(0 <= secondary_indices_[k] &&
                       secondary_indices_[k] < __n_global_nodes);
                k++;
            }
        }
    }
    else if (this->isGlobalColumnOrder())
    {
        for (Index i = 0; i < n_owned_nodes_; ++i)
        {
            const auto& comp = crsRowStencil[i];
            for (auto it = comp.begin(); it != comp.end(); ++it)
            {
                const stk::mesh::Entity node = (*it);
                primary_indices_[k] = bulkData.global_id(node);
                secondary_indices_[k] = bulkData.local_id(node);
                assert(0 <= primary_indices_[k] &&
                       primary_indices_[k] < __n_global_nodes);
                assert(0 <= secondary_indices_[k] &&
                       secondary_indices_[k] < __n_active_nodes);
                k++;
            }
        }
    }
    else
    {
        errorMsg("nodeGraph::buildGraph_: only matrices with local or "
                 "global column index order are supported");
    }

    this->sortPrimaryIndices_();

#ifndef NDEBUG
    {
        // graph sanity check:
        assert(n_owned_nodes_ == meshPtr_->nNodes());
        assert(n_ghost_nodes_ == meshPtr_->nShadowNodes());
        assert(__n_active_nodes == meshPtr_->nActiveNodes());

        // local column order in mesh API:
        // | owned nodes | shadow nodes | useless nodes (not part of graph) |
        // 0 < nNodes < nActiveNodes < nAllNodes

        // CRSNodeGraph API:
        // | owned nodes | ghost nodes |
        // 0 < nOwnedNodes < nAllNodes
        // nAllNodes in graph == owned nodes + ghost nodes (==shadow nodes)

        // a.) local index consistency
        std::set<Index> local_idx;
        std::set<Index> global_idx;
        for (const Index i : this->localIndices())
        {
            local_idx.insert(i);
        }
        for (const Index i : this->globalIndices())
        {
            global_idx.insert(i);
        }

        std::vector<Index> missing_owned, missing_ghost;
        for (Index i = 0; i < n_owned_nodes_; i++) // owned nodes
        {
            if (local_idx.erase(i) == 0)
            {
                missing_owned.push_back(i);
            }
        }
        for (Index i = n_owned_nodes_; i < __n_active_nodes; i++) // ghosts
        {
            if (local_idx.erase(i) == 0)
            {
                missing_ghost.push_back(i);
            }
        }
        assert(missing_owned.empty()); // required

        // may not be satisfied due to full/reduced stencils and if only a
        // subset of nodes in a particular element are used in the graph.
        for (int rank = 0; rank < messager::nProcs(); rank++)
        {
            if (rank == messager::myProcNo())
            {
                if (!missing_ghost.empty())
                {
                    std::cout << "WARNING (rank=" << rank
                              << "): there are active nodes in STK mesh not "
                                 "used in node graph (local node indices:";
                    for (const Index i : missing_ghost)
                    {
                        std::cout << " " << i;
                    }
                    std::cout << ")" << std::endl;
                }
            }
            MPI_Barrier(messager::comm());
        }

        // b.) any local column index must be less than total active nodes in
        // layout (if this fails it means the graph contains local id's that
        // belong to useless nodes).
        for (const label col_idx : this->localIndices())
        {
            assert(col_idx < __n_active_nodes);
        }
    }
#endif /* NDEBUG */
}

void nodeGraph::rebuildGraph()
{
    // sanity check: typically you want to rebuild something that has been built
    // previously. The assertion itself would not be necessary however.
    assert(this->isBuilt());

    // reset first
    resetGraph_();

    // build again
    buildGraph();
}

// Access

mesh& nodeGraph::meshRef()
{
    return *meshPtr_;
}

const mesh& nodeGraph::meshRef() const
{
    return *meshPtr_;
}

std::ostream& operator<<(std::ostream& os, const nodeGraph& graph)
{
    if (messager::master())
    {
        os << std::endl << "Node graph" << std::endl;
    }

    const stk::mesh::MetaData& metaData = graph.meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData = graph.meshRef().bulkDataRef();

    STKScalarField* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, graph.meshRef().getCoordinateFieldName());

    const auto& offsets = graph.offsets();
    const auto& indices = graph.indices();

    for (label iProc = 0; iProc < messager::nProcs(); iProc++)
    {
        if (messager::myProcNo() == iProc)
        {
            if (messager::parallel())
            {
                os << "Proc: " << iProc << std::endl;
            }
            os << "{" << std::endl;

            // select all locally owned nodes relevant to the field;
            stk::mesh::BucketVector const& nodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK,
                graph.meshRef().locallyOwnedInteriorPartsSelector());

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& b = **ib;
                const stk::mesh::Bucket::size_type length = b.size();
                scalar* value = stk::mesh::field_data(*coordsSTKFieldPtr, b);
                for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
                {
                    stk::mesh::Entity node = b[k];

                    stk::mesh::EntityId lid = bulkData.local_id(node);
                    stk::mesh::EntityId id = bulkData.identifier(node);

                    std::cout << "lid: " << lid << " [id=" << id
                              << "], coords: (";

                    for (label i = 0; i < 2; i++)
                    {
                        os << std::scientific << std::setprecision(14)
                           << value[SPATIAL_DIM * k + i] << ", ";
                    }
                    os << std::scientific << std::setprecision(14)
                       << value[SPATIAL_DIM * k + (SPATIAL_DIM - 1)];

                    std::span<const label> row =
                        std::span<const label>(indices).subspan(
                            offsets[lid], offsets[lid + 1] - offsets[lid]);

                    std::cout << "), connected to: " << row.size() << "( ";

                    for (auto colIdx : row)
                    {
                        std::cout << colIdx << " ";
                    }

                    std::cout << ")" << std::endl;
                }
            }
            os << "}" << std::endl;
        }
        messager::barrier();
    }

    return os;
}

} // namespace accel
