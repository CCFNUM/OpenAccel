// File       : nodeGraph.cpp
// Created    : Wed Jan 22 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "nodeGraph.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#include "ipInfo.h"
#include "mesh.h"
#include "messager.h"
#include "zone.h"

#include <algorithm>
#include <chrono>
#include <unordered_map>

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
    // subset graph (restricted to a zone vector) takes a dedicated path
    if (!subsetZones_.empty())
    {
        buildSubsetGraph_();
        return;
    }

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

    // Additional connections from non-conformal boundaries (fully-implicit
    // matrix)
    if (meshPtr_->hasInterfaces())
    {
        // conformal unions last (pass 1): row-merge needs equal pair stencils
        for (int pass = 0; pass < 2; ++pass)
            for (label iInterface = 0; iInterface < meshPtr_->nInterfaces();
                 iInterface++)
            {
                const auto& interf = meshPtr_->interfaceRef(iInterface);
                if (interf.isConformalTreatment() != (1 == pass))
                    continue;

                if (interf.isConformalTreatment())
                {
                    // loop over matching node pairs
                    const auto& nodePairs = interf.matchingNodePairVector();

                    // add src's element-neighbour columns to an owned row;
                    // mirror the base-loop stencil (reduced edge-pairs vs full)
                    // for symmetry
                    const bool reducedStencil =
                        (this->getLayout() & GraphLayout::Stencil__Reduced);
                    auto addElemNbrs =
                        [&](stk::mesh::Entity src, ulabel destRow)
                    {
                        const stk::mesh::Entity* elems =
                            bulkData.begin_elements(src);
                        const unsigned ne = bulkData.num_elements(src);
                        for (unsigned e = 0; e < ne; ++e)
                        {
                            const stk::mesh::Entity elem = elems[e];
                            const stk::mesh::Entity* en =
                                bulkData.begin_nodes(elem);
                            const unsigned nn = bulkData.num_nodes(elem);
                            if (!reducedStencil)
                            {
                                for (unsigned k = 0; k < nn; ++k)
                                    crsRowStencil[destRow].insert(en[k]);
                                continue;
                            }
                            // reduced: connect src only to its SCS-edge
                            // partners
                            MasterElement* meSCS =
                                MasterElementRepo::get_surface_master_element(
                                    bulkData.bucket(elem).topology());
                            const label numScsIp = meSCS->numIntPoints_;
                            const label* lrscv = meSCS->adjacentNodes();
                            for (label ip = 0; ip < numScsIp; ++ip)
                            {
                                const stk::mesh::Entity nodeL =
                                    en[lrscv[2 * ip]];
                                const stk::mesh::Entity nodeR =
                                    en[lrscv[2 * ip + 1]];
                                if (nodeL == src || nodeR == src)
                                {
                                    crsRowStencil[destRow].insert(nodeL);
                                    crsRowStencil[destRow].insert(nodeR);
                                }
                            }
                        }
                    };

                    const ulabel nOwned = static_cast<ulabel>(n_owned_nodes_);
                    for (const auto& nodePair : nodePairs)
                    {
                        const auto& node1 = nodePair.first;
                        const auto& node2 = nodePair.second;
                        const ulabel lid1 = bulkData.local_id(node1);
                        const ulabel lid2 = bulkData.local_id(node2);

                        // union each node's stencil into the other's owned row;
                        // co-owned uses prebuilt stencils, split rebuilds from
                        // elems
                        if (lid1 < nOwned)
                        {
                            if (lid2 < nOwned)
                                for (const auto& n : crsRowStencil[lid2])
                                    crsRowStencil[lid1].insert(n);
                            else
                                addElemNbrs(node2, lid1);
                        }
                        if (lid2 < nOwned)
                        {
                            if (lid1 < nOwned)
                                for (const auto& n : crsRowStencil[lid1])
                                    crsRowStencil[lid2].insert(n);
                            else
                                addElemNbrs(node1, lid2);
                        }
                    }
                }
                else
                {
                    // Note: In case of a reduced stencil, the node on the
                    // current side will always be implicitly connected to all
                    // nodes on the current and opposing sides

                    const auto addInterfaceStencil =
                        [&](const interfaceSideInfo& sideInfo)
                    {
                        for (const auto& faceIpInfoVec : sideInfo.ipInfoVec())
                        {
                            for (const ipInfo* ip : faceIpInfoVec)
                            {
                                if (ip->isExposed_)
                                    continue;

                                stk::mesh::Entity const*
                                    current_elem_node_rels =
                                        bulkData.begin_nodes(
                                            ip->currentElement_);
                                const label current_num_elem_nodes =
                                    bulkData.num_nodes(ip->currentElement_);

                                stk::mesh::Entity const*
                                    opposing_elem_node_rels =
                                        bulkData.begin_nodes(
                                            ip->opposingElement_);
                                const label opposing_num_elem_nodes =
                                    bulkData.num_nodes(ip->opposingElement_);

                                for (label nc = 0; nc < current_num_elem_nodes;
                                     ++nc)
                                {
                                    stk::mesh::EntityId rowLid =
                                        bulkData.local_id(
                                            current_elem_node_rels[nc]);
                                    if (rowLid >=
                                        static_cast<ulabel>(n_owned_nodes_))
                                        continue;

                                    // Cross-side coupling (every current x
                                    // opposing pair).
                                    for (label no = 0;
                                         no < opposing_num_elem_nodes;
                                         ++no)
                                    {
                                        crsRowStencil[rowLid].insert(
                                            opposing_elem_node_rels[no]);
                                    }

                                    // Same-side intra-element coupling.
                                    for (label nc2 = 0;
                                         nc2 < current_num_elem_nodes;
                                         ++nc2)
                                    {
                                        crsRowStencil[rowLid].insert(
                                            current_elem_node_rels[nc2]);
                                    }
                                }
                            }
                        }
                    };

                    addInterfaceStencil(interf.masterInfoRef());
                    addInterfaceStencil(interf.slaveInfoRef());
                }
            }

        // closure: pair rows holding a slave column also get its master; keep
        // pair sets equal
        for (label iI = 0; iI < meshPtr_->nInterfaces(); ++iI)
        {
            const interface& interf = meshPtr_->interfaceRef(iI);
            if (!interf.isConformalTreatment())
                continue;

            std::unordered_map<size_t, stk::mesh::Entity> slaveToMaster;
            for (const auto& nodePair : interf.matchingNodePairVector())
            {
                slaveToMaster.emplace(nodePair.second.local_offset(),
                                      nodePair.first);
            }

            const ulabel nOwned = static_cast<ulabel>(n_owned_nodes_);
            const auto closeRow = [&](stk::mesh::Entity n)
            {
                const ulabel lid = bulkData.local_id(n);
                if (lid >= nOwned)
                    return;
                const std::vector<stk::mesh::Entity> cols(
                    crsRowStencil[lid].begin(), crsRowStencil[lid].end());
                for (const stk::mesh::Entity c : cols)
                {
                    const auto it = slaveToMaster.find(c.local_offset());
                    if (it != slaveToMaster.end())
                        crsRowStencil[lid].insert(it->second);
                }
            };
            for (const auto& nodePair : interf.matchingNodePairVector())
            {
                closeRow(nodePair.first);
                closeRow(nodePair.second);
                // re-equalize the pair sets after closure insertions
                const ulabel lid1 = bulkData.local_id(nodePair.first);
                const ulabel lid2 = bulkData.local_id(nodePair.second);
                if (lid1 < nOwned && lid2 < nOwned)
                {
                    for (const auto& c : crsRowStencil[lid2])
                        crsRowStencil[lid1].insert(c);
                    for (const auto& c : crsRowStencil[lid1])
                        crsRowStencil[lid2].insert(c);
                }
            }
        }
    }

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

void nodeGraph::buildSubsetGraph_()
{
    const auto buildStart = std::chrono::steady_clock::now();

    stk::mesh::MetaData& metaData = meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = meshPtr_->bulkDataRef();

    const bool parallel = messager::parallel();

    // selector over the subset zones' interior parts (owned + ghosted)
    stk::mesh::PartVector interiorParts;
    for (const zone* z : subsetZones_)
    {
        const stk::mesh::PartVector zoneParts = z->interiorParts();
        interiorParts.insert(
            interiorParts.end(), zoneParts.begin(), zoneParts.end());
    }
    const stk::mesh::Selector subsetSel = stk::mesh::selectUnion(interiorParts);

    // remap mesh-local-id -> subset row (-1 = not in this sub-graph)
    local_to_row_.assign(meshPtr_->nNodes() + meshPtr_->nShadowNodes(),
                         static_cast<int64_t>(-1));

    // row -> entity (owned rows 0..M-1 first, ghost rows appended on demand)
    std::vector<stk::mesh::Entity> rowEntity;

    // 1.) contiguous rows for owned subset nodes, walked in local-id order so
    // the rows inherit the mesh node ordering (bucket order would undo it)
    {
        const stk::mesh::Selector ownedSubsetSel =
            metaData.locally_owned_part() & subsetSel;
        const auto& localNodeIDToEntity = meshPtr_->localNodeIDToEntity();
        for (label lid = 0; lid < meshPtr_->nNodes(); ++lid)
        {
            const stk::mesh::Entity n = localNodeIDToEntity[lid];
            if (!ownedSubsetSel(bulkData.bucket(n)))
                continue;

            local_to_row_[lid] = static_cast<int64_t>(rowEntity.size());
            rowEntity.push_back(n);
        }
    }
    n_owned_nodes_ = static_cast<Index>(rowEntity.size());

    // pre-mark non-owned subset nodes (-2): no selector calls in the hot path
    {
        const stk::mesh::Selector ghostSubsetSel =
            subsetSel & !metaData.locally_owned_part();
        for (auto* bp :
             bulkData.get_buckets(stk::topology::NODE_RANK, ghostSubsetSel))
        {
            for (size_t i = 0; i < bp->size(); ++i)
            {
                const ulabel lid = bulkData.local_id((*bp)[i]);
                if (lid < local_to_row_.size())
                    local_to_row_[lid] = static_cast<int64_t>(-2);
            }
        }
    }

    // ghost rows assigned on demand for subset nodes owned by other ranks
    const auto rowOf = [&](stk::mesh::Entity n) -> int64_t
    {
        const ulabel lid = bulkData.local_id(n);
        // out-of-range lid (DG-ghosted solid node) not in sub-graph
        if (lid >= local_to_row_.size())
            return -1;
        const int64_t r = local_to_row_[lid];
        if (r >= 0)
            return r;
        if (r == -2) // pre-marked subset ghost: assign the next ghost row
        {
            const int64_t gr = static_cast<int64_t>(rowEntity.size());
            local_to_row_[lid] = gr;
            rowEntity.push_back(n);
            return gr;
        }
        return -1; // outside subset (e.g. solid side of fluid-solid) -> dropped
    };

    // stencil over owned rows; columns may be owned or ghost subset nodes
    // (sorted-unique vectors: far cheaper than one tree node per column)
    std::vector<std::vector<Index>> crsRowStencil(n_owned_nodes_);
    for (auto& row : crsRowStencil)
        row.reserve(32);
    const auto insertSorted = [](std::vector<Index>& row, const Index c)
    {
        const auto it = std::lower_bound(row.begin(), row.end(), c);
        if (it == row.end() || *it != c)
            row.insert(it, c);
    };
    const auto addPair = [&](stk::mesh::Entity a, stk::mesh::Entity b)
    {
        const ulabel la = bulkData.local_id(a);
        if (la >= local_to_row_.size())
            return;
        const int64_t r1 = local_to_row_[la];
        if (r1 < 0 || r1 >= static_cast<int64_t>(n_owned_nodes_))
            return; // a is not an owned subset row
        const int64_t r2 = rowOf(b);
        if (r2 >= 0)
            insertSorted(crsRowStencil[r1], static_cast<Index>(r2));
    };

    // intra-element coupling (owned + ghosted subset elements)
    {
        const bool reducedStencil =
            (this->getLayout() & GraphLayout::Stencil__Reduced);
        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, subsetSel);
        for (auto* ebp : elementBuckets)
        {
            // reduced: connect only the SCS edge partners (mirrors buildGraph_)
            MasterElement* meSCS =
                reducedStencil ? MasterElementRepo::get_surface_master_element(
                                     ebp->topology())
                               : nullptr;
            const label numScsIp = reducedStencil ? meSCS->numIntPoints_ : 0;
            const label* lrscv = reducedStencil ? meSCS->adjacentNodes() : 0;
            for (size_t ie = 0; ie < ebp->size(); ++ie)
            {
                stk::mesh::Entity const* nr = ebp->begin_nodes(ie);
                if (reducedStencil)
                {
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        const stk::mesh::Entity nodeL = nr[lrscv[2 * ip]];
                        const stk::mesh::Entity nodeR = nr[lrscv[2 * ip + 1]];
                        addPair(nodeL, nodeL);
                        addPair(nodeL, nodeR);
                        addPair(nodeR, nodeL);
                        addPair(nodeR, nodeR);
                    }
                    continue;
                }
                const label nn = bulkData.num_nodes((*ebp)[ie]);
                for (label a = 0; a < nn; ++a)
                {
                    for (label b = 0; b < nn; ++b)
                    {
                        addPair(nr[a], nr[b]);
                    }
                }
            }
        }
    }

    // interface coupling within the subset (pairs outside the subset dropped)
    if (meshPtr_->hasInterfaces())
    {
        // conformal unions last (pass 1): row-merge needs equal pair stencils
        for (int pass = 0; pass < 2; ++pass)
            for (label iI = 0; iI < meshPtr_->nInterfaces(); ++iI)
            {
                const interface& interf = meshPtr_->interfaceRef(iI);
                if (interf.isConformalTreatment() != (1 == pass))
                    continue;
                if (interf.isConformalTreatment())
                {
                    // row-merge: union pair stencils (mirrors full-graph build)
                    const bool reducedStencil =
                        (this->getLayout() & GraphLayout::Stencil__Reduced);
                    const auto addColToRow =
                        [&](int64_t r, stk::mesh::Entity col)
                    {
                        const int64_t rc = rowOf(col);
                        if (rc >= 0)
                            insertSorted(crsRowStencil[r],
                                         static_cast<Index>(rc));
                    };
                    const auto addElemNbrs =
                        [&](stk::mesh::Entity src, int64_t r)
                    {
                        const stk::mesh::Entity* elems =
                            bulkData.begin_elements(src);
                        const unsigned ne = bulkData.num_elements(src);
                        for (unsigned e = 0; e < ne; ++e)
                        {
                            const stk::mesh::Entity elem = elems[e];
                            const stk::mesh::Entity* en =
                                bulkData.begin_nodes(elem);
                            const unsigned nn = bulkData.num_nodes(elem);
                            if (!reducedStencil)
                            {
                                for (unsigned k = 0; k < nn; ++k)
                                    addColToRow(r, en[k]);
                                continue;
                            }
                            // reduced: connect src only to its SCS-edge
                            // partners
                            MasterElement* meSCS =
                                MasterElementRepo::get_surface_master_element(
                                    bulkData.bucket(elem).topology());
                            const label numScsIp = meSCS->numIntPoints_;
                            const label* lrscv = meSCS->adjacentNodes();
                            for (label ip = 0; ip < numScsIp; ++ip)
                            {
                                const stk::mesh::Entity nodeL =
                                    en[lrscv[2 * ip]];
                                const stk::mesh::Entity nodeR =
                                    en[lrscv[2 * ip + 1]];
                                if (nodeL == src || nodeR == src)
                                {
                                    addColToRow(r, nodeL);
                                    addColToRow(r, nodeR);
                                }
                            }
                        }
                    };

                    const int64_t nOwned = static_cast<int64_t>(n_owned_nodes_);
                    for (const auto& nodePair : interf.matchingNodePairVector())
                    {
                        const auto& node1 = nodePair.first;
                        const auto& node2 = nodePair.second;
                        const int64_t r1 = rowOf(node1);
                        const int64_t r2 = rowOf(node2);
                        const bool r1Owned = (r1 >= 0 && r1 < nOwned);
                        const bool r2Owned = (r2 >= 0 && r2 < nOwned);
                        if (r1Owned)
                        {
                            if (r2Owned)
                                for (const Index c : crsRowStencil[r2])
                                    insertSorted(crsRowStencil[r1], c);
                            else
                                addElemNbrs(node2, r1);
                        }
                        if (r2Owned)
                        {
                            if (r1Owned)
                                for (const Index c : crsRowStencil[r1])
                                    insertSorted(crsRowStencil[r2], c);
                            else
                                addElemNbrs(node1, r2);
                        }
                    }
                    continue;
                }

                const auto addSide = [&](const interfaceSideInfo& side)
                {
                    for (const auto& faceIpInfoVec : side.ipInfoVec())
                    {
                        for (const ipInfo* ip : faceIpInfoVec)
                        {
                            if (ip->isExposed_)
                                continue;
                            stk::mesh::Entity const* cn =
                                bulkData.begin_nodes(ip->currentElement_);
                            const label ncur =
                                bulkData.num_nodes(ip->currentElement_);
                            stk::mesh::Entity const* on =
                                bulkData.begin_nodes(ip->opposingElement_);
                            const label nopp =
                                bulkData.num_nodes(ip->opposingElement_);
                            for (label a = 0; a < ncur; ++a)
                            {
                                for (label b = 0; b < nopp; ++b)
                                {
                                    addPair(cn[a], on[b]);
                                }
                                for (label b = 0; b < ncur; ++b)
                                {
                                    addPair(cn[a], cn[b]);
                                }
                            }
                        }
                    }
                };
                addSide(interf.masterInfoRef());
                addSide(interf.slaveInfoRef());
            }

        // closure: pair rows holding a slave column also get its master; keep
        // pair sets equal
        for (label iI = 0; iI < meshPtr_->nInterfaces(); ++iI)
        {
            const interface& interf = meshPtr_->interfaceRef(iI);
            if (!interf.isConformalTreatment())
                continue;

            std::unordered_map<Index, Index> slaveToMasterRow;
            for (const auto& nodePair : interf.matchingNodePairVector())
            {
                const int64_t r1 = rowOf(nodePair.first);
                const int64_t r2 = rowOf(nodePair.second);
                if (r1 >= 0 && r2 >= 0)
                    slaveToMasterRow.emplace(static_cast<Index>(r2),
                                             static_cast<Index>(r1));
            }

            const int64_t nOwned = static_cast<int64_t>(n_owned_nodes_);
            const auto closeRow = [&](const int64_t r)
            {
                if (r < 0 || r >= nOwned)
                    return;
                const std::vector<Index> cols = crsRowStencil[r]; // copy
                for (const Index c : cols)
                {
                    const auto it = slaveToMasterRow.find(c);
                    if (it != slaveToMasterRow.end())
                        insertSorted(crsRowStencil[r], it->second);
                }
            };
            for (const auto& nodePair : interf.matchingNodePairVector())
            {
                const int64_t r1 = rowOf(nodePair.first);
                const int64_t r2 = rowOf(nodePair.second);
                closeRow(r1);
                closeRow(r2);
                // re-equalize the pair sets after closure insertions
                if (r1 >= 0 && r1 < nOwned && r2 >= 0 && r2 < nOwned)
                {
                    for (const Index c : crsRowStencil[r2])
                        insertSorted(crsRowStencil[r1], c);
                    for (const Index c : crsRowStencil[r1])
                        insertSorted(crsRowStencil[r2], c);
                }
            }
        }
    }

    // untouched ghost marks revert to -1 (localToRow contract: row or -1)
    for (auto& r : local_to_row_)
    {
        if (r == -2)
            r = static_cast<int64_t>(-1);
    }

    n_ghost_nodes_ = static_cast<Index>(rowEntity.size()) - n_owned_nodes_;

    // contiguous-per-rank global subset ids: offset = exclusive prefix sum of M
    Index offset = 0;
    if (parallel)
    {
        Index m = n_owned_nodes_;
        MPI_Exscan(&m, &offset, 1, MPIDataType<Index>::type(), MPI_SUM, comm_);
    }

    // owned id = offset+row; ghost id from its owner via the aux field exchange
    local_to_global_.assign(local_to_row_.size(), static_cast<int64_t>(-1));
    stk::mesh::Field<label>* auxField =
        parallel
            ? metaData.get_field<label>(stk::topology::NODE_RANK, mesh::aux)
            : nullptr;
    if (parallel)
    {
        if (!auxField)
            errorMsg("nodeGraph::buildSubsetGraph_: global-id field 'aux' not "
                     "found; cannot number subset ghosts in parallel");
        for (Index r = 0; r < n_owned_nodes_; ++r)
        {
            label* fd = stk::mesh::field_data(*auxField, rowEntity[r]);
            if (!fd)
                errorMsg(
                    "nodeGraph::buildSubsetGraph_: aux not allocated on an "
                    "owned subset node");
            *fd = static_cast<label>(offset + r);
        }
        stk::mesh::communicate_field_data(bulkData, {auxField});
    }
    for (Index r = 0; r < static_cast<Index>(rowEntity.size()); ++r)
    {
        const stk::mesh::Entity n = rowEntity[r];
        int64_t gid;
        if (r < n_owned_nodes_)
        {
            gid = static_cast<int64_t>(offset) + r;
        }
        else
        {
            const label* fd = stk::mesh::field_data(*auxField, n);
            if (!fd)
                errorMsg("nodeGraph::buildSubsetGraph_: aux not allocated on a "
                         "subset ghost node");
            gid = static_cast<int64_t>(*fd);
        }
        local_to_global_[bulkData.local_id(n)] = gid;
    }

    // finalize CRS: local order primary=row/secondary=gid; global order swaps
    const bool localOrder = this->isLocalColumnOrder();
    row_ptr_.resize(n_owned_nodes_ + 1);
    row_ptr_[0] = 0;
    Index nnz = 0;
    for (Index i = 0; i < n_owned_nodes_; ++i)
    {
        nnz += static_cast<Index>(crsRowStencil[i].size());
        row_ptr_[i + 1] = nnz;
    }
    primary_indices_.resize(nnz);
    secondary_indices_.resize(nnz);
    Index k = 0;
    for (Index i = 0; i < n_owned_nodes_; ++i)
    {
        for (const Index col : crsRowStencil[i])
        {
            const Index globalcol = static_cast<Index>(
                local_to_global_[bulkData.local_id(rowEntity[col])]);
            primary_indices_[k] = localOrder ? col : globalcol;
            secondary_indices_[k] = localOrder ? globalcol : col;
            ++k;
        }
    }

    // restore aux to full global ids (for updateLocalNodeIDs_ on moving mesh);
    // only subset-zone copies were overwritten (owners + exchanged ghosts)
    if (parallel)
    {
        const stk::mesh::BucketVector& restoreBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, subsetSel);
        for (auto* bp : restoreBuckets)
        {
            label* a = stk::mesh::field_data(*auxField, *bp);
            for (size_t i = 0; i < bp->size(); ++i)
            {
                a[i] = static_cast<label>(bulkData.global_id((*bp)[i]));
            }
        }
    }

    // local-order columns come out of the sorted stencils already in order
    if (!this->isLocalColumnOrder())
        this->sortPrimaryIndices_();

    if (messager::master())
    {
        const auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - buildStart);
        std::cout << "  [subset-graph] build: " << n_owned_nodes_ << " rows, "
                  << nnz << " nnz, " << dt.count() << " ms" << std::endl;
    }
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
