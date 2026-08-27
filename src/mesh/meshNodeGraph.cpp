// File       : meshNodeGraph.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "controls.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#include "ipInfo.h"
#include "mesh.h"
#include "messager.h"
#include "zone.h"

namespace accel
{

void mesh::lazyInitializeNodeGraph_()
{
    // For consistency, node local and global id's must be adjusted
    // to align with CRS needs
    initializeLocalNodeIDs_();

    // node graph (connectivity) is built on-demand, see
    // mesh::getGlobalOrderGraphPtr and mesh:getLocalOrderGraphPtr.
}

std::unique_ptr<nodeGraph>
mesh::createNodeGraph_(const ::linearSolver::GraphLayout layout)
{
    return std::make_unique<nodeGraph>(
        messager::comm(), this, layout | stencil_);
}

std::vector<label> mesh::zoneKey_(const std::vector<const zone*>& zones) const
{
    std::vector<label> key;
    key.reserve(zones.size());
    for (const zone* z : zones)
        key.push_back(z->index());
    std::sort(key.begin(), key.end());
    key.erase(std::unique(key.begin(), key.end()), key.end());
    return key;
}

bool mesh::useFullGraphForZones_(const std::vector<const zone*>& zones) const
{
    // expert override: disable subset node graphs entirely
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.forceFullNodeGraph_)
        return true;

    const std::vector<label> key = zoneKey_(zones);
    // empty, or spanning every zone => the full mesh graph already fits
    if (key.empty() || key.size() >= zoneVector_.size())
        return true;

    // row-merge cannot couple across the subset boundary: mixed -> full graph
    const auto inKey = [&key](const label z)
    { return std::binary_search(key.begin(), key.end(), z); };
    for (label i = 0; i < nInterfaces(); ++i)
    {
        const interface& interf = interfaceRef(i);
        if (!interf.isConformalTreatment())
            continue;
        if (inKey(interf.masterZoneIndex()) != inKey(interf.slaveZoneIndex()))
        {
            if (messager::master())
                std::cout << "  [subset-graph] conformal interface '"
                          << interf.name()
                          << "' crosses the zone subset boundary: using full "
                             "graph"
                          << std::endl;
            return true;
        }
    }
    return false;
}

mesh::customGraphEntry&
mesh::findOrCreateCustomGraphEntry_(const std::vector<const zone*>& zones)
{
    const std::vector<label> key = zoneKey_(zones);
    std::string ks;
    for (const label z : key)
        ks += (ks.empty() ? "" : ",") + std::to_string(z);
    for (auto& e : customGraphCache_)
    {
        if (e.zoneKey == key)
        {
            if (messager::master())
                std::cout << "  [subset-graph] reuse shared graph for zones {"
                          << ks << "}" << std::endl;
            return e;
        }
    }
    if (messager::master())
        std::cout << "  [subset-graph] new graph for zones {" << ks << "}"
                  << std::endl;
    customGraphCache_.push_back(customGraphEntry{});
    customGraphCache_.back().zoneKey = key;
    return customGraphCache_.back();
}

void mesh::updateNodeGraph_()
{
    // Need to re-number the nodes only in parallel
    // because there could be new ghosted nodes available
    updateLocalNodeIDs_();

    // Rebuild node connectivities graph
    if (globalOrderGraphPtr_)
    {
        globalOrderGraphPtr_->rebuildGraph();
    }
    if (localOrderGraphPtr_)
    {
        localOrderGraphPtr_->rebuildGraph();
    }

    // rebuild any cached subset graphs (parallel too: re-numbers ghosts)
    for (auto& e : customGraphCache_)
    {
        if (e.local)
            e.local->rebuildGraph();
        if (e.global)
            e.global->rebuildGraph();
    }

    // stencils changed: cached conformal row-merge position maps are stale
    for (label i = 0; i < nInterfaces(); ++i)
    {
        interfaceRef(i).clearConformalRowToRowMaps();
    }
}

void mesh::buildOwnedNodeAdjacency_(std::vector<label>& rowPtr,
                                    std::vector<label>& colIdx) const
{
    const auto& bulkData = this->bulkDataRef();
    const label nOwned = nNodes_;
    const stk::mesh::Selector interiorSel =
        this->universalInteriorPartsSelector();

    // couplings that do not follow from element connectivity: matched pairs of
    // a conformal interface, and the current/opposing elements of a
    // non-conformal one. Mirrors what nodeGraph::buildGraph_ adds to the rows
    std::vector<std::vector<label>> interfaceNbrs(nOwned);
    const auto addInterfaceEdge =
        [&](const stk::mesh::Entity a, const stk::mesh::Entity b)
    {
        const ulabel la = bulkData.local_id(a);
        const ulabel lb = bulkData.local_id(b);
        if (la >= static_cast<ulabel>(nOwned) ||
            lb >= static_cast<ulabel>(nOwned) || la == lb)
            return;
        interfaceNbrs[la].push_back(static_cast<label>(lb));
        interfaceNbrs[lb].push_back(static_cast<label>(la));
    };

    if (hasInterfaces_)
    {
        for (label iInterface = 0; iInterface < nInterfaces(); ++iInterface)
        {
            const interface& interf = interfaceRef(iInterface);

            if (interf.isConformalTreatment())
            {
                for (const auto& nodePair : interf.matchingNodePairVector())
                {
                    addInterfaceEdge(nodePair.first, nodePair.second);
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

                        const stk::mesh::Entity* cn =
                            bulkData.begin_nodes(ip->currentElement_);
                        const label nCur =
                            bulkData.num_nodes(ip->currentElement_);
                        const stk::mesh::Entity* on =
                            bulkData.begin_nodes(ip->opposingElement_);
                        const label nOpp =
                            bulkData.num_nodes(ip->opposingElement_);

                        for (label a = 0; a < nCur; ++a)
                        {
                            for (label b = 0; b < nOpp; ++b)
                            {
                                addInterfaceEdge(cn[a], on[b]);
                            }
                        }
                    }
                }
            };
            addSide(interf.masterInfoRef());
            addSide(interf.slaveInfoRef());
        }
    }

    for (auto& nbrs : interfaceNbrs)
    {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // two passes (count, then fill) so no oversized temporary is needed;
    // mark[] carries the owner row to de-duplicate without a per-row set
    rowPtr.assign(nOwned + 1, 0);
    std::vector<label> mark(nOwned, -1);

    const auto sweep = [&](const bool counting)
    {
        for (label u = 0; u < nOwned; ++u)
        {
            label k = counting ? 0 : rowPtr[u];

            const stk::mesh::Entity node = localNodeIDToEntity_[u];
            const stk::mesh::Entity* elems = bulkData.begin_elements(node);
            const unsigned nElems = bulkData.num_elements(node);

            const auto visit = [&](const label v)
            {
                if (v == u || mark[v] == u)
                    return;
                mark[v] = u;
                if (counting)
                    ++k;
                else
                    colIdx[k++] = v;
            };

            for (unsigned e = 0; e < nElems; ++e)
            {
                const stk::mesh::Entity elem = elems[e];
                if (!interiorSel(bulkData.bucket(elem)))
                    continue;

                const stk::mesh::Entity* en = bulkData.begin_nodes(elem);
                const unsigned nn = bulkData.num_nodes(elem);
                for (unsigned i = 0; i < nn; ++i)
                {
                    const ulabel v = bulkData.local_id(en[i]);
                    if (v < static_cast<ulabel>(nOwned))
                        visit(static_cast<label>(v));
                }
            }

            for (const label v : interfaceNbrs[u])
            {
                visit(v);
            }

            if (counting)
                rowPtr[u + 1] = k;
        }
    };

    sweep(true);
    for (label u = 0; u < nOwned; ++u)
    {
        rowPtr[u + 1] += rowPtr[u];
    }
    colIdx.resize(rowPtr[nOwned]);
    std::fill(mark.begin(), mark.end(), -1);
    sweep(false);
}

std::vector<label> mesh::computeOwnedNodePermutation_() const
{
    const label n = nNodes_;

    // perm[bucket-order id] = new local id; identity until RCM says otherwise
    std::vector<label> perm(n);
    for (label i = 0; i < n; ++i)
    {
        perm[i] = i;
    }

    std::vector<label> rowPtr, colIdx;
    buildOwnedNodeAdjacency_(rowPtr, colIdx);

    // zone of every owned node. Zones are node-disjoint (they only meet through
    // interfaces), so the min() is just a defensive tie-break. Nodes outside
    // every zone land in the trailing group
    const label nGroups = this->nZones() + 1;
    std::vector<label> zoneOf(n, this->nZones());
    {
        const auto& bulkData = this->bulkDataRef();
        for (label iZone = 0; iZone < this->nZones(); ++iZone)
        {
            const stk::mesh::Selector zoneSel =
                this->metaDataRef().locally_owned_part() &
                stk::mesh::selectUnion(this->zonePtr(iZone)->interiorParts());

            for (const auto* bucket :
                 bulkData.get_buckets(stk::topology::NODE_RANK, zoneSel))
            {
                for (label i = 0; i < static_cast<label>(bucket->size()); ++i)
                {
                    const ulabel lid = bulkData.local_id((*bucket)[i]);
                    if (lid < static_cast<ulabel>(n) &&
                        iZone < zoneOf[static_cast<label>(lid)])
                    {
                        zoneOf[static_cast<label>(lid)] = iZone;
                    }
                }
            }
        }
    }

    std::vector<label> deg(n);
    for (label i = 0; i < n; ++i)
    {
        deg[i] = rowPtr[i + 1] - rowPtr[i];
    }

    std::vector<char> numbered(n, 0);
    std::vector<label> levelOf(n, -1);
    std::vector<label> bfsQueue;
    bfsQueue.reserve(n);

    // level-set BFS over the not-yet-numbered nodes of zone g; returns the
    // eccentricity of src and collects the nodes of its last level
    const auto bfs =
        [&](const label src, const label g, std::vector<label>& lastLevel)
    {
        bfsQueue.clear();
        bfsQueue.push_back(src);
        levelOf[src] = 0;
        label maxLevel = 0;

        for (size_t head = 0; head < bfsQueue.size(); ++head)
        {
            const label u = bfsQueue[head];
            for (label k = rowPtr[u]; k < rowPtr[u + 1]; ++k)
            {
                const label v = colIdx[k];
                if (numbered[v] || zoneOf[v] != g || levelOf[v] >= 0)
                    continue;
                levelOf[v] = levelOf[u] + 1;
                maxLevel = std::max(maxLevel, levelOf[v]);
                bfsQueue.push_back(v);
            }
        }

        lastLevel.clear();
        for (const label u : bfsQueue)
        {
            if (levelOf[u] == maxLevel)
                lastLevel.push_back(u);
        }
        for (const label u : bfsQueue)
        {
            levelOf[u] = -1; // reset the scratch for the next sweep
        }
        return maxLevel;
    };

    // George-Liu: walk to a node of (near) maximum eccentricity to start from
    const auto pseudoPeripheralNode = [&](label root, const label g)
    {
        std::vector<label> lastLevel;
        label ecc = bfs(root, g, lastLevel);

        for (label iter = 0; iter < 5 && !lastLevel.empty(); ++iter)
        {
            label candidate = lastLevel[0];
            for (const label u : lastLevel)
            {
                if (deg[u] < deg[candidate])
                    candidate = u;
            }

            std::vector<label> candidateLastLevel;
            const label candidateEcc = bfs(candidate, g, candidateLastLevel);
            if (candidateEcc <= ecc)
                break;

            ecc = candidateEcc;
            root = candidate;
            lastLevel.swap(candidateLastLevel);
        }
        return root;
    };

    std::vector<label> order; // new local id -> bucket-order id
    order.reserve(n);

    for (label g = 0; g < nGroups; ++g)
    {
        const size_t zoneStart = order.size();

        for (label seed = 0; seed < n; ++seed)
        {
            if (numbered[seed] || zoneOf[seed] != g)
                continue;

            // Cuthill-McKee over this component: neighbours of a node enter the
            // queue by ascending degree
            const label root = pseudoPeripheralNode(seed, g);
            numbered[root] = 1;
            order.push_back(root);

            for (size_t head = order.size() - 1; head < order.size(); ++head)
            {
                const label u = order[head];
                const size_t firstNbr = order.size();

                for (label k = rowPtr[u]; k < rowPtr[u + 1]; ++k)
                {
                    const label v = colIdx[k];
                    if (numbered[v] || zoneOf[v] != g)
                        continue;
                    numbered[v] = 1;
                    order.push_back(v);
                }

                std::sort(order.begin() + firstNbr,
                          order.end(),
                          [&deg](const label a, const label b)
                { return deg[a] < deg[b] || (deg[a] == deg[b] && a < b); });
            }
        }

        // reverse this zone's segment only: RCM, and the zone stays contiguous
        std::reverse(order.begin() + zoneStart, order.end());
    }

    assert(static_cast<label>(order.size()) == n);
    for (label newID = 0; newID < n; ++newID)
    {
        perm[order[newID]] = newID;
    }

    // Report what the renumbering bought, over the owned rows. Rows are only
    // re-ordered inside a zone, so the in-zone figures are what RCM controls;
    // couplings across zones (zone interfaces) jump between two contiguous
    // zone ranges and set a floor on the overall bandwidth
    label bw[2] = {0, 0};        // in-zone: before, after
    label bwAll[2] = {0, 0};     // including the cross-zone couplings
    int64_t profile[2] = {0, 0}; // in-zone: before, after
    bool anyCrossZone = false;

    for (label i = 0; i < n; ++i)
    {
        label rowMax[2] = {0, 0};
        for (label k = rowPtr[i]; k < rowPtr[i + 1]; ++k)
        {
            const label j = colIdx[k];
            const label dOld = std::abs(i - j);
            const label dNew = std::abs(perm[i] - perm[j]);

            bwAll[0] = std::max(bwAll[0], dOld);
            bwAll[1] = std::max(bwAll[1], dNew);

            if (zoneOf[i] != zoneOf[j])
            {
                anyCrossZone = true;
                continue;
            }
            rowMax[0] = std::max(rowMax[0], dOld);
            rowMax[1] = std::max(rowMax[1], dNew);
        }
        bw[0] = std::max(bw[0], rowMax[0]);
        bw[1] = std::max(bw[1], rowMax[1]);
        profile[0] += rowMax[0];
        profile[1] += rowMax[1];
    }

    if (messager::parallel())
    {
        int crossZone = anyCrossZone ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE,
                      bw,
                      2,
                      ::linearSolver::MPIDataType<label>::type(),
                      MPI_MAX,
                      messager::comm());
        MPI_Allreduce(MPI_IN_PLACE,
                      bwAll,
                      2,
                      ::linearSolver::MPIDataType<label>::type(),
                      MPI_MAX,
                      messager::comm());
        MPI_Allreduce(
            MPI_IN_PLACE, &crossZone, 1, MPI_INT, MPI_MAX, messager::comm());
        MPI_Allreduce(MPI_IN_PLACE,
                      profile,
                      2,
                      ::linearSolver::MPIDataType<int64_t>::type(),
                      MPI_SUM,
                      messager::comm());
        anyCrossZone = (crossZone == 1);
    }

    if (messager::master())
    {
        std::cout << "\n[node-ordering] reverse Cuthill-McKee (per zone):\n"
                  << "  bandwidth " << bw[0] << " -> " << bw[1] << ", profile "
                  << profile[0] << " -> " << profile[1] << std::endl;

        if (anyCrossZone)
        {
            std::cout << "[node-ordering] zones are numbered contiguously;\n"
                      << "  interface couplings across zones hold the overall "
                      << "bandwidth at " << bwAll[1] << " (was " << bwAll[0]
                      << ")" << std::endl;
        }
    }

    return perm;
}

void mesh::applyOwnedNodePermutation_(const std::vector<label>& perm)
{
    auto& bulkData = this->bulkDataRef();
    const label nOwned = nNodes_;

    assert(static_cast<label>(perm.size()) == nOwned);
    assert(static_cast<label>(localNodeIDToEntity_.size()) >= nOwned);

    // snapshot first: the loop below invalidates the id it reads from
    std::vector<stk::mesh::Entity> ownedNodes(
        localNodeIDToEntity_.begin(), localNodeIDToEntity_.begin() + nOwned);

    for (label oldID = 0; oldID < nOwned; ++oldID)
    {
        const stk::mesh::Entity node = ownedNodes[oldID];
        const label newID = perm[oldID];
        bulkData.set_local_id(node, newID);
        localNodeIDToEntity_[newID] = node;
    }

    // the global id of an owned node must stay offset + local id: the linear
    // system takes the row's own global index from it (see
    // CRSNodeGraph::localToGlobalRow)
    if (!messager::parallel())
    {
        for (label lid = 0; lid < nOwned; ++lid)
        {
            bulkData.set_global_id(localNodeIDToEntity_[lid], lid);
        }
        return;
    }

    label offset = 0;
    label nOwnedNodes = nOwned;
    MPI_Exscan(&nOwnedNodes,
               &offset,
               1,
               ::linearSolver::MPIDataType<label>::type(),
               MPI_SUM,
               messager::comm());

    for (label lid = 0; lid < nOwned; ++lid)
    {
        label* globalIdentity = stk::mesh::field_data(
            *globalIdentityFieldPtr_, localNodeIDToEntity_[lid]);
        *globalIdentity = offset + lid;
    }

    // push the new owner ids out to the shared and ghosted copies
    stk::mesh::communicate_field_data(bulkData, {globalIdentityFieldPtr_});

    for (const auto* bucket : bulkData.get_buckets(
             stk::topology::NODE_RANK, this->universalInteriorPartsSelector()))
    {
        const label* globalIdentityb =
            stk::mesh::field_data(*globalIdentityFieldPtr_, *bucket);
        for (label i = 0; i < static_cast<label>(bucket->size()); ++i)
        {
            bulkData.set_global_id((*bucket)[i], globalIdentityb[i]);
        }
    }
}

void mesh::rebuildLocalNodeIDToEntity_()
{
    const auto& bulkData = this->bulkDataRef();

    localNodeIDToEntity_.assign(nAllNodes_, stk::mesh::Entity::InvalidEntity);

    for (const auto* bucket : bulkData.get_buckets(
             stk::topology::NODE_RANK, this->universalInteriorPartsSelector()))
    {
        for (label i = 0; i < static_cast<label>(bucket->size()); ++i)
        {
            const stk::mesh::Entity node = (*bucket)[i];
            localNodeIDToEntity_[bulkData.local_id(node)] = node;
        }
    }
}

void mesh::initializeLocalNodeIDs_()
{
    // Set new local id for nodes in case of parallel run
    if (!messager::parallel())
    {
        // Important: setting new id's for the nodes in a serial number
        // is only a practice of convenience to make the code in parallel
        // and serial re-usable, however, in thoery no need to do that:
        // in fact we set the global id here to be equal to the local one

        auto& bulkData = this->bulkDataRef();

        // Count number of nodes: for parallel this is number of owned nodes
        // only for the proc
        nNodes_ = stk::mesh::count_entities(
            bulkData,
            stk::topology::NODE_RANK,
            this->locallyOwnedInteriorPartsSelector());

        // set size of local id to node entity map
        localNodeIDToEntity_.resize(nNodes_, stk::mesh::Entity::InvalidEntity);

        nShadowNodes_ = 0;
        nUselessNodes_ = 0;
        nActiveNodes_ = nNodes_;
        nAllNodes_ = nNodes_;

        // Update the global id in stk: need to initialize first
        bulkData.initialize_global_ids();

        const auto& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK,
                                 this->locallyOwnedInteriorPartsSelector());

        // reassign contiguous ids over the active nodes; ignored blocks would
        // otherwise leave gaps in the stk default local ids
        label nodeID_new = 0;
        for (label iBucket = 0;
             iBucket < static_cast<label>(nodeBuckets.size());
             ++iBucket)
        {
            const stk::mesh::Bucket& nodeBucket = *nodeBuckets[iBucket];

            for (label iNode = 0; iNode < static_cast<label>(nodeBucket.size());
                 ++iNode)
            {
                const auto& node = nodeBucket[iNode];

                bulkData.set_local_id(node, nodeID_new);
                bulkData.set_global_id(node, nodeID_new);
                localNodeIDToEntity_[nodeID_new] = node;
                nodeID_new++;
            }
        }

#ifndef NDEBUG
        for (const auto& entity : localNodeIDToEntity_)
        {
            if (entity == stk::mesh::Entity::InvalidEntity)
            {
                errorMsg("localNodeIDToEntity_ not properly populated");
            }
        }
#endif /* NDEBUG */
    }
    else
    {
        // get some MPI info
        label rank = messager::myProcNo();

        auto& bulkData = this->bulkDataRef();
        auto& metaData = this->metaDataRef();

        // Important: Multiple things to be done here:
        // 1- Identifiers of nodes are not contiguous by default
        // meaning that the nodes that belong to the same proc might
        // not have consecutive numbering. Therefore, the global
        // identity container provided by bulkData must be populated
        // and used whenever a global ID is required. The population of the
        // global id array must be contiguous.
        // 2- local identities should also be re-set such that owned
        // nodes are first, then shared not owned, then aura. These
        // arrangements are essential for the matrix assembly
        // 3- There will be aura nodes that are not to be part of the
        // node-graph. These are nodes that are not neighbours to
        // any owned nodes. It happens usually in a ghosted element
        // that is ghosted from a lower rank proc

        // Define new global identities
        {
            // Define the global identity field, fill and synchronize. Use the
            // auxiliary stk field for this purpose
            globalIdentityFieldPtr_ =
                &metaData.declare_field<label>(stk::topology::NODE_RANK, aux);

            // put field on all interior parts

            // Put the stk field on interior mesh parts
            for (const stk::mesh::Part* part : this->interiorActiveParts())
            {
                if (!globalIdentityFieldPtr_->defined_on(*part))
                {
                    stk::mesh::put_field_on_mesh(
                        *globalIdentityFieldPtr_, *part, 1, nullptr);
                }
            }

            // perform the calulation
            {
                // Calculate offset
                label nOwnedNodes = stk::mesh::count_entities(
                    bulkDataRef(),
                    stk::topology::NODE_RANK,
                    this->locallyOwnedInteriorPartsSelector());
                label offset = 0;
                MPI_Exscan(&nOwnedNodes,
                           &offset,
                           1,
                           ::linearSolver::MPIDataType<label>::type(),
                           MPI_SUM,
                           messager::comm());

                // Define new id, counter
                int64_t nodeID_new = offset;

                stk::mesh::BucketVector const& nodeBuckets =
                    bulkData.get_buckets(
                        stk::topology::NODE_RANK,
                        this->locallyOwnedInteriorPartsSelector());
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;

                    const label nNodesPerBucket =
                        static_cast<label>(nodeBucket.size());

                    // field chunks in bucket
                    label* globalIdentityb = stk::mesh::field_data(
                        *globalIdentityFieldPtr_, nodeBucket);

                    for (label iNode = 0; iNode < nNodesPerBucket; ++iNode)
                    {
                        globalIdentityb[iNode] = nodeID_new++;
                    }
                }

                // Synchronize
                stk::mesh::communicate_field_data(bulkData,
                                                  {globalIdentityFieldPtr_});
            }

            // Update the global id in stk: need to initialize first
            bulkData.initialize_global_ids();

            const auto& nodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK,
                                     this->universalInteriorPartsSelector());

            label nBuckets = static_cast<label>(nodeBuckets.size());
            for (label iBucket = 0; iBucket < nBuckets; ++iBucket)
            {
                const stk::mesh::Bucket& nodeBucket = *nodeBuckets[iBucket];

                const label nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                label* globalIdentityb =
                    stk::mesh::field_data(*globalIdentityFieldPtr_, nodeBucket);

                for (label iNode = 0; iNode < nNodesPerBucket; ++iNode)
                {
                    const auto& node = nodeBucket[iNode];

                    bulkData.set_global_id(node, globalIdentityb[iNode]);
                }
            }
        }

        // Set size of local id to entity map. Get first the number of all nodes
        // in this process
        {
            label nAllNodes = stk::mesh::count_entities(
                bulkData,
                stk::topology::NODE_RANK,
                this->universalInteriorPartsSelector());

            localNodeIDToEntity_.resize(nAllNodes,
                                        stk::mesh::Entity::InvalidEntity);
        }

        // Define new local identities
        {
            // Define new id, counter
            label nodeID_new = 0;

            // Sizes
            label nOwnedNodes = 0;
            label nSharedAndNotOwnedNodes = 0;
            label nGhostedNodes = 0;

            // Owned nodes (unchanged)
            const auto& ownedNodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK,
                metaData.locally_owned_part() &
                    stk::mesh::selectUnion(this->interiorActiveParts()));

            label nBuckets = static_cast<label>(ownedNodeBuckets.size());
            for (label iBucket = 0; iBucket < nBuckets; ++iBucket)
            {
                const stk::mesh::Bucket& theBucketRef =
                    *ownedNodeBuckets[iBucket];
                for (label iNode = 0;
                     iNode < static_cast<label>(theBucketRef.size());
                     ++iNode)
                {
                    const auto& node = theBucketRef[iNode];
                    bulkData.set_local_id(node, nodeID_new);
                    localNodeIDToEntity_[nodeID_new] = node;
                    nodeID_new++;
                    nOwnedNodes++;
                }
            }

            // Shared not-owned nodes
            const auto& sharedNotOwnedNodeBucketsPtrArray =
                bulkData.get_buckets(
                    stk::topology::NODE_RANK,
                    (!metaData.locally_owned_part() &
                     metaData.globally_shared_part()) &
                        stk::mesh::selectUnion(this->interiorActiveParts()));

            for (label iBucket = 0;
                 iBucket <
                 static_cast<label>(sharedNotOwnedNodeBucketsPtrArray.size());
                 ++iBucket)
            {
                const stk::mesh::Bucket& theBucketRef =
                    *sharedNotOwnedNodeBucketsPtrArray[iBucket];
                for (label iNode = 0;
                     iNode < static_cast<label>(theBucketRef.size());
                     ++iNode)
                {
                    const auto& node = theBucketRef[iNode];
                    bulkData.set_local_id(node, nodeID_new);
                    localNodeIDToEntity_[nodeID_new] = node;
                    nodeID_new++;
                    nSharedAndNotOwnedNodes++;
                }
            }

            // Ghosted nodes (i.e. aura + custom ghosted)
            const auto& ghostedNodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK,
                (!metaData.locally_owned_part() &
                 !metaData.globally_shared_part()) &
                    stk::mesh::selectUnion(this->interiorActiveParts()));

            for (label iBucket = 0;
                 iBucket < static_cast<label>(ghostedNodeBuckets.size());
                 ++iBucket)
            {
                const stk::mesh::Bucket& theBucketRef =
                    *ghostedNodeBuckets[iBucket];
                for (label iNode = 0;
                     iNode < static_cast<label>(theBucketRef.size());
                     ++iNode)
                {
                    const auto& node = theBucketRef[iNode];
                    bulkData.set_local_id(node, nodeID_new);
                    localNodeIDToEntity_[nodeID_new] = node;
                    nodeID_new++;
                    nGhostedNodes++;
                }
            }

            // Store number of nodes
            nNodes_ = nOwnedNodes;
            nAllNodes_ = nOwnedNodes + nSharedAndNotOwnedNodes + nGhostedNodes;

            assert(nNodes_ == label(stk::mesh::count_entities(
                                  bulkData,
                                  stk::topology::NODE_RANK,
                                  metaData.locally_owned_part() &
                                      stk::mesh::selectUnion(
                                          this->interiorActiveParts()))));
            assert(nAllNodes_ == label(stk::mesh::count_entities(
                                     bulkData,
                                     stk::topology::NODE_RANK,
                                     metaData.universal_part() &
                                         stk::mesh::selectUnion(
                                             this->interiorActiveParts()))));

            // Assign nodes to 1 if they belong to useful elements. Useful
            // elements are those which contain at least one owned node. The
            // objective of this task is to extract useless aura elements as
            // well as any other element which cannot be detected by the owner
            // processor built-in functions. We do it manual.
            std::vector<label> activeNodeFlag(nAllNodes_, 0);

            // set the owned + shared not owned to 1: we are sure about these
            for (label i = 0; i < nOwnedNodes + nSharedAndNotOwnedNodes; i++)
            {
                activeNodeFlag[i] = 1;
            }

            const auto& allElementBuckets = bulkData.get_buckets(
                stk::topology::ELEMENT_RANK,
                !metaData.locally_owned_part() &
                    !metaData.globally_shared_part() &
                    stk::mesh::selectUnion(this->interiorActiveParts()));

            for (size_t iElementBucket = 0;
                 iElementBucket < allElementBuckets.size();
                 ++iElementBucket)
            {
                const stk::mesh::Bucket& theBucketRef =
                    *allElementBuckets[iElementBucket];

                for (label iElement = 0;
                     iElement < static_cast<label>(theBucketRef.size());
                     ++iElement)
                {
                    const auto& element = theBucketRef[iElement];

                    // this aura element will not have any owned node .. we skip
                    // (Trilinos STK convention: the owner of a shared element
                    // is the rank with the _lowest_ id among the ranks the
                    // element is shared with.)
                    if (bulkData.parallel_owner_rank(element) < rank)
                        continue;

                    stk::mesh::Entity const* nodeRels =
                        bulkData.begin_nodes(element);
                    label numNodes = bulkData.num_nodes(element);

                    bool notFound = true;
                    for (label ni = 0; ni < numNodes; ++ni)
                    {
                        stk::mesh::Entity node = nodeRels[ni];
                        if (bulkData.parallel_owner_rank(node) == rank)
                        {
                            notFound = false;
                            break;
                        }
                    }

                    if (!notFound)
                    {
                        // useless element
                        for (label ni = 0; ni < numNodes; ++ni)
                        {
                            stk::mesh::Entity node = nodeRels[ni];
                            stk::mesh::EntityId nodeID =
                                bulkData.local_id(node);
                            if (activeNodeFlag[nodeID] == 0)
                            {
                                activeNodeFlag[nodeID] = 1;
                            }
                        }
                    }
                }
            }

            // Loop over all elements that are ghosted as aura at non-conformal
            // interfaces. The nodes of these elements are ALWAYS useful
            if (hasInterfaces_)
            {
                for (label iInterface = 0; iInterface < nInterfaces();
                     iInterface++)
                {
                    if (interfaceRef(iInterface).interfaceGhosting_)
                    {
                        std::vector<stk::mesh::EntityKey> recvList;
                        interfaceRef(iInterface)
                            .interfaceGhosting_->receive_list(recvList);

                        for (const stk::mesh::EntityKey& key : recvList)
                        {
                            stk::mesh::Entity entity = bulkData.get_entity(key);
                            stk::mesh::EntityRank entityRank =
                                bulkData.entity_rank(entity);
                            if (entityRank ==
                                stk::mesh::EntityRank::ELEMENT_RANK)
                            {
                                stk::mesh::Entity const* nodeRels =
                                    bulkData.begin_nodes(entity);
                                label numNodes = bulkData.num_nodes(entity);
                                for (label ni = 0; ni < numNodes; ++ni)
                                {
                                    stk::mesh::Entity node = nodeRels[ni];
                                    stk::mesh::EntityId nodeID =
                                        bulkData.local_id(node);
                                    activeNodeFlag[nodeID] = 1;
                                }
                            }
                        }
                    }
                }
            }

            // Set some sizes
            nActiveNodes_ = 0;
            for (label i = 0; i < static_cast<label>(activeNodeFlag.size());
                 i++)
            {
                if (activeNodeFlag[i] == 1)
                {
                    nActiveNodes_++;
                }
            }

            nShadowNodes_ = nActiveNodes_ - nNodes_;
            nUselessNodes_ = nAllNodes_ - nActiveNodes_;

            // Re-map all ghosted nodes (shadow nodes) + shift useless to the
            // end (useless will be excluded throughout the code)
            label k = nNodes_;
            label l = nActiveNodes_;

            const auto& notOwnedNodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK,
                ((metaData.universal_part() & !metaData.locally_owned_part()) &
                 stk::mesh::selectUnion(this->interiorActiveParts())));

            for (size_t iNodeBucket = 0;
                 iNodeBucket < notOwnedNodeBuckets.size();
                 ++iNodeBucket)
            {
                const stk::mesh::Bucket& theBucketRef =
                    *notOwnedNodeBuckets[iNodeBucket];

                for (label iNode = 0;
                     iNode < static_cast<label>(theBucketRef.size());
                     ++iNode)
                {
                    const auto& node = theBucketRef[iNode];
                    if (activeNodeFlag[bulkData.local_id(node)] == 1)
                    {
                        localNodeIDToEntity_[k] = node;
                        bulkData.set_local_id(node, k++);
                    }
                    else
                    {
                        localNodeIDToEntity_[l] = node;
                        bulkData.set_local_id(node, l++);
                    }
                }
            }

            assert(k == nActiveNodes_);
            assert(l == nAllNodes_);
        }

#ifndef NDEBUG
        for (const auto& entity : localNodeIDToEntity_)
        {
            if (entity == stk::mesh::Entity::InvalidEntity)
            {
                errorMsg("localNodeIDToEntity_ not properly populated");
            }
        }
#endif /* NDEBUG */
    }

    // Up to here the local id of a node is its position in the STK bucket
    // order, i.e. whatever order the mesh file happened to use. Re-order the
    // owned nodes (the matrix rows) with reverse Cuthill-McKee to shrink the
    // bandwidth. Shared, ghosted and useless nodes keep the slots assigned
    // above, so the | owned | shadow | useless | layout is untouched
    if (this->controlsRef().isRenumbered())
    {
        ownedNodePermutation_ = computeOwnedNodePermutation_();
        applyOwnedNodePermutation_(ownedNodePermutation_);
    }
}

void mesh::updateLocalNodeIDs_()
{
    if (!messager::parallel())
        return;

    // get some MPI info
    label rank = messager::myProcNo();

    auto& bulkData = this->bulkDataRef();
    auto& metaData = this->metaDataRef();

    // 1.) update global IDs
    stk::mesh::communicate_field_data(bulkData, {globalIdentityFieldPtr_});

    bulkData.initialize_global_ids(); // clear and resize
    for (const auto* bucket : bulkData.get_buckets(
             stk::topology::NODE_RANK, this->universalInteriorPartsSelector()))
    {
        const label nNodesPerBucket = static_cast<label>(bucket->size());
        label* globalIdentityb =
            stk::mesh::field_data(*globalIdentityFieldPtr_, *bucket);

        for (label iNode = 0; iNode < nNodesPerBucket; ++iNode)
        {
            const auto& node = (*bucket)[iNode];
            bulkData.set_global_id(node, globalIdentityb[iNode]);
        }
    }

    // 2.) compute new local IDs
    // Define new id, counter
    label nodeID_new = 0;

    // Sizes
    label nOwnedNodes = 0;
    label nSharedAndNotOwnedNodes = 0;
    label nGhostedNodes = 0;

    // Owned nodes (unchanged)
    const auto& ownedNodeBuckets = bulkData.get_buckets(
        stk::topology::NODE_RANK,
        metaData.locally_owned_part() &
            stk::mesh::selectUnion(this->interiorActiveParts()));

    label nBuckets = static_cast<label>(ownedNodeBuckets.size());
    for (label iBucket = 0; iBucket < nBuckets; ++iBucket)
    {
        const stk::mesh::Bucket& theBucketRef = *ownedNodeBuckets[iBucket];
        for (label iNode = 0; iNode < static_cast<label>(theBucketRef.size());
             ++iNode)
        {
            const auto& node = theBucketRef[iNode];
            bulkData.set_local_id(node, nodeID_new);
            nodeID_new++;
            nOwnedNodes++;
        }
    }

    // Shared not-owned nodes
    const auto& sharedNotOwnedNodeBucketsPtrArray = bulkData.get_buckets(
        stk::topology::NODE_RANK,
        (!metaData.locally_owned_part() & metaData.globally_shared_part()) &
            stk::mesh::selectUnion(this->interiorActiveParts()));

    for (label iBucket = 0;
         iBucket < static_cast<label>(sharedNotOwnedNodeBucketsPtrArray.size());
         ++iBucket)
    {
        const stk::mesh::Bucket& theBucketRef =
            *sharedNotOwnedNodeBucketsPtrArray[iBucket];
        for (label iNode = 0; iNode < static_cast<label>(theBucketRef.size());
             ++iNode)
        {
            const auto& node = theBucketRef[iNode];
            bulkData.set_local_id(node, nodeID_new);
            nodeID_new++;
            nSharedAndNotOwnedNodes++;
        }
    }

    // Ghosted nodes (i.e. aura + custom ghosted)
    const auto& ghostedNodeBuckets = bulkData.get_buckets(
        stk::topology::NODE_RANK,
        (!metaData.locally_owned_part() & !metaData.globally_shared_part()) &
            stk::mesh::selectUnion(this->interiorActiveParts()));

    for (label iBucket = 0;
         iBucket < static_cast<label>(ghostedNodeBuckets.size());
         ++iBucket)
    {
        const stk::mesh::Bucket& theBucketRef = *ghostedNodeBuckets[iBucket];
        for (label iNode = 0; iNode < static_cast<label>(theBucketRef.size());
             ++iNode)
        {
            const auto& node = theBucketRef[iNode];
            bulkData.set_local_id(node, nodeID_new);
            nodeID_new++;
            nGhostedNodes++;
        }
    }

    // Store number of nodes
    nNodes_ = nOwnedNodes;
    nAllNodes_ = nOwnedNodes + nSharedAndNotOwnedNodes + nGhostedNodes;

    assert(nNodes_ ==
           label(stk::mesh::count_entities(
               bulkData,
               stk::topology::NODE_RANK,
               metaData.locally_owned_part() &
                   stk::mesh::selectUnion(this->interiorActiveParts()))));
    assert(nAllNodes_ ==
           label(stk::mesh::count_entities(
               bulkData,
               stk::topology::NODE_RANK,
               metaData.universal_part() &
                   stk::mesh::selectUnion(this->interiorActiveParts()))));

    // Assign nodes to 1 if they belong to useful elements. Useful
    // elements are those which contain at least one owned node. The
    // objective of this task is to extract useless aura elements as
    // well as any other element which cannot be detected by the owner
    // processor built-in functions. We do it manual.
    std::vector<label> activeNodeFlag(nAllNodes_, 0);

    // set the owned + shared not owned to 1: we are sure about these
    for (label i = 0; i < nOwnedNodes + nSharedAndNotOwnedNodes; i++)
    {
        activeNodeFlag[i] = 1;
    }

    const auto& allElementBuckets = bulkData.get_buckets(
        stk::topology::ELEMENT_RANK,
        !metaData.locally_owned_part() & !metaData.globally_shared_part() &
            stk::mesh::selectUnion(this->interiorActiveParts()));

    for (size_t iElementBucket = 0; iElementBucket < allElementBuckets.size();
         ++iElementBucket)
    {
        const stk::mesh::Bucket& theBucketRef =
            *allElementBuckets[iElementBucket];

        for (label iElement = 0;
             iElement < static_cast<label>(theBucketRef.size());
             ++iElement)
        {
            const auto& element = theBucketRef[iElement];

            // this aura element will not have any owned node .. we skip
            // (Trilinos STK convention: the owner of a shared element
            // is the rank with the _lowest_ id among the ranks the
            // element is shared with.)
            if (bulkData.parallel_owner_rank(element) < rank)
                continue;

            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            bool notFound = true;
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];
                if (bulkData.parallel_owner_rank(node) == rank)
                {
                    notFound = false;
                    break;
                }
            }

            if (!notFound)
            {
                // useless element
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];
                    stk::mesh::EntityId nodeID = bulkData.local_id(node);
                    if (activeNodeFlag[nodeID] == 0)
                    {
                        activeNodeFlag[nodeID] = 1;
                    }
                }
            }
        }
    }

    // Loop over all elements that are ghosted as aura at non-conformal
    // interfaces. The nodes of these elements are ALWAYS useful
    if (hasInterfaces_)
    {
        for (label iInterface = 0; iInterface < nInterfaces(); iInterface++)
        {
            if (interfaceRef(iInterface).interfaceGhosting_)
            {
                std::vector<stk::mesh::EntityKey> recvList;
                interfaceRef(iInterface)
                    .interfaceGhosting_->receive_list(recvList);

                for (const stk::mesh::EntityKey& key : recvList)
                {
                    stk::mesh::Entity entity = bulkData.get_entity(key);
                    stk::mesh::EntityRank entityRank =
                        bulkData.entity_rank(entity);
                    if (entityRank == stk::mesh::EntityRank::ELEMENT_RANK)
                    {
                        stk::mesh::Entity const* nodeRels =
                            bulkData.begin_nodes(entity);
                        label numNodes = bulkData.num_nodes(entity);
                        for (label ni = 0; ni < numNodes; ++ni)
                        {
                            stk::mesh::Entity node = nodeRels[ni];
                            stk::mesh::EntityId nodeID =
                                bulkData.local_id(node);
                            activeNodeFlag[nodeID] = 1;
                        }
                    }
                }
            }
        }
    }

    // Set some sizes
    nActiveNodes_ = 0;
    for (label i = 0; i < static_cast<label>(activeNodeFlag.size()); i++)
    {
        if (activeNodeFlag[i] == 1)
        {
            nActiveNodes_++;
        }
    }

    nShadowNodes_ = nActiveNodes_ - nNodes_;
    nUselessNodes_ = nAllNodes_ - nActiveNodes_;

    // Re-map all ghosted nodes (shadow nodes) + shift useless to the
    // end (useless will be excluded throughout the code)
    label k = nNodes_;
    label l = nActiveNodes_;

    const auto& notOwnedNodeBuckets = bulkData.get_buckets(
        stk::topology::NODE_RANK,
        ((metaData.universal_part() & !metaData.locally_owned_part()) &
         stk::mesh::selectUnion(this->interiorActiveParts())));

    for (size_t iNodeBucket = 0; iNodeBucket < notOwnedNodeBuckets.size();
         ++iNodeBucket)
    {
        const stk::mesh::Bucket& theBucketRef =
            *notOwnedNodeBuckets[iNodeBucket];

        for (label iNode = 0; iNode < static_cast<label>(theBucketRef.size());
             ++iNode)
        {
            const auto& node = theBucketRef[iNode];
            if (activeNodeFlag[bulkData.local_id(node)] == 1)
            {
                bulkData.set_local_id(node, k++);
            }
            else
            {
                bulkData.set_local_id(node, l++);
            }
        }
    }

    assert(k == nActiveNodes_);
    assert(l == nAllNodes_);

    // the local ids just changed, and ghosting may have brought in new nodes
    rebuildLocalNodeIDToEntity_();

    // re-apply the bandwidth-reducing order, which the bucket-order pass above
    // has just undone. Mesh motion does not change the topology, so the
    // permutation is only recomputed if the owned node set itself changed.
    // The decision is reduced: computeOwnedNodePermutation_ is collective
    if (this->controlsRef().isRenumbered())
    {
        int staleOrdering =
            (static_cast<label>(ownedNodePermutation_.size()) != nNodes_) ? 1
                                                                          : 0;
        MPI_Allreduce(MPI_IN_PLACE,
                      &staleOrdering,
                      1,
                      MPI_INT,
                      MPI_MAX,
                      messager::comm());

        if (staleOrdering)
        {
            ownedNodePermutation_ = computeOwnedNodePermutation_();
        }
        applyOwnedNodePermutation_(ownedNodePermutation_);
    }
}

} // namespace accel
