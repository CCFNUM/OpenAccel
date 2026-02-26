// File       : meshNodeGraph.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#ifdef HAS_INTERFACE
#include "interface.h"
#include "interfaceSideInfo.h"
#endif /* HAS_INTERFACE */
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

        label nNodesPerBucket = static_cast<label>(nodeBuckets.size());
        for (label iBucket = 0; iBucket < nNodesPerBucket; ++iBucket)
        {
            const stk::mesh::Bucket& nodeBucket = *nodeBuckets[iBucket];

            const label nNodesPerBucket = static_cast<label>(nodeBucket.size());

            for (label iNode = 0; iNode < nNodesPerBucket; ++iNode)
            {
                const auto& node = nodeBucket[iNode];

                // Set global ID to local ID. FIXME: For deactivated
                // element blocks, the local id must also be redefined.
                bulkData.set_global_id(node, bulkData.local_id(node));

                // set the map value for this node
                localNodeIDToEntity_[bulkData.local_id(node)] = node;
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

#ifdef HAS_INTERFACE
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
#endif /* HAS_INTERFACE */

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

#ifdef HAS_INTERFACE
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
#endif /* HAS_INTERFACE */

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
}

} // namespace accel
