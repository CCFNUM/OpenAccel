// File       : interface.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

#include "interface.h"
#include "messager.h"
#include "surfaceComparator.h"

namespace accel
{

void addDownwardRelations(const stk::mesh::BulkData& bulk,
                          std::vector<stk::mesh::EntityKey>& entityKeys);

void keepElemsNotAlreadyGhosted(const stk::mesh::BulkData& bulk,
                                const stk::mesh::EntityProcVec& alreadyGhosted,
                                stk::mesh::EntityProcVec& elemsToGhost);

void fillSendGhostsToRemoveFromGhosting(
    const stk::mesh::EntityProcVec& curSendGhosts,
    const stk::mesh::EntityProcVec& intersection,
    stk::mesh::EntityProcVec& sendGhostsToRemove);

void communicateToFillRecvGhostsToRemove(
    const stk::mesh::BulkData& bulk,
    const stk::mesh::EntityProcVec& sendGhostsToRemove,
    std::vector<stk::mesh::EntityKey>& recvGhostsToRemove);

void keepOnlyElems(const stk::mesh::BulkData& bulk,
                   stk::mesh::EntityProcVec& entityProcs);

void computePreciseGhostingLists(
    const stk::mesh::BulkData& bulk,
    stk::mesh::EntityProcVec& elemsToGhost,
    stk::mesh::EntityProcVec& curSendGhosts,
    std::vector<stk::mesh::EntityKey>& recvGhostsToRemove);

interface::interface(mesh* meshPtr, label index, std::string name)
    : meshPtr_(meshPtr), name_(name), index_(index), interfaceGhosting_(nullptr)
{
}

// Methods

void interface::setup()
{
    masterInfoPtr_->setup();
    slaveInfoPtr_->setup();
}

void interface::initialize()
{
    // determine master-slave geometric relations
    determineGeometricRelations_();

    if (!isConformalTreatment())
    {
        // initialize interface sides
        masterInfoPtr_->initialize();
        slaveInfoPtr_->initialize();

        // initialize ghostings for the interface
        if (messager::parallel())
        {
            initializeGhostings_();
        }

        // complete search
        masterInfoPtr_->completeSearch();
        slaveInfoPtr_->completeSearch();

        // for conformal interfaces, determine the opposing gauss point
        // ids so that the transfer can directly copy ip-to-ip
        if (isConformal_)
        {
            masterInfoPtr_->determineOpposingGaussPointIds();
            slaveInfoPtr_->determineOpposingGaussPointIds();
        }
    }
}

void interface::update()
{
    // update geometric relations
    determineGeometricRelations_();

    if (!isConformalTreatment())
    {
        // initialize interface sides
        masterInfoPtr_->update();
        slaveInfoPtr_->update();

        // update ghostings for the interface
        if (messager::parallel())
        {
            updateGhostings_();
        }

        // complete search
        masterInfoPtr_->completeSearch();
        slaveInfoPtr_->completeSearch();

        // for conformal interfaces, determine the opposing gauss point
        // ids so that the transfer can directly copy ip-to-ip
        if (isConformal_)
        {
            masterInfoPtr_->determineOpposingGaussPointIds();
            slaveInfoPtr_->determineOpposingGaussPointIds();
        }
    }
}

// Operations

void interface::completeSearch()
{
    this->masterInfoRef().completeSearch();
    this->slaveInfoRef().completeSearch();
}

void interface::provideDiagnosis()
{
    this->masterInfoRef().provideDiagnosis();
    this->slaveInfoRef().provideDiagnosis();
}

void interface::errorCheck()
{
    this->masterInfoRef().errorCheck();
    this->slaveInfoRef().errorCheck();
}

void interface::initializeGhostings_()
{
    // ensure no duplications in ghosted elements
    stk::util::sort_and_unique(elemsToGhost_);

    // check for ghosting need
    size_t nGhostedElements = elemsToGhost_.size();
    size_t nGlobalGhostedElements = 0;
    stk::all_reduce_sum(
        MPI_COMM_WORLD, &nGhostedElements, &nGlobalGhostedElements, 1);

    if (nGlobalGhostedElements > 0)
    {
        if (messager::master())
            std::cout << "  DG algorithm will ghost " << nGlobalGhostedElements
                      << " entities: " << std::endl;

        stk::mesh::BulkData& bulkData = meshPtr_->bulkDataRef();

        bulkData.modification_begin();

        // create new ghosting
        std::string theGhostName = "interface_ghosting" + name_;
        interfaceGhosting_ = &bulkData.create_ghosting(theGhostName);

        bulkData.change_ghosting(*interfaceGhosting_, elemsToGhost_);

        bulkData.modification_end();

        ops::populateGhostCommProcs(
            bulkData, *interfaceGhosting_, ghostCommProcs_);

        assert(interfaceGhosting_);

        // ensure that the coordinates for the ghosted elements (required
        // for the fine search) are up-to-date
        const auto* coordsSTKFieldPtr =
            meshPtr_->metaDataRef().get_field<scalar>(
                stk::topology::NODE_RANK, meshPtr_->getCoordinateFieldName());

        std::vector<const stk::mesh::FieldBase*> fieldVec = {coordsSTKFieldPtr};
        stk::mesh::communicate_field_data(*interfaceGhosting_, fieldVec);
    }
}

void interface::updateGhostings_()
{
    // ensure no duplications in ghosted elements
    stk::util::sort_and_unique(elemsToGhost_);

    std::vector<stk::mesh::EntityKey> recvGhostsToRemove;

    if (interfaceGhosting_)
    {
        stk::mesh::EntityProcVec currentSendGhosts;

        interfaceGhosting_->send_list(currentSendGhosts);

        // We want both elemsToGhost_ to only contain elements not already
        // ghosted, and a list of receive-ghosts that no longer need to be
        // ghosted.
        computePreciseGhostingLists(meshPtr_->bulkDataRef(),
                                    elemsToGhost_,
                                    currentSendGhosts,
                                    recvGhostsToRemove);
    }

    // check for ghosting need
    size_t local[2] = {elemsToGhost_.size(), recvGhostsToRemove.size()};
    size_t global[2] = {0, 0};
    stk::all_reduce_sum(MPI_COMM_WORLD, local, global, 2);

    if (global[0] > 0 || global[1] > 0)
    {
        if (messager::master())
            std::cout << "DG algorithm will ghost a new number of entities: "
                      << global[0] << " and remove " << global[1]
                      << " entities from ghosting." << std::endl;

        stk::mesh::BulkData& bulkData = meshPtr_->bulkDataRef();

        bulkData.modification_begin();

        // unlikely to happen, but must be checked for
        if (!interfaceGhosting_)
        {
            // create new ghosting
            std::string theGhostName = "interface_ghosting" + name_;
            interfaceGhosting_ = &bulkData.create_ghosting(theGhostName);
        }

        bulkData.change_ghosting(
            *interfaceGhosting_, elemsToGhost_, recvGhostsToRemove);

        bulkData.modification_end();

        ops::populateGhostCommProcs(
            bulkData, *interfaceGhosting_, ghostCommProcs_);
    }

    // ensure that the coordinates for the ghosted elements (required for
    // the fine search) are up-to-date
    if (interfaceGhosting_)
    {
        const auto* coordsSTKFieldPtr =
            meshPtr_->metaDataRef().get_field<scalar>(
                stk::topology::NODE_RANK, meshPtr_->getCoordinateFieldName());

        std::vector<const stk::mesh::FieldBase*> fieldVec = {coordsSTKFieldPtr};
        stk::mesh::communicate_field_data(*interfaceGhosting_, fieldVec);
    }
}

void interface::determineGeometricRelations_()
{
    // create a surface comparator object to determine translation
    // vector and/or rotation matrix
    auto cmp = utils::surfaceComparator(masterInfoPtr_->currentPartVec_,
                                        slaveInfoPtr_->currentPartVec_);

    utils::vector translationVector;
    utils::matrix rotationMatrix;

    switch (option_)
    {
        case interfaceModelOption::rotationalPeriodicity:
            {
                scalar angle =
                    cmp.determineSeparationAngle(rotationAxis_, axisLocation_);

                if (messager::master())
                {
                    std::cout << std::endl
                              << this->name() << " => Rotational Periodicity"
                              << std::endl;
                    std::cout << "\tRotation angle (slave to master): "
                              << -angle * 180.0 / M_PI << " deg" << std::endl;
                }

                // calculate translation vector and rotation matrix
                translationVector = axisLocation_;
                rotationMatrix =
                    utils::getRotationMatrix(-angle, rotationAxis_);

                utils::vector result =
                    utils::transformVector(rotationMatrix, translationVector);
                translationVector -= result;
            }
            break;

        case interfaceModelOption::translationalPeriodicity:
            {
                // identity only
                rotationMatrix = utils::matrix::Identity();

                // get translation vector = sepration vector from slave to
                // master
                translationVector = cmp.determineSeparationVector();

                if (messager::master())
                {
                    std::cout << std::endl
                              << this->name() << " => Translational Periodicity"
                              << std::endl;
                    std::cout << "\tTranslation vector (slave to master): "
                              << translationVector << std::endl;
                }
            }
            break;

        case interfaceModelOption::generalConnection:
            {
                // identity only
                rotationMatrix = utils::matrix::Identity();

                if (messager::master())
                {
                    std::cout << std::endl
                              << this->name() << " => General Connection"
                              << std::endl;
                }
            }
            break;
    }

    // copy translation vector and rotation matrix to master and slave sides
    this->masterInfoRef().translationVector_ = translationVector;
    this->masterInfoRef().rotationMatrix_ = rotationMatrix;
    this->slaveInfoRef().translationVector_ = -translationVector;
    this->slaveInfoRef().rotationMatrix_ = rotationMatrix.transpose();

#if 0
    // check if master and slave are totally overlapping
    isOverlap_ =
        cmp.checkOverlap(overlapTolerance_, translationVector, rotationMatrix);

    if (messager::master() && !isOverlap_)
    {
        std::cout << "\t* nonoverlap configuration detected" << std::endl;
    }
#endif

    // check if the interface is conformal and fill the matching node pairs
    isConformal_ = cmp.checkConformality(matchingNodePairVector_,
                                         conformalityTolerance_,
                                         translationVector,
                                         rotationMatrix);

    if (messager::master() && isConformal_)
    {
        if (isForceNonconformalTreatment_)
        {
            std::cout
                << "\t* conformal configuration detected, however, a "
                   "non-conformal treatment \n\t* of the interface will be "
                   "processed for this interface"
                << std::endl;
        }
        else
        {
            if (messager::parallel())
            {
                errorMsg(
                    "conformal treatment is not enabled for parallel runs.");
            }
            else
            {
                std::cout << "\t* conformal configuration detected"
                          << std::endl;
            }
        }
    }
}

void interface::populateConformalRowToRowMapping_()
{
#ifndef NDEBUG
    // sanity check
    if (!isConformalTreatment())
    {
        errorMsg("Wrong call for populateConformalRowToRowMapping_ method for "
                 "non-conformal interface");
    }
#endif

    stk::mesh::BulkData& bulkData = meshPtr_->bulkDataRef();
    stk::mesh::MetaData& metaData = meshPtr_->metaDataRef();

    // get local order graph
    const auto* localOrderGraphPtr = this->meshPtr()->getLocalOrderGraphPtr();

    // get local id to node mapper from mesh
    const auto& localNodeIDToEntity = this->meshPtr()->localNodeIDToEntity();

    // set size
    conformalRowToRowMap_.resize(matchingNodePairVector_.size());

    label iPair = 0;
    for (const auto& nodePair : matchingNodePairVector_)
    {
        // get required local data for the matching pair

        // data for stencil 1 (of node 1)
        const auto& node1 = nodePair.first;
        const label& lid1 = bulkData.local_id(node1);
        const auto cols = localOrderGraphPtr->rowLocalIndices(lid1);

        // data for stencil 2 (of node 2)
        const auto& node2 = nodePair.second;
        const label& lid2 = bulkData.local_id(node2);

        // set size of the local mapper
        conformalRowToRowMap_[iPair].resize(cols.size(), -1);

        // populate the mapper
        for (label iCol = 0; iCol < cols.size(); iCol++)
        {
            // get the node that corresponds to the current column of stencil 2
            const auto& node_s2 = localNodeIDToEntity[cols[iCol]];

            // search for the matching node
            auto it = std::find_if(matchingNodePairVector_.begin(),
                                   matchingNodePairVector_.end(),
                                   [&](const auto& p)
            { return p.second == node_s2; });

            if (it != matchingNodePairVector_.end())
            {
                // get matching node in stencil 1
                const auto& node_s1 = it->first;
                label lid_s1 = bulkData.local_id(node_s1);

                // get the corresponding column order in the stencil of node 1
                label iCol_s1;
                for (iCol_s1 = 0; iCol_s1 < cols.size(); iCol_s1++)
                {
                    if (cols[iCol_s1] == lid_s1)
                    {
                        break;
                    }
                }

                // store in the mapper
                conformalRowToRowMap_[iPair][iCol] = iCol_s1;
            }
            else
            {
                // same
                label iCol_s1 = iCol;

                // the current column corresponds to a node
                // that is not sitting on the interface
                conformalRowToRowMap_[iPair][iCol] = iCol_s1;
            }
        }
#ifndef NDEBUG
        // diagnostics
        for (auto& m : conformalRowToRowMap_[iPair])
        {
            if (m < 0)
            {
                errorMsg("mapper is not properly populated");
            }
        }
#endif /* NDEBUG */

        // increment the pair counter
        iPair++;
    }
}

// Access

mesh* interface::meshPtr()
{
    return meshPtr_;
}

const mesh* interface::meshPtr() const
{
    return meshPtr_;
}

mesh& interface::meshRef()
{
    return *meshPtr_;
}

const mesh& interface::meshRef() const
{
    return *meshPtr_;
}

void addDownwardRelations(const stk::mesh::BulkData& bulk,
                          std::vector<stk::mesh::EntityKey>& entityKeys)
{
    size_t numEntities = entityKeys.size();
    for (size_t i = 0; i < numEntities; ++i)
    {
        stk::mesh::Entity ent = bulk.get_entity(entityKeys[i]);
        if (bulk.is_valid(ent))
        {
            stk::mesh::EntityRank thisRank = bulk.entity_rank(ent);

            for (stk::mesh::EntityRank irank = stk::topology::NODE_RANK;
                 irank < thisRank;
                 ++irank)
            {
                unsigned num = bulk.num_connectivity(ent, irank);
                const stk::mesh::Entity* downwardEntities =
                    bulk.begin(ent, irank);

                for (unsigned j = 0; j < num; ++j)
                {
                    stk::mesh::EntityKey key =
                        bulk.entity_key(downwardEntities[j]);
                    const stk::mesh::Bucket& bkt =
                        bulk.bucket(downwardEntities[j]);
                    if (!bkt.shared())
                    {
                        entityKeys.push_back(key);
                    }
                }
            }
        }
    }
}

void keepElemsNotAlreadyGhosted(const stk::mesh::BulkData& bulk,
                                const stk::mesh::EntityProcVec& alreadyGhosted,
                                stk::mesh::EntityProcVec& elemsToGhost)
{
    if (!alreadyGhosted.empty())
    {
        size_t numKept = 0;
        size_t num = elemsToGhost.size();
        for (size_t i = 0; i < num; ++i)
        {
            if (!std::binary_search(alreadyGhosted.begin(),
                                    alreadyGhosted.end(),
                                    elemsToGhost[i]))
            {
                elemsToGhost[numKept++] = elemsToGhost[i];
            }
        }
        elemsToGhost.resize(numKept);
    }
}

void fillSendGhostsToRemoveFromGhosting(
    const stk::mesh::EntityProcVec& curSendGhosts,
    const stk::mesh::EntityProcVec& intersection,
    stk::mesh::EntityProcVec& sendGhostsToRemove)
{
    sendGhostsToRemove.reserve(curSendGhosts.size() - intersection.size());
    for (size_t i = 0; i < curSendGhosts.size(); ++i)
    {
        if (!std::binary_search(
                intersection.begin(), intersection.end(), curSendGhosts[i]))
        {
            sendGhostsToRemove.push_back(curSendGhosts[i]);
        }
    }
}

void communicateToFillRecvGhostsToRemove(
    const stk::mesh::BulkData& bulk,
    const stk::mesh::EntityProcVec& sendGhostsToRemove,
    std::vector<stk::mesh::EntityKey>& recvGhostsToRemove)
{
    stk::CommSparse commSparse(bulk.parallel());
    stk::pack_and_communicate(commSparse,
                              [&]()
    {
        for (const stk::mesh::EntityProc& entityProc : sendGhostsToRemove)
        {
            stk::mesh::EntityKey key = bulk.entity_key(entityProc.first);
            stk::CommBuffer& buf = commSparse.send_buffer(entityProc.second);
            buf.pack<stk::mesh::EntityKey>(key);
        }
    });

    int numProcs = bulk.parallel_size();
    for (int p = 0; p < numProcs; ++p)
    {
        if (p == bulk.parallel_rank())
        {
            continue;
        }
        stk::CommBuffer& buf = commSparse.recv_buffer(p);
        while (buf.remaining())
        {
            stk::mesh::EntityKey key;
            buf.unpack<stk::mesh::EntityKey>(key);
            recvGhostsToRemove.push_back(key);
        }
    }

    addDownwardRelations(bulk, recvGhostsToRemove);
}

void keepOnlyElems(const stk::mesh::BulkData& bulk,
                   stk::mesh::EntityProcVec& entityProcs)
{
    size_t elemCounter = 0;
    for (size_t i = 0; i < entityProcs.size(); ++i)
    {
        if (bulk.entity_rank(entityProcs[i].first) == stk::topology::ELEM_RANK)
        {
            entityProcs[elemCounter++] = entityProcs[i];
        }
    }
    entityProcs.resize(elemCounter);
}

void computePreciseGhostingLists(
    const stk::mesh::BulkData& bulk,
    stk::mesh::EntityProcVec& elemsToGhost,
    stk::mesh::EntityProcVec& curSendGhosts,
    std::vector<stk::mesh::EntityKey>& recvGhostsToRemove)
{
    keepOnlyElems(bulk, curSendGhosts);
    stk::util::sort_and_unique(curSendGhosts);
    stk::util::sort_and_unique(elemsToGhost);

    stk::mesh::EntityProcVec intersection;
    std::set_intersection(curSendGhosts.begin(),
                          curSendGhosts.end(),
                          elemsToGhost.begin(),
                          elemsToGhost.end(),
                          std::back_inserter(intersection));

    keepElemsNotAlreadyGhosted(bulk, intersection, elemsToGhost);

    stk::mesh::EntityProcVec sendGhostsToRemove;
    fillSendGhostsToRemoveFromGhosting(
        curSendGhosts, intersection, sendGhostsToRemove);

    communicateToFillRecvGhostsToRemove(
        bulk, sendGhostsToRemove, recvGhostsToRemove);
}

} // namespace accel

#endif /* HAS_INTERFACE */
