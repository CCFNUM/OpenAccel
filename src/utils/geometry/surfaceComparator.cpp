// File       : surfaceComparator.h
// Created    : Tue Sep 16 2025 14:15:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "surfaceComparator.h"
#include "messager.h"

namespace accel
{

#if SPATIAL_DIM == 3

namespace utils
{

surfaceComparator::surfaceComparator(stk::mesh::PartVector surface1Parts,
                                     stk::mesh::PartVector surface2Parts)
    : bulkData_(surface1Parts.back()->mesh_bulk_data()),
      metaData_(surface1Parts.back()->mesh_meta_data()),
      surface1Parts_(surface1Parts), surface2Parts_(surface2Parts),
      coordsSTKFieldPtr_(
          metaData_.get_field<scalar>(stk::topology::NODE_RANK,
                                      metaData_.coordinate_field_name()))
{
}

bool surfaceComparator::checkOverlap(scalar overlapCheckSearchTolerance,
                                     vector sepVec,
                                     matrix rotMat)
{
    // Get boundary/edge nodes for both surfaces
    auto getBoundaryNodes = [&](const stk::mesh::PartVector& surfaceParts,
                                std::set<stk::mesh::Entity>& boundaryNodes)
    {
        // Map to track edge usage count
        std::map<std::pair<stk::mesh::EntityId, stk::mesh::EntityId>, int>
            edgeCount;

        stk::mesh::BucketVector const& sideBuckets =
            bulkData_.get_buckets(metaData_.side_rank(),
                                  metaData_.locally_owned_part() &
                                      stk::mesh::selectUnion(surfaceParts));

        // Count edge occurrences across all sides
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;
            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                stk::mesh::Entity side = sideBucket[iSide];
                const stk::mesh::Entity* nodes = bulkData_.begin_nodes(side);
                const unsigned numNodes = bulkData_.num_nodes(side);

                for (unsigned i = 0; i < numNodes; ++i)
                {
                    unsigned next = (i + 1) % numNodes;
                    stk::mesh::EntityId id1 = bulkData_.identifier(nodes[i]);
                    stk::mesh::EntityId id2 = bulkData_.identifier(nodes[next]);

                    // Create ordered pair for edge
                    auto edge = (id1 < id2) ? std::make_pair(id1, id2)
                                            : std::make_pair(id2, id1);
                    edgeCount[edge]++;
                }
            }
        }

        // Boundary edges appear only once
        for (const auto& edgePair : edgeCount)
        {
            if (edgePair.second == 1)
            {
                stk::mesh::Entity node1 = bulkData_.get_entity(
                    stk::topology::NODE_RANK, edgePair.first.first);
                stk::mesh::Entity node2 = bulkData_.get_entity(
                    stk::topology::NODE_RANK, edgePair.first.second);

                if (bulkData_.is_valid(node1))
                    boundaryNodes.insert(node1);
                if (bulkData_.is_valid(node2))
                    boundaryNodes.insert(node2);
            }
        }
    };

    // Get boundary nodes for both surfaces
    std::set<stk::mesh::Entity> surface1BoundaryNodes;
    std::set<stk::mesh::Entity> surface2BoundaryNodes;

    getBoundaryNodes(surface1Parts_, surface1BoundaryNodes);
    getBoundaryNodes(surface2Parts_, surface2BoundaryNodes);

    // Create bounding spheres for boundary nodes
    std::vector<sphereBoundingBox> sphereBoundingBoxSurface1Vec;
    std::vector<sphereBoundingBox> sphereBoundingBoxSurface2Vec;

    Point nodeCenter;

    // Surface1 boundary nodes
    for (auto node : surface1BoundaryNodes)
    {
        theEntityKey theIdent(bulkData_.entity_key(node), messager::myProcNo());

        const scalar* coords = stk::mesh::field_data(*coordsSTKFieldPtr_, node);
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            nodeCenter[j] = coords[j];
        }

        sphereBoundingBox theSphere(
            Sphere(nodeCenter, overlapCheckSearchTolerance), theIdent);
        sphereBoundingBoxSurface1Vec.push_back(theSphere);
    }

    // Surface2 boundary nodes: rotate and translate to align with surface1
    utils::vector rotCoords;

    for (auto node : surface2BoundaryNodes)
    {
        theEntityKey theIdent(bulkData_.entity_key(node), messager::myProcNo());

        const utils::vectorViewC coords(
            stk::mesh::field_data(*coordsSTKFieldPtr_, node));
        // rotate coordinates
        rotCoords = utils::transformVector(rotMat, coords);

        // create the bounding point that is rotated and translated
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            const scalar xj = rotCoords(j);
            nodeCenter[j] = xj + sepVec[j];
        }

        sphereBoundingBox theSphere(
            Sphere(nodeCenter, overlapCheckSearchTolerance), theIdent);
        sphereBoundingBoxSurface2Vec.push_back(theSphere);
    }

    // Perform coarse search
    stk::search::SearchMethod searchMethod = stk::search::KDTREE;
    if (conformalityCheckSearchMethodName_ != "stk_kdtree")
    {
        errorMsg("searchMethod only supports stk_kdtree");
    }

    searchKeyVector searchKeyPair;
    stk::search::coarse_search(sphereBoundingBoxSurface1Vec,
                               sphereBoundingBoxSurface2Vec,
                               searchMethod,
                               messager::comm(),
                               searchKeyPair);

    // Count unique matching boundary nodes
    std::set<stk::mesh::EntityKey> matchedSurface1Nodes;
    std::set<stk::mesh::EntityKey> matchedSurface2Nodes;

    for (size_t i = 0, size = searchKeyPair.size(); i < size; ++i)
    {
        matchedSurface1Nodes.insert(searchKeyPair[i].first.id());
        matchedSurface2Nodes.insert(searchKeyPair[i].second.id());
    }

    // Count owned matched boundary nodes
    label nMatchedSurface1BoundaryNodes = 0;
    label nMatchedSurface2BoundaryNodes = 0;

    for (const auto& key : matchedSurface1Nodes)
    {
        stk::mesh::Entity node = bulkData_.get_entity(key);
        if (bulkData_.is_valid(node) &&
            bulkData_.parallel_owner_rank(node) == messager::myProcNo())
        {
            nMatchedSurface1BoundaryNodes++;
        }
    }

    for (const auto& key : matchedSurface2Nodes)
    {
        stk::mesh::Entity node = bulkData_.get_entity(key);
        if (bulkData_.is_valid(node) &&
            bulkData_.parallel_owner_rank(node) == messager::myProcNo())
        {
            nMatchedSurface2BoundaryNodes++;
        }
    }

    // Reduce across all processors
    messager::sumReduce(nMatchedSurface1BoundaryNodes);
    messager::sumReduce(nMatchedSurface2BoundaryNodes);

    label nBoundaryNodesOnSurface1 = sphereBoundingBoxSurface1Vec.size();
    messager::sumReduce(nBoundaryNodesOnSurface1);

    label nBoundaryNodesOnSurface2 = sphereBoundingBoxSurface2Vec.size();
    messager::sumReduce(nBoundaryNodesOnSurface2);

    // For full overlap, all boundary nodes from both surfaces must match
    bool isOverlapping =
        (nMatchedSurface1BoundaryNodes == nBoundaryNodesOnSurface1) &&
        (nMatchedSurface2BoundaryNodes == nBoundaryNodesOnSurface2);

    return isOverlapping;
}

bool surfaceComparator::checkConformality(
    std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
        matchingNodePairVector,
    scalar conformalityCheckSearchTolerance,
    vector sepVec,
    matrix rotMat,
    bool allowPartialOverlap)
{
    // Check that both surfaces feature the same side topologies.
    // Each rank only owns a subset of sides, so we reduce topology
    // flags globally before comparing.
    const label nTopos = stk::topology::NUM_TOPOLOGIES;

    auto collectTopoFlags =
        [&](const stk::mesh::PartVector& parts, std::vector<label>& flags)
    {
        flags.assign(nTopos, 0);
        stk::mesh::BucketVector const& sideBuckets = bulkData_.get_buckets(
            metaData_.side_rank(),
            metaData_.locally_owned_part() & stk::mesh::selectUnion(parts));
        for (const stk::mesh::Bucket* bp : sideBuckets)
        {
            flags[bp->topology().value()] = 1;
        }
    };

    std::vector<label> localFlags1, localFlags2;
    collectTopoFlags(surface1Parts_, localFlags1);
    collectTopoFlags(surface2Parts_, localFlags2);

    std::vector<label> globalFlags1(nTopos), globalFlags2(nTopos);
    stk::all_reduce_max(
        messager::comm(), localFlags1.data(), globalFlags1.data(), nTopos);
    stk::all_reduce_max(
        messager::comm(), localFlags2.data(), globalFlags2.data(), nTopos);

    if (globalFlags1 != globalFlags2)
    {
        return false;
    }

    // required data structures; master/slave
    std::vector<sphereBoundingBox> sphereBoundingBoxSurface1Vec;
    std::vector<sphereBoundingBox> sphereBoundingBoxSurface2Vec;

    // Point
    Point masterCenter, slaveCenter;

    // Master: setup sphereBoundingBoxMasterVec,
    stk::mesh::BucketVector const& surface1NodeBuckets =
        bulkData_.get_buckets(stk::topology::NODE_RANK,
                              metaData_.locally_owned_part() &
                                  stk::mesh::selectUnion(surface1Parts_));

    for (stk::mesh::BucketVector::const_iterator ib =
             surface1NodeBuckets.begin();
         ib != surface1NodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        // point to data
        const scalar* coords = stk::mesh::field_data(*coordsSTKFieldPtr_, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            stk::mesh::Entity node = b[k];

            // setup ident
            theEntityKey theIdent(bulkData_.entity_key(node),
                                  messager::myProcNo());

            // define offset for all nodal fields that are of SPATIAL_DIM
            const size_t offSet = k * SPATIAL_DIM;

            // sum local coords for translation; define localCoords for bounding
            // point
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar cxj = coords[offSet + j];
                masterCenter[j] = cxj;
            }

            // create the bounding point sphere and push back
            sphereBoundingBox theSphere(
                Sphere(masterCenter, conformalityCheckSearchTolerance),
                theIdent);
            sphereBoundingBoxSurface1Vec.push_back(theSphere);
        }
    }

    // Slave: setup sphereBoundingBoxSlaveVec; translate slave onto master
    stk::mesh::BucketVector const& surface2NodeBuckets =
        bulkData_.get_buckets(stk::topology::NODE_RANK,
                              metaData_.locally_owned_part() &
                                  stk::mesh::selectUnion(surface2Parts_));

    // auxiliary
    utils::vector rotCoords;

    for (stk::mesh::BucketVector::const_iterator ib =
             surface2NodeBuckets.begin();
         ib != surface2NodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        // point to data
        const scalar* beginCoords =
            stk::mesh::field_data(*coordsSTKFieldPtr_, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            stk::mesh::Entity node = b[k];

            // setup ident
            theEntityKey theIdent(bulkData_.entity_key(node),
                                  messager::myProcNo());

            // define offset for all nodal fields that are of
            // SPATIAL_DIM
            const size_t offset = k * SPATIAL_DIM;

            // rotate coordinates
            const utils::vectorViewC coords(&beginCoords[offset]);
            rotCoords = utils::transformVector(rotMat, coords);

            // create the bounding point that is translated, then
            // push back
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar xj = rotCoords[j];
                slaveCenter[j] = xj + sepVec[j];
            }
            sphereBoundingBox theSphere(
                Sphere(slaveCenter, conformalityCheckSearchTolerance),
                theIdent);
            sphereBoundingBoxSurface2Vec.push_back(theSphere);
        }
    }

    // determine search method for this pair; default is stk_kdtree
    stk::search::SearchMethod searchMethod = stk::search::KDTREE;
    if (conformalityCheckSearchMethodName_ != "stk_kdtree")
    {
        errorMsg("searchMethod only supports stk_kdtree");
    }

    // will want to stuff product of search to a single vector
    searchKeyVector searchKeyPair;
    stk::search::coarse_search(sphereBoundingBoxSurface2Vec,
                               sphereBoundingBoxSurface1Vec,
                               searchMethod,
                               messager::comm(),
                               searchKeyPair);

    // manage ghosting: create, ghost, popoulate and destroy
    std::vector<stk::mesh::EntityProc> sendNodes;
    for (size_t i = 0, size = searchKeyPair.size(); i < size; ++i)
    {
        unsigned domainProc = searchKeyPair[i].first.proc();
        unsigned rangeProc = searchKeyPair[i].second.proc();

        if ((messager::myProcNo() != domainProc) &&
            (messager::myProcNo() == rangeProc))
        {
            stk::mesh::Entity rangeNode =
                bulkData_.get_entity(searchKeyPair[i].second.id());
            sendNodes.push_back(stk::mesh::EntityProc(rangeNode, domainProc));
        }
        else if ((messager::myProcNo() == domainProc) &&
                 (messager::myProcNo() != rangeProc))
        {
            stk::mesh::Entity domainNode =
                bulkData_.get_entity(searchKeyPair[i].first.id());
            sendNodes.push_back(stk::mesh::EntityProc(domainNode, rangeProc));
        }
    }

    addRangeNodesToSharersOfDomainNodes_(searchKeyPair, sendNodes);

    size_t numNodes = sendNodes.size();
    size_t g_numNodes = 0;
    stk::all_reduce_sum(messager::comm(), &numNodes, &g_numNodes, 1);

    stk::mesh::Ghosting* conformalGhosting{nullptr};
    if (g_numNodes > 0)
    {
        std::vector<label> ghostCommProcs;

        // check if we need to ghost
        bulkData_.modification_begin();

        // adjust ghosting name to the current surfaces' names
        conformalGhosting = &bulkData_.create_ghosting(
            "conformal_ghosting_" + surface1Parts_.back()->name() + "_" +
            surface2Parts_.back()->name());

        bulkData_.change_ghosting(*conformalGhosting, sendNodes);
        bulkData_.modification_end();

        ops::populateGhostCommProcs(
            bulkData_, *conformalGhosting, ghostCommProcs);
    }

    // now populate pair vector
    for (size_t i = 0, size = searchKeyPair.size(); i < size; ++i)
    {
        stk::mesh::Entity domainNode =
            bulkData_.get_entity(searchKeyPair[i].first.id());
        stk::mesh::Entity rangeNode =
            bulkData_.get_entity(searchKeyPair[i].second.id());

        // unique surface1 node to surface2 node pair
        std::pair<stk::mesh::Entity, stk::mesh::Entity> theFirstPair =
            std::make_pair(rangeNode, domainNode);
        matchingNodePairVector.push_back(theFirstPair);
    }

    // Count unique matched nodes on both surfaces using EntityKey sets.
    // The coarse search was called with surface2 as domain (first)
    // and surface1 as range (second).
    std::set<stk::mesh::EntityKey> ownedMatchedSurface1Keys;
    std::set<stk::mesh::EntityKey> ownedMatchedSurface2Keys;

    for (size_t i = 0, size = searchKeyPair.size(); i < size; ++i)
    {
        if (searchKeyPair[i].second.proc() == messager::myProcNo())
        {
            ownedMatchedSurface1Keys.insert(searchKeyPair[i].second.id());
        }
        if (searchKeyPair[i].first.proc() == messager::myProcNo())
        {
            ownedMatchedSurface2Keys.insert(searchKeyPair[i].first.id());
        }
    }

    label nOwnedMatchingSurface1Nodes =
        static_cast<label>(ownedMatchedSurface1Keys.size());
    label nOwnedMatchingSurface2Nodes =
        static_cast<label>(ownedMatchedSurface2Keys.size());

    messager::sumReduce(nOwnedMatchingSurface1Nodes);
    messager::sumReduce(nOwnedMatchingSurface2Nodes);

    label nNodesOnSurface1 = sphereBoundingBoxSurface1Vec.size();
    messager::sumReduce(nNodesOnSurface1);

    label nNodesOnSurface2 = sphereBoundingBoxSurface2Vec.size();
    messager::sumReduce(nNodesOnSurface2);

    label isConformal = false;
    if (allowPartialOverlap)
    {
        // The smaller surface must be fully matched, and the number of
        // matched nodes on the larger surface must equal the smaller
        // surface's total (ensuring a 1-to-1 mapping in the overlap).
        if (nNodesOnSurface1 <= nNodesOnSurface2)
        {
            isConformal = (nOwnedMatchingSurface1Nodes == nNodesOnSurface1) &&
                          (nOwnedMatchingSurface2Nodes == nNodesOnSurface1);
        }
        else
        {
            isConformal = (nOwnedMatchingSurface2Nodes == nNodesOnSurface2) &&
                          (nOwnedMatchingSurface1Nodes == nNodesOnSurface2);
        }
    }
    else
    {
        // Full overlap: all nodes on both surfaces must match 1-to-1.
        isConformal = (nNodesOnSurface1 == nNodesOnSurface2) &&
                      (nOwnedMatchingSurface1Nodes == nNodesOnSurface1) &&
                      (nOwnedMatchingSurface2Nodes == nNodesOnSurface2);
    }

    // must be destroyed: only used to communicate send-nodes
    if (g_numNodes > 0)
    {
        bulkData_.modification_begin();
        bulkData_.destroy_ghosting(*conformalGhosting);
        bulkData_.modification_end();
    }

    return isConformal;
}

vector surfaceComparator::determineSeparationVector()
{
    // now get the centroids and populate the displacement
    vector centroid_1(vector::Zero());
    vector centroid_2(vector::Zero());

    calcCentroid_(surface1Parts_, centroid_1);
    calcCentroid_(surface2Parts_, centroid_2);

    // translation vector: displacement from side2 to side1
    return centroid_1 - centroid_2;
}

scalar surfaceComparator::determineSeparationAngle(vector rotationAxis,
                                                   vector axisLocation)
{
    basis B(rotationAxis);

    auto SA = collectSamples_(surface1Parts_, B, axisLocation);
    auto SB = collectSamples_(surface2Parts_, B, axisLocation);

    // Compute global moments m_k for k=1..maxK
    std::vector<std::complex<scalar>> mA(maxK_ + 1, 0.0),
        mB(maxK_ + 1, 0.0); // 0 unused
    for (const auto& s : SA)
    {
        std::complex<scalar> z(s.x, s.y);
        std::complex<scalar> zk = 1.0;
        for (label k = 1; k <= maxK_; ++k)
        {
            zk *= z;
            mA[k] += s.w * zk;
        }
    }
    for (const auto& s : SB)
    {
        std::complex<scalar> z(s.x, s.y);
        std::complex<scalar> zk = 1.0;
        for (label k = 1; k <= maxK_; ++k)
        {
            zk *= z;
            mB[k] += s.w * zk;
        }
    }

    // Reduce across ranks
    allreduceComplexSum_(mA);
    allreduceComplexSum_(mB);

    // Try k = 1..maxK, pick the first with decent magnitude
    for (label k = 1; k <= maxK_; ++k)
    {
        scalar mag = std::abs(mA[k]) * std::abs(mB[k]);
        if (mag > SMALL)
        {
            scalar ang = std::arg(mB[k]) - std::arg(mA[k]); // in (-pi,pi]
            // wrap
            while (ang <= -M_PI)
                ang += 2 * M_PI;
            while (ang > M_PI)
                ang -= 2 * M_PI;
            return ang / k;
        }
    }
    throw std::runtime_error(
        "estimate_angle_with_axis_point: ill-conditioned (moments near zero).");
}

void surfaceComparator::readNodeXYZ_(stk::mesh::Entity node, vector& out)
{
    const scalar* xyz = stk::mesh::field_data(*coordsSTKFieldPtr_, node);
    out = vector(xyz[0], xyz[1], xyz[2]);
}

void surfaceComparator::decomposeFaceToTries_(stk::mesh::Entity face,
                                              std::vector<faceTriCentroid>& out)
{
    out.clear();
    const stk::mesh::Entity* nodes = bulkData_.begin_nodes(face);
    const unsigned n = bulkData_.num_nodes(face);
    if (n < 3)
        return;

    vector v0;
    readNodeXYZ_(nodes[0], v0);

    for (unsigned i = 1; i + 1 < n; ++i)
    {
        vector v1, v2;
        readNodeXYZ_(nodes[i], v1);
        readNodeXYZ_(nodes[i + 1], v2);

        vector a = v1 - v0, b = v2 - v0;
        scalar area = scalar(0.5) * (a.cross(b)).norm();
        if (area <= 0.0)
            continue;

        vector ctr = (v0 + v1 + v2) * scalar(1.0 / 3.0);
        out.push_back({ctr, area});
    }
}

std::vector<surfaceComparator::sample2D>
surfaceComparator::collectSamples_(const stk::mesh::PartVector& parts,
                                   const basis& B,
                                   const vector& p)
{
    std::vector<sample2D> S;
    std::vector<faceTriCentroid> tris;

    stk::mesh::Selector selOwnedFaces =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(parts);
    std::vector<stk::mesh::Entity> faces;
    stk::mesh::get_entities(
        bulkData_, stk::topology::FACE_RANK, selOwnedFaces, faces);

    S.reserve(faces.size() * 2);

    for (auto f : faces)
    {
        tris.clear();
        decomposeFaceToTries_(f, tris);
        for (const auto& t : tris)
        {
            scalar X, Y;
            B.projectToPlane(t.c, p, X, Y);
            scalar r = std::hypot(X, Y);
            if (r == 0.0)
                continue;
            scalar phi = std::atan2(Y, X);
            S.push_back({X, Y, r, phi, t.A});
        }
    }
    return S;
}

void surfaceComparator::allreduceComplexSum_(
    std::vector<std::complex<scalar>>& v)
{
    const label N = static_cast<label>(v.size());
    std::vector<scalar> send(2 * N), recv(2 * N);
    for (label i = 0; i < N; ++i)
    {
        send[2 * i] = v[i].real();
        send[2 * i + 1] = v[i].imag();
    }
    stk::all_reduce_sum(messager::comm(), send.data(), recv.data(), 2 * N);
    for (label i = 0; i < N; ++i)
        v[i] = std::complex<scalar>(recv[2 * i], recv[2 * i + 1]);
}

void surfaceComparator::calcCentroid_(const stk::mesh::PartVector& surfaceParts,
                                      vector& centroid)
{
    vector lcentroid(vector::Zero());
    scalar area = 0.0;

    stk::mesh::BucketVector const& sideBuckets = bulkData_.get_buckets(
        metaData_.side_rank(),
        metaData_.locally_owned_part() & stk::mesh::selectUnion(surfaceParts));
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // define scratch field
        std::vector<scalar> ws_coordinates(nodesPerSide * SPATIAL_DIM);
        std::vector<scalar> ws_scs_areav(numScsBip * SPATIAL_DIM);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // calculate face centroid
            stk::mesh::Entity const* sideNodeRels = bulkData_.begin_nodes(side);
            label numSideNodes = bulkData_.num_nodes(side);

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            label num_nodes = sideBucket.num_nodes(iSide);
            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];
                scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr_, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_coordinates[offSet + j] = coords[j];
                }
            }

            // compute scs integration point areavec
            scalar scs_error = 0.0;
            meFC->determinant(
                1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

            scalar sideCentroid[SPATIAL_DIM] = {0};
            scalar sideAreaMag = 0.0;

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label faceOffSet = ip * SPATIAL_DIM;

                // calc vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = ws_scs_areav[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                sideAreaMag += amag;
            }

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    sideCentroid[i] += coords[i];
                }
            }

            // side centroid
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                sideCentroid[i] /= static_cast<scalar>(numSideNodes);
            }

            // accumulate to centroid
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                lcentroid[i] += sideCentroid[i] * sideAreaMag;
            }

            // accumulate to area
            area += sideAreaMag;
        }
    }

    scalar garea = 0.0;

    // sync
    stk::all_reduce_sum(messager::comm(), &area, &garea, 1);
    stk::all_reduce_sum(
        messager::comm(), lcentroid.data(), centroid.data(), SPATIAL_DIM);

    // average
    for (label i = 0; i < SPATIAL_DIM; ++i)
    {
        centroid[i] /= garea;
    }
}

void surfaceComparator::addRangeNodesToSharersOfDomainNodes_(
    const searchKeyVector& searchKeyPair,
    std::vector<stk::mesh::EntityProc>& sendNodes)
{
    stk::CommSparse commSparse(messager::comm());

    auto packingLambda = [&]()
    {
        label theRank = messager::myProcNo();
        std::vector<label> sharingProcs;
        for (size_t i = 0; i < searchKeyPair.size(); ++i)
        {
            label domainProc = searchKeyPair[i].first.proc();
            label rangeProc = searchKeyPair[i].second.proc();

            stk::mesh::EntityKey domainKey = searchKeyPair[i].first.id();
            stk::mesh::Entity domainNode = bulkData_.get_entity(domainKey);

            if ((theRank == domainProc) &&
                bulkData_.bucket(domainNode).shared())
            {
                stk::mesh::EntityId rangeId = searchKeyPair[i].second.id().id();
                stk::CommBuffer& sbuf = commSparse.send_buffer(rangeProc);

                bulkData_.comm_shared_procs(domainKey, sharingProcs);
                if (theRank == rangeProc)
                {
                    stk::mesh::Entity rangeNode =
                        bulkData_.get_entity(stk::topology::NODE_RANK, rangeId);
                    for (label p : sharingProcs)
                    {
                        sendNodes.push_back(
                            stk::mesh::EntityProc(rangeNode, p));
                    }
                }
                else
                {
                    for (label p : sharingProcs)
                    {
                        if (p != theRank && p != rangeProc)
                        {
                            sbuf.pack(rangeId);
                            sbuf.pack(p);
                        }
                    }
                }
            }
        }
    };

    stk::pack_and_communicate(commSparse, packingLambda);

    stk::unpack_communications(commSparse,
                               [&](label p)
    {
        stk::CommBuffer& rbuf = commSparse.recv_buffer(p);
        stk::mesh::EntityId gid;
        rbuf.unpack(gid);
        label proc;
        rbuf.unpack(proc);
        sendNodes.push_back(stk::mesh::EntityProc(
            bulkData_.get_entity(stk::topology::NODE_RANK, gid), proc));
    });
}

} // namespace utils

#endif /* SPATIAL_DIM = 3 */

} // namespace accel
