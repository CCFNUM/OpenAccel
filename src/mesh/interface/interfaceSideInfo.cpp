// File       : interfaceSideInfo.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

// code
#include "interfaceSideInfo.h"
#include "dgInfo.h"
#include "interface.h"
#include "macros.h"
#include "mesh.h"
#include "messager.h"
#include "zone.h"

namespace accel
{

// compare operator
struct compareGaussPoint
{
    compareGaussPoint()
    {
    }

    bool operator()(const std::pair<theKey, theKey>& p, const uint64_t i)
    {
        return (p.first.id() < i);
    }

    bool operator()(const uint64_t i, const std::pair<theKey, theKey>& p)
    {
        return (i < p.first.id());
    }
};

// Unique search keys predicate: low to high id
struct sortSearchKeysPredicate
{
    sortSearchKeysPredicate()
    {
    }

    bool operator()(const std::pair<theKey, theKey>& p1,
                    const std::pair<theKey, theKey>& p2)
    {
        return (p1.first.id() < p2.first.id());
    }
};

interfaceSideInfo::interfaceSideInfo(
    interface* interfPtr,
    bool isMasterSide,
    const stk::mesh::PartVector currentPartVec,
    const stk::mesh::PartVector opposingPartVec,
    interfaceModelOption option,
    const scalar expandBoxPercentage,
    const std::string& searchMethodName,
    const bool clipIsoParametricCoords,
    const scalar searchTolerance,
    const bool dynamicSearchTolAlg,
    const bool useShifted,
    const std::string debugName)
    : currentPartVec_(currentPartVec), opposingPartVec_(opposingPartVec),
      interfPtr_(interfPtr), name_(debugName), isMasterSide_(isMasterSide),
      expandBoxPercentage_(expandBoxPercentage),
      searchMethod_(stk::search::KDTREE),
      clipIsoParametricCoords_(clipIsoParametricCoords),
      searchTolerance_(searchTolerance),
      dynamicSearchTolAlg_(dynamicSearchTolAlg), useShifted_(useShifted),
      dataHandler_(new dataHandler), interfaceModelOption_(option)
{
    // determine search method for this pair
    if (searchMethodName != "stk_kdtree")
    {
        if (messager::master())
        {
            std::cout
                << "interfaceSideInfo::search_method only supports stk_kdtree"
                << std::endl;
        }
    }
}

interfaceSideInfo::~interfaceSideInfo()
{
    deleteDgInfo();
}

static void deleteFaceInfo(std::vector<dgInfo*>& face_info)
{
    for (dgInfo* dg : face_info)
    {
        if (dg)
        {
            delete dg;
        }
    }
    face_info.clear();
}

void interfaceSideInfo::deleteDgInfo()
{
    std::vector<std::vector<dgInfo*>>::iterator ii;
    for (ii = dgInfoVec_.begin(); ii != dgInfoVec_.end(); ++ii)
    {
        std::vector<dgInfo*>& faceDgInfoVec = (*ii);
        deleteFaceInfo(faceDgInfoVec);
    }
    dgInfoVec_.clear();
}

const interface* interfaceSideInfo::interfPtr() const
{
    return interfPtr_;
}

const zone* interfaceSideInfo::zonePtr() const
{
    label zoneIndex = isMasterSide_ ? interfPtr_->masterZoneIndex()
                                    : interfPtr_->slaveZoneIndex();
    return interfPtr_->meshPtr()->zonePtr(zoneIndex);
}

dataHandler& interfaceSideInfo::dataHandlerRef()
{
    return *dataHandler_.get();
}

const dataHandler& interfaceSideInfo::dataHandlerRef() const
{
    return *dataHandler_.get();
}

void interfaceSideInfo::setup()
{
}

void interfaceSideInfo::initialize()
{
    // construct DG info vector
    constructDgInfo();

    // construct the points and boxes required for the search
    constructBoundingPoints();
    constructBoundingBoxes();

    // perform the coarse search
    doSearch_();

    // determine elements required to be ghosted
    determineElemsToGhost();
}

void interfaceSideInfo::update()
{
    // clear all search-related info
    reset();

    // construct DG info vector
    constructDgInfo();

    // construct the points and boxes required for the search
    constructBoundingPoints();
    constructBoundingBoxes();

    // perform the coarse search
    doSearch_();

    // determine elements required to be ghosted
    determineElemsToGhost();
}

void interfaceSideInfo::reset()
{
    // clear all search related info
    boundingSphereVec_.clear();
    boundingFaceElementBoxVec_.clear();
    searchKeyPair_.clear();
    deleteDgInfo();

    // set back to false
    hasNonoverlap_ = false;
}

void interfaceSideInfo::constructDgInfo()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    stk::mesh::Selector selSides =
        (metaData.locally_owned_part() | metaData.aura_part()) &
        stk::mesh::selectUnion(currentPartVec_);

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selSides);

    // need to keep track of some sort of local id for each gauss point...
    uint64_t localGaussPointId = 0;
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;

        const stk::mesh::Bucket::size_type length = b.size();

        // extract connected element topology
        b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology currentElemTopo = parentTopo[0];

        // volume and surface master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(currentElemTopo);
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(b.topology());

        // master element-specific values
        const label numScsBip = meFC->numIntPoints_;

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // get face, global and local id
            stk::mesh::Entity side = b[k];
            uint64_t globalFaceId = bulkData.identifier(side);

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label currentFaceOrdinal = face_elem_ords[0];

            std::vector<dgInfo*> faceDgInfoVec(numScsBip);
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                dgInfo* dgInfo_ptr = new dgInfo(messager::myProcNo(),
                                                globalFaceId,
                                                localGaussPointId++,
                                                ip,
                                                side,
                                                element,
                                                currentFaceOrdinal,
                                                meFC,
                                                meSCS,
                                                currentElemTopo,
                                                searchTolerance_);
                faceDgInfoVec[ip] = dgInfo_ptr;
            }

            // push them all back
            dgInfoVec_.push_back(faceDgInfoVec);
        }
    }
}

void interfaceSideInfo::constructBoundingPoints()
{
    // Note: When searching for opposing faces, it has to be ensured that
    // for a shifted ip, the nodes on the master side do not conform with those
    // on the slave side, otherwise, only regular ip's to be used in the search
    // process, and will be later restored

    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    // hold the point location for integration points
    Point currentIpCoords;
    Point currentRegularIpCoords;

    // nodal fields to gather
    std::vector<scalar> ws_face_coordinates;

    // master element
    std::vector<scalar> ws_face_shape_function;
    std::vector<scalar> ws_regular_face_shape_function;

    // fields
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    std::vector<std::vector<dgInfo*>>::iterator ii;
    for (ii = dgInfoVec_.begin(); ii != dgInfoVec_.end(); ++ii)
    {
        std::vector<dgInfo*>& theVec = (*ii);

        //=======================================================
        // all ips on this face use a common face master element
        //            gather common operations once
        //=======================================================
        const dgInfo* firstDgInfo = theVec[0];
        MasterElement* meFC = firstDgInfo->meFCCurrent_;

        // master element-specific values
        const label numScsBip = meFC->numIntPoints_;
        const label nodesPerSide = meFC->nodesPerElement_;

        // algorithm related; face
        ws_face_coordinates.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_regular_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_face_coordinates = &ws_face_coordinates[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];
        scalar* p_regular_face_shape_function =
            &ws_regular_face_shape_function[0];

        // populate shape function
        if (useShifted_)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        // populate regular shape function
        meFC->shape_fcn(&p_regular_face_shape_function[0]);

        // gather nodal data off of face
        stk::mesh::Entity const* sideNodeRels =
            bulkData.begin_nodes(firstDgInfo->currentFace_);
        const label numSideNodes =
            bulkData.num_nodes(firstDgInfo->currentFace_);

        // sanity check on num nodes
        STK_ThrowAssert(numSideNodes == nodesPerSide);
        for (label ni = 0; ni < numSideNodes; ++ni)
        {
            stk::mesh::Entity node = sideNodeRels[ni];
            scalar* coords = stk::mesh::field_data(*coordsSTKFieldPtr, node);
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                p_face_coordinates[ni * SPATIAL_DIM + j] = coords[j];
            }
        }

        // now loop over all ips on this face
        for (size_t k = 0; k < theVec.size(); ++k)
        {
            dgInfo* dgInfo = theVec[k];

            // extract point radius; set to small if dynamic alg is not
            // activated
            const scalar pointRadius =
                dynamicSearchTolAlg_
                    ? dgInfo->nearestDistance_ * dgInfo->nearestDistanceSafety_
                    : 1.0e-16;

            // local and current ip
            const uint64_t localIp = dgInfo->localGaussPointId_;
            const label currentFaceIp = dgInfo->currentGaussPointId_;

            // compute coordinates
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                currentIpCoords[j] = 0.0;
                currentRegularIpCoords[j] = 0.0;
            }

            // interpolate to gauss point
            for (label ic = 0; ic < nodesPerSide; ++ic)
            {
                const scalar r =
                    p_face_shape_function[currentFaceIp * nodesPerSide + ic];
                const scalar rg =
                    p_regular_face_shape_function[currentFaceIp * nodesPerSide +
                                                  ic];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    currentIpCoords[j] +=
                        r * p_face_coordinates[ic * SPATIAL_DIM + j];
                    currentRegularIpCoords[j] +=
                        rg * p_face_coordinates[ic * SPATIAL_DIM + j];
                }
            }

            // extract isoparametric coords on current face from meFC
            const scalar* intgLoc =
                useShifted_ ? &meFC->intgLocShift_[0] : &meFC->intgLoc_[0];

            // copy these coordinates
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                dgInfo->currentGaussPointCoords_[j] = currentIpCoords[j];
            }

            // save face iso-parametric coordinates; extract conversion factor
            // from CVFEM to isInElement
            const scalar conversionFac = meFC->scaleToStandardIsoFac_;
            for (label j = 0; j < SPATIAL_DIM - 1; ++j)
            {
                dgInfo->currentIsoParCoords_[j] =
                    conversionFac *
                    intgLoc[currentFaceIp * (SPATIAL_DIM - 1) + j];
            }

            // setup ident for this point; use local integration point id
            stk::search::IdentProc<uint64_t, label> theIdent(
                localIp, messager::myProcNo());

            // create the bounding sphere and push back
            if (interfPtr_->isConformal())
            {
                boundingSphere theSphere(
                    Sphere(currentRegularIpCoords, pointRadius), theIdent);
                boundingSphereVec_.push_back(theSphere);
            }
            else
            {
                boundingSphere theSphere(Sphere(currentIpCoords, pointRadius),
                                         theIdent);
                boundingSphereVec_.push_back(theSphere);
            }
        }
    }
}

void interfaceSideInfo::constructBoundingBoxes()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    // specify dynamic tolerance algorithm factor
    const scalar dynamicFac = dynamicSearchTolAlg_ ? 0.0 : 1.0;

    // fields
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    // points
    Point minCorner, maxCorner;

    stk::mesh::Selector selOwned = metaData.locally_owned_part() &
                                   stk::mesh::selectUnion(opposingPartVec_);

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwned);

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // initialize max and min
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                minCorner[j] = +1.0e16;
                maxCorner[j] = -1.0e16;
            }

            // extract elem_node_relations
            stk::mesh::Entity const* face_node_rels =
                bulkData.begin_nodes(side);
            const label num_nodes = bulkData.num_nodes(side);

            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = face_node_rels[ni];

                const utils::vectorViewC coords(
                    stk::mesh::field_data(*coordsSTKFieldPtr, node));

                // Calculate min/max
                switch (interfaceModelOption_)
                {
                    case interfaceModelOption::rotationalPeriodicity:
                        {
                            utils::vector newCoords =
                                utils::transformVector(rotationMatrix_, coords);
                            newCoords += translationVector_;

                            // check max/min
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                minCorner[j] =
                                    std::min(minCorner[j], newCoords[j]);
                                maxCorner[j] =
                                    std::max(maxCorner[j], newCoords[j]);
                            }
                        }
                        break;

                    case interfaceModelOption::translationalPeriodicity:
                        {
                            utils::vector newCoords(coords);
                            newCoords += translationVector_;

                            // check max/min
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                minCorner[j] =
                                    std::min(minCorner[j], newCoords[j]);
                                maxCorner[j] =
                                    std::max(maxCorner[j], newCoords[j]);
                            }
                        }
                        break;

                    case interfaceModelOption::generalConnection:
                        {
                            // check max/min
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                minCorner[j] =
                                    std::min(minCorner[j], coords[j]);
                                maxCorner[j] =
                                    std::max(maxCorner[j], coords[j]);
                            }
                        }
                        break;
                }
            }

            // setup ident
            stk::search::IdentProc<uint64_t, label> theIdent(
                bulkData.identifier(side), messager::myProcNo());

            // expand the box by both % and search tolerance
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const scalar theMin = minCorner[i];
                const scalar theMax = maxCorner[i];
                const scalar increment =
                    expandBoxPercentage_ * (theMax - theMin) +
                    searchTolerance_ * dynamicFac;

                minCorner[i] -= increment;
                maxCorner[i] += increment;
            }

            // create the bounding point box and push back
            boundingElementBox theBox(Box(minCorner, maxCorner), theIdent);
            boundingFaceElementBoxVec_.push_back(theBox);
        }
    }
}

void interfaceSideInfo::doSearch_()
{
    stk::search::coarse_search(boundingSphereVec_,
                               boundingFaceElementBoxVec_,
                               searchMethod_,
                               messager::comm(),
                               searchKeyPair_);

    if (dynamicSearchTolAlg_)
    {
        repeatSearchIfNeeded_(boundingSphereVec_, searchKeyPair_);
    }

    // sort based on local gauss point id
    std::sort(searchKeyPair_.begin(),
              searchKeyPair_.end(),
              sortSearchKeysPredicate());
}

void interfaceSideInfo::repeatSearchIfNeeded_(
    const std::vector<boundingSphere>& boundingSphereVec,
    std::vector<std::pair<theKey, theKey>>& searchKeyPair) const
{
    unsigned num_iterations = 0;

    std::vector<boundingSphere> SphereVec(boundingSphereVec);
    deleteRangePointsFound_(SphereVec, searchKeyPair);

    label any_not_empty, not_empty = !SphereVec.empty();
    stk::all_reduce_sum(MPI_COMM_WORLD, &not_empty, &any_not_empty, 1);

    while (any_not_empty)
    {
        for (auto& ii : SphereVec)
            ii.first.set_radius(2 * ii.first.radius());
        std::vector<std::pair<theKey, theKey>> KeyPair;
        stk::search::coarse_search(SphereVec,
                                   boundingFaceElementBoxVec_,
                                   searchMethod_,
                                   MPI_COMM_WORLD,
                                   KeyPair);

        searchKeyPair.reserve(searchKeyPair.size() + KeyPair.size());
        searchKeyPair.insert(
            searchKeyPair.end(), KeyPair.begin(), KeyPair.end());

        deleteRangePointsFound_(SphereVec, searchKeyPair);
        not_empty = !SphereVec.empty();
        stk::all_reduce_sum(MPI_COMM_WORLD, &not_empty, &any_not_empty, 1);
        ++num_iterations;

        if (10 < num_iterations)
        {
            if (messager::master())
                std::cout
                    << "interfaceSideInfo::repeatSearchIfNeeded issue with "
                    << name_ << std::endl;
            if (messager::master())
                std::cout << "Increased search tolerance 10 times and still "
                             "failed to find match."
                          << std::endl;
            throw std::runtime_error(
                "Could be an internal logic error. Try turning off dynamic "
                "search tolerance algorithm...");
        }
    }
}

void interfaceSideInfo::deleteRangePointsFound_(
    std::vector<boundingSphere>& SphereVec,
    const std::vector<std::pair<theKey, theKey>>& searchKeyPair) const
{
    struct compare
    {
        bool operator()(const boundingSphere& a, const theKey& b) const
        {
            return a.second < b;
        }

        bool operator()(const theKey& a, const boundingSphere& b) const
        {
            return a < b.second;
        }

        bool operator()(const boundingSphere& a, const boundingSphere& b) const
        {
            return a.second.id() < b.second.id();
        }
    };

    if (!std::is_sorted(SphereVec.begin(), SphereVec.end(), compare()))
        std::sort(SphereVec.begin(), SphereVec.end(), compare());

    std::vector<theKey> keys_found;
    keys_found.reserve(searchKeyPair.size());
    for (const auto& ii : searchKeyPair)
    {
        keys_found.push_back(ii.first);
    }
    {
        std::sort(keys_found.begin(), keys_found.end());
        const auto it = std::unique(keys_found.begin(), keys_found.end());
        keys_found.resize(it - keys_found.begin());
    }
    std::vector<boundingSphere> difference(SphereVec.size());
    {
        const auto it = std::set_difference(SphereVec.begin(),
                                            SphereVec.end(),
                                            keys_found.begin(),
                                            keys_found.end(),
                                            difference.begin(),
                                            compare());
        difference.resize(it - difference.begin());
    }
    swap(difference, SphereVec);
}

void interfaceSideInfo::determineElemsToGhost()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    auto& elementsToGhost = interfPtr_->elemsToGhost_;

    std::vector<std::pair<theKey, theKey>>::const_iterator ii;
    for (ii = searchKeyPair_.begin(); ii != searchKeyPair_.end(); ++ii)
    {
        const uint64_t theBox = ii->second.id();
        unsigned theRank = messager::myProcNo();
        const unsigned pt_proc = ii->first.proc();
        const unsigned box_proc = ii->second.proc();

        if ((box_proc == theRank) && (pt_proc != theRank))
        {
            // Send box to pt proc

            // find the face element
            stk::mesh::Entity side =
                bulkData.get_entity(metaData.side_rank(), theBox);

            if (bulkData.is_valid(side))
            {
                // extract the connected element
                const stk::mesh::Entity* face_elem_rels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);
                stk::mesh::Entity element = face_elem_rels[0];

                // deal with elements to push back to be ghosted; downward
                // relations come for the ride...
                stk::mesh::EntityProc theElemPair(element, pt_proc);
                elementsToGhost.push_back(theElemPair);
            }
            else
            {
                // extract the connected element
                const stk::mesh::Entity* face_elem_rels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);
                stk::mesh::Entity element = face_elem_rels[0];

                messager::print("error at " + name_ + " element id: " +
                                std::to_string(bulkData.identifier(element)));
            }
        }
    }
}

void interfaceSideInfo::completeSearch()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    // dynamic algorithm requires normal distance between point and ip
    // must be always 3
    scalar bestElemIpCoords[3];

    // fields
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    std::vector<scalar> currentGaussPointCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoordsOrg(SPATIAL_DIM);

    // invert the process... Loop over dgInfoVec_ and query searchKeyPair_ for
    // this information
    std::vector<dgInfo*> exposedDgInfoVec;
    std::vector<std::vector<dgInfo*>>::iterator ii;
    for (ii = dgInfoVec_.begin(); ii != dgInfoVec_.end(); ++ii)
    {
        std::vector<dgInfo*>& theVec = (*ii);
        for (size_t k = 0; k < theVec.size(); ++k)
        {
            dgInfo* dgInfo = theVec[k];
            const uint64_t localGaussPointId = dgInfo->localGaussPointId_;

            // set initial nearestDistance and save off nearest distance under
            // dgInfo
            scalar nearestDistance = std::numeric_limits<scalar>::max();
            const scalar nearestDistanceSaved = dgInfo->nearestDistance_;

            std::pair<std::vector<std::pair<theKey, theKey>>::const_iterator,
                      std::vector<std::pair<theKey, theKey>>::const_iterator>
                p2 = std::equal_range(searchKeyPair_.begin(),
                                      searchKeyPair_.end(),
                                      localGaussPointId,
                                      compareGaussPoint());

            if (p2.first == p2.second)
            {
                // no intersections for this gauss point: exposed
                dgInfo->gaussPointExposed_ = true;
                exposedDgInfoVec.push_back(dgInfo);
            }
            else
            {
                for (std::vector<std::pair<theKey, theKey>>::const_iterator jj =
                         p2.first;
                     jj != p2.second;
                     ++jj)
                {
                    const uint64_t theBox = jj->second.id();
                    const unsigned theRank = messager::myProcNo();
                    const unsigned pt_proc = jj->first.proc();

                    // check if I own the point...
                    if (theRank == pt_proc)
                    {
                        // yes, I own the point... However, what about the face
                        // element? Who owns that?

                        // proceed as required; all elements should have already
                        // been ghosted via the coarse search
                        stk::mesh::Entity opposingFace =
                            bulkData.get_entity(metaData.side_rank(), theBox);
                        if (!(bulkData.is_valid(opposingFace)))
                            errorMsg("Invalid opposing face " +
                                     std::to_string(theBox) + " in rank " +
                                     std::to_string(theRank));

                        label opposingFaceIsGhosted =
                            bulkData.bucket(opposingFace).owned() ? 0 : 1;

                        // extract the gauss point coordinates
                        currentGaussPointCoords =
                            dgInfo->currentGaussPointCoords_;

                        // now load the face elemental nodal coords
                        stk::mesh::Entity const* sideNodeRels =
                            bulkData.begin_nodes(opposingFace);
                        label numNodes = bulkData.num_nodes(opposingFace);

                        std::vector<scalar> theOpposingFaceCoords(SPATIAL_DIM *
                                                                  numNodes);

                        // In case of a periodic interface, we must keep a copy
                        // of the real iso parametric coords of the opposing
                        // face. To do this, the real opposing face nodes must
                        // be determined
                        std::vector<scalar> theOpposingFaceCoordsOrg(
                            SPATIAL_DIM * numNodes);

                        // Apply transformation to opposing face nodes
                        switch (interfaceModelOption_)
                        {
                            case interfaceModelOption::rotationalPeriodicity:
                                {
                                    for (label ni = 0; ni < numNodes; ++ni)
                                    {
                                        stk::mesh::Entity node =
                                            sideNodeRels[ni];
                                        const utils::vectorViewC coords(
                                            stk::mesh::field_data(
                                                *coordsSTKFieldPtr, node));
                                        utils::vector new_coords =
                                            utils::transformVector(
                                                rotationMatrix_, coords);
                                        new_coords += translationVector_;

                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            const label offSet =
                                                j * numNodes + ni;
                                            theOpposingFaceCoords[offSet] =
                                                new_coords(j);

                                            theOpposingFaceCoordsOrg[offSet] =
                                                coords(j);
                                        }
                                    }
                                }
                                break;

                            case interfaceModelOption::translationalPeriodicity:
                                {
                                    for (label ni = 0; ni < numNodes; ++ni)
                                    {
                                        stk::mesh::Entity node =
                                            sideNodeRels[ni];
                                        const utils::vectorViewC coords(
                                            stk::mesh::field_data(
                                                *coordsSTKFieldPtr, node));
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            const label offSet =
                                                j * numNodes + ni;
                                            theOpposingFaceCoords[offSet] =
                                                coords(j) +
                                                translationVector_(j);

                                            theOpposingFaceCoordsOrg[offSet] =
                                                coords(j);
                                        }
                                    }
                                }
                                break;

                            case interfaceModelOption::generalConnection:
                                {
                                    for (label ni = 0; ni < numNodes; ++ni)
                                    {
                                        stk::mesh::Entity node =
                                            sideNodeRels[ni];
                                        const utils::vectorViewC coords(
                                            stk::mesh::field_data(
                                                *coordsSTKFieldPtr, node));
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            const label offSet =
                                                j * numNodes + ni;
                                            theOpposingFaceCoords[offSet] =
                                                coords(j);

                                            theOpposingFaceCoordsOrg[offSet] =
                                                coords(j);
                                        }
                                    }
                                }
                                break;
                        }

                        // extract the topo from this face element...
                        const stk::topology theOpposingFaceTopo =
                            bulkData.bucket(opposingFace).topology();
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theOpposingFaceTopo);

                        // extract the connected element to the opposing face
                        const stk::mesh::Entity* opposing_face_elem_rels =
                            bulkData.begin_elements(opposingFace);
                        STK_ThrowAssert(bulkData.num_elements(opposingFace) ==
                                        1);
                        stk::mesh::Entity opposingElement =
                            opposing_face_elem_rels[0];

                        // extract the opposing element topo and associated
                        // master element
                        const stk::topology theOpposingElementTopo =
                            bulkData.bucket(opposingElement).topology();
                        MasterElement* meSCS =
                            MasterElementRepo::get_surface_master_element(
                                theOpposingElementTopo);

                        // possible reuse
                        dgInfo->allOpposingFaceIds_.push_back(
                            bulkData.identifier(opposingFace));

                        // find distance between true current gauss point coords
                        // (the point) and the candidate bounding box
                        const scalar nearDistance =
                            meFC->isInElement(&theOpposingFaceCoords[0],
                                              &currentGaussPointCoords[0],
                                              &opposingIsoParCoords[0]);

                        // check if this is the best candidate
                        if (nearDistance < dgInfo->bestX_)
                        {
                            // save the opposing face element and master element
                            dgInfo->opposingFace_ = opposingFace;
                            dgInfo->meFCOpposing_ = meFC;

                            if (dynamicSearchTolAlg_)
                            {
                                // find the projected normal distance between
                                // point and centroid; all we need is an
                                // approximation
                                meFC->interpolatePoint(
                                    SPATIAL_DIM,
                                    &opposingIsoParCoords[0],
                                    &theOpposingFaceCoords[0],
                                    &bestElemIpCoords[0]);
                                scalar theDistance = 0.0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    scalar dxj = currentGaussPointCoords[j] -
                                                 bestElemIpCoords[j];
                                    theDistance += dxj * dxj;
                                }
                                theDistance = std::sqrt(theDistance);
                                nearestDistance =
                                    std::min(nearestDistance, theDistance);

                                // If the nearest distance between the surfaces
                                // at this point is smaller then the current
                                // distance can be reduced a bit.  Otherwise
                                // make sure the current distance is increased
                                // as needed.
                                if (nearestDistance < dgInfo->nearestDistance_)
                                {
                                    const scalar relax = 0.8;
                                    dgInfo->nearestDistance_ =
                                        relax * nearestDistanceSaved +
                                        (1.0 - relax) * nearestDistance;
                                }
                                else
                                {
                                    dgInfo->nearestDistance_ = nearestDistance;
                                }
                            }

                            // save off ordinal for opposing face
                            const stk::mesh::ConnectivityOrdinal*
                                face_elem_ords =
                                    bulkData.begin_element_ordinals(
                                        opposingFace);
                            dgInfo->opposingFaceOrdinal_ = face_elem_ords[0];

                            // save off all required opposing information
                            dgInfo->opposingElement_ = opposingElement;
                            dgInfo->meSCSOpposing_ = meSCS;
                            dgInfo->opposingElementTopo_ =
                                theOpposingElementTopo;
                            dgInfo->opposingIsoParCoords_ =
                                opposingIsoParCoords;
                            dgInfo->bestX_ = nearDistance;
                            dgInfo->opposingFaceIsGhosted_ =
                                opposingFaceIsGhosted;
                        }
                    }
                    else
                    {
                        // not this proc's issue
                    }
                }
            }
        }
    }

    // check for reuse and also provide diagnostics on sizes for opposing
    // surface set
    size_t totalOpposingFaceSize = 0;
    size_t totalDgInfoSize = 0;
    size_t maxOpposingSize = 0;
    size_t minOpposingSize = 1e6;

    for (size_t iv = 0; iv < dgInfoVec_.size(); ++iv)
    {
        std::vector<dgInfo*>& theVec = dgInfoVec_[iv];
        for (size_t k = 0; k < theVec.size(); ++k)
        {
            // extract the info object; new and old
            dgInfo* dgInfo = theVec[k];

            if (dgInfo->gaussPointExposed_)
            {
                // do nothing
            }
            else
            {
                // counts
                size_t opposingCount = dgInfo->allOpposingFaceIds_.size();
                totalDgInfoSize++;
                totalOpposingFaceSize += opposingCount;
                maxOpposingSize = std::max(maxOpposingSize, opposingCount);
                minOpposingSize = std::min(minOpposingSize, opposingCount);

                // hack ... this assumes that no geometric errors exist BODGEE
                if (!bulkData.is_valid(dgInfo->opposingFace_))
                {
                    dgInfo->gaussPointExposed_ = true;
                    exposedDgInfoVec.push_back(dgInfo);
                }
            }
        }
    }

    // global sum
    if (messager::master())
    {
        std::cout << "  > DgInfo size overview for: " << name_ << std::endl;
    }

    // set boolean of exposed dg's to true
    if (exposedDgInfoVec.size() > 0)
    {
        hasNonoverlap_ = true;
    }

    label hasNonoverlap = hasNonoverlap_ ? 1 : 0;
    label numberOfExposedDg = exposedDgInfoVec.size();
    messager::maxReduce(hasNonoverlap);
    messager::sumReduce(numberOfExposedDg);

    // finally, provide mean opposing face count
    size_t g_total[2] = {};
    size_t g_minOpposingSize;
    size_t g_maxOpposingSize;
    size_t l_total[2] = {totalDgInfoSize, totalOpposingFaceSize};

    stk::all_reduce_sum(MPI_COMM_WORLD, l_total, g_total, 2);
    stk::all_reduce_min(
        MPI_COMM_WORLD, &minOpposingSize, &g_minOpposingSize, 1);
    stk::all_reduce_max(
        MPI_COMM_WORLD, &maxOpposingSize, &g_maxOpposingSize, 1);

    if (messager::master())
    {
        std::cout << "    Min/Max/Average opposing face size: "
                  << g_minOpposingSize << "/" << g_maxOpposingSize << "/"
                  << std::fixed << std::setprecision(2)
                  << static_cast<double>(g_total[1]) / g_total[0] << std::endl
                  << std::endl;

        if (hasNonoverlap)
        {
            std::cout << "  ** The interface has an exposed portion "
                         "(non-overlap). Count of exposed dg's: "
                      << numberOfExposedDg << std::endl
                      << std::endl;
        }
    }
}

void interfaceSideInfo::determineOpposingGaussPointIds()
{
    for (auto& faceDgInfoVec : dgInfoVec_)
    {
        for (dgInfo* dg : faceDgInfoVec)
        {
            if (dg->gaussPointExposed_)
                continue;

            MasterElement* meFCOpposing = dg->meFCOpposing_;
            const label npe = meFCOpposing->nodesPerElement_;
            const std::vector<double>& intgLoc =
                useShifted_ ? meFCOpposing->intgLocShift_
                            : meFCOpposing->intgLoc_;
            const label nDimFace = static_cast<label>(intgLoc.size()) / npe;
            const std::vector<scalar>& oppIsoParCoords =
                dg->opposingIsoParCoords_;

            // Find the opposing IP whose parametric location is closest
            // to the opposingIsoParCoords (accounts for different node
            // orderings between current and opposing faces)
            label bestIp = 0;
            scalar bestDist2 = std::numeric_limits<scalar>::max();
            for (label ip = 0; ip < npe; ++ip)
            {
                scalar dist2 = 0.0;
                for (label d = 0; d < nDimFace; ++d)
                {
                    const scalar diff =
                        oppIsoParCoords[d] - intgLoc[ip * nDimFace + d];
                    dist2 += diff * diff;
                }
                if (dist2 < bestDist2)
                {
                    bestDist2 = dist2;
                    bestIp = ip;
                }
            }
            dg->opposingGaussPointId_ = bestIp;
        }
    }
}

void interfaceSideInfo::provideDiagnosis()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    std::vector<scalar> currentGaussPointCoords(SPATIAL_DIM);
    std::vector<scalar> opposingGaussPointCoords(SPATIAL_DIM);

    std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

    std::cout << std::endl;
    std::cout << "Non Conformal Alg review for surface: " << name_ << std::endl;
    std::cout << "===================================== " << std::endl;
    std::vector<std::vector<dgInfo*>>::iterator ii;
    for (ii = dgInfoVec_.begin(); ii != dgInfoVec_.end(); ++ii)
    {
        std::vector<dgInfo*>& theVec = (*ii);
        for (size_t k = 0; k < theVec.size(); ++k)
        {
            dgInfo* dgInfo = theVec[k];

            // first, dump info
            dgInfo->dumpInfo();

            // now proceed to detailed face/element current/opposing checks
            const uint64_t localGaussPointId = dgInfo->localGaussPointId_;
            const uint64_t currentGaussPointId = dgInfo->currentGaussPointId_;

            // extract current face
            stk::mesh::Entity currentFace = dgInfo->currentFace_;

            // extract the gauss point isopar/geometric coordinates for current
            currentGaussPointCoords = dgInfo->currentGaussPointCoords_;
            currentIsoParCoords = dgInfo->currentIsoParCoords_;

            // extract the master element for current; with npe
            MasterElement* meFCCurrent = dgInfo->meFCCurrent_;
            const label currentNodesPerFace = meFCCurrent->nodesPerElement_;

            // face:node relations
            stk::mesh::Entity const* current_face_node_rels =
                bulkData.begin_nodes(currentFace);
            label current_face_num_nodes = bulkData.num_nodes(currentFace);

            // gather nodal coordinates
            std::vector<scalar> currentFaceNodalCoords(SPATIAL_DIM *
                                                       currentNodesPerFace);
            for (label ni = 0; ni < current_face_num_nodes; ++ni)
            {
                stk::mesh::Entity node = current_face_node_rels[ni];
                const scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    currentFaceNodalCoords[j * currentNodesPerFace + ni] =
                        coords[j];
                }
            }

            // interpolate to current ip
            std::vector<scalar> checkCurrentFaceGaussPointCoords(SPATIAL_DIM);
            meFCCurrent->interpolatePoint(SPATIAL_DIM,
                                          &currentIsoParCoords[0],
                                          &currentFaceNodalCoords[0],
                                          &checkCurrentFaceGaussPointCoords[0]);

            // extract current element
            stk::mesh::Entity currentElement = dgInfo->currentElement_;

            // best X
            const scalar bX = dgInfo->bestX_;

            // best opposing face
            stk::mesh::Entity theBestFace = dgInfo->opposingFace_;

            // extract the gauss point isopar coordiantes for opposing
            opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

            // extract the master element for opposing; with npe
            MasterElement* meFCOpposing = dgInfo->meFCOpposing_;
            const label opposingNodesPerFace = meFCOpposing->nodesPerElement_;

            // face:node relations
            stk::mesh::Entity const* opposing_face_node_rels =
                bulkData.begin_nodes(theBestFace);
            label opposing_face_num_nodes = bulkData.num_nodes(theBestFace);

            // gather nodal coordinates
            std::vector<scalar> opposingFaceNodalCoords(SPATIAL_DIM *
                                                        opposingNodesPerFace);

            // Apply transformation to minCorner and maxCorner
            switch (interfaceModelOption_)
            {
                case interfaceModelOption::rotationalPeriodicity:
                    {
                        for (label ni = 0; ni < opposing_face_num_nodes; ++ni)
                        {
                            stk::mesh::Entity node =
                                opposing_face_node_rels[ni];

                            const utils::vectorViewC coords(
                                stk::mesh::field_data(*coordsSTKFieldPtr,
                                                      node));
                            utils::vector new_coords =
                                utils::transformVector(rotationMatrix_, coords);
                            new_coords += translationVector_;

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                opposingFaceNodalCoords
                                    [j * opposingNodesPerFace + ni] =
                                        new_coords(j);
                            }
                        }
                    }
                    break;

                case interfaceModelOption::translationalPeriodicity:
                    {
                        for (label ni = 0; ni < opposing_face_num_nodes; ++ni)
                        {
                            stk::mesh::Entity node =
                                opposing_face_node_rels[ni];
                            const utils::vectorViewC coords(
                                stk::mesh::field_data(*coordsSTKFieldPtr,
                                                      node));
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                opposingFaceNodalCoords
                                    [j * opposingNodesPerFace + ni] =
                                        coords(j) + translationVector_(j);
                            }
                        }
                    }
                    break;

                case interfaceModelOption::generalConnection:
                    {
                        for (label ni = 0; ni < opposing_face_num_nodes; ++ni)
                        {
                            stk::mesh::Entity node =
                                opposing_face_node_rels[ni];
                            const utils::vectorViewC coords(
                                stk::mesh::field_data(*coordsSTKFieldPtr,
                                                      node));
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                opposingFaceNodalCoords
                                    [j * opposingNodesPerFace + ni] = coords(j);
                            }
                        }
                    }
                    break;
            }

            // interpolate to opposing ip
            std::vector<scalar> checkOpposingFaceGaussPointCoords(SPATIAL_DIM);
            meFCOpposing->interpolatePoint(
                SPATIAL_DIM,
                &opposingIsoParCoords[0],
                &opposingFaceNodalCoords[0],
                &(checkOpposingFaceGaussPointCoords[0]));

            // global id for opposing element
            const uint64_t opElemId =
                bulkData.identifier(dgInfo->opposingElement_);

            // compute a norm between the curent nd opposing coordinate checks
            scalar distanceNorm = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                distanceNorm +=
                    std::pow(checkCurrentFaceGaussPointCoords[j] -
                                 checkOpposingFaceGaussPointCoords[j],
                             2);
            }
            distanceNorm = std::sqrt(distanceNorm);

            // provide output...
            std::cout << "Gauss Point Lid: " << localGaussPointId << " Review "
                      << std::endl;
            std::cout << "  encapsulated by nodes with Gid: (";
            for (label ni = 0; ni < current_face_num_nodes; ++ni)
            {
                stk::mesh::Entity node = current_face_node_rels[ni];
                std::cout << bulkData.identifier(node) << " ";
            }
            std::cout << ")" << std::endl;
            std::cout << "Current Gauss Point id: " << currentGaussPointId
                      << " (nearest node) "
                      << bulkData.identifier(
                             current_face_node_rels[currentGaussPointId])
                      << std::endl;

            std::cout << "  Current element Gid: "
                      << bulkData.identifier(currentElement)
                      << " (face ordinal: " << dgInfo->currentFaceOrdinal_
                      << ")" << std::endl;

            std::cout << "  has ip coordinates: ";
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                std::cout << currentGaussPointCoords[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "  The best X is: " << bX << std::endl;
            std::cout << "  Opposing element Gid: " << opElemId
                      << " (face ordinal: " << dgInfo->opposingFaceOrdinal_
                      << ")" << std::endl;
            std::cout << "  encapsulated by Gid: (";
            for (label ni = 0; ni < opposing_face_num_nodes; ++ni)
            {
                stk::mesh::Entity node = opposing_face_node_rels[ni];
                std::cout << bulkData.identifier(node) << " ";
            }
            std::cout << ")" << std::endl;
            std::cout << "  INTERNAL CHECK.... does current Gp and opposing "
                         "found Gp match coordinates? What error?"
                      << std::endl;
            std::cout << "  current and opposing ip coordinates:        "
                      << std::endl;
            for (label i = 0; i < SPATIAL_DIM; ++i)
                std::cout << "      " << i << " "
                          << checkCurrentFaceGaussPointCoords[i] << " "
                          << checkOpposingFaceGaussPointCoords[i] << std::endl;
            std::cout << "  current and opposing ip isoPar coordinates: "
                      << std::endl;
            for (label i = 0; i < SPATIAL_DIM - 1; ++i)
            {
                std::cout << "      " << i << " " << currentIsoParCoords[i]
                          << " " << opposingIsoParCoords[i] << std::endl;
            }
            std::cout << std::endl;
            std::cout << " in the end, the Error Distance Norm is: "
                      << distanceNorm << std::endl;
            std::cout << "-----------------------------------------------------"
                         "--------------"
                      << std::endl;
        }
    }
}

void interfaceSideInfo::dumpSearchResults()
{
    const stk::mesh::MetaData& metaData = interfPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData = interfPtr_->meshRef().bulkDataRef();

    for (label proci = 0; proci < messager::nProcs(); proci++)
    {
        if (proci == messager::myProcNo())
        {
            std::vector<std::pair<theKey, theKey>>::const_iterator ii;
            for (ii = searchKeyPair_.begin(); ii != searchKeyPair_.end(); ++ii)
            {
                const uint64_t theBox = ii->second.id();
                const unsigned pt_proc = ii->first.proc();
                const unsigned box_proc = ii->second.proc();

                // find the face element
                stk::mesh::Entity side =
                    bulkData.get_entity(metaData.side_rank(), theBox);

                if (bulkData.is_valid(side))
                {
                    // extract the connected element
                    const stk::mesh::Entity* face_elem_rels =
                        bulkData.begin_elements(side);
                    STK_ThrowAssert(bulkData.num_elements(side) == 1);
                    stk::mesh::Entity element = face_elem_rels[0];

                    std::cout
                        << "[" << proci << "]\t" << "ip (id: " << ii->first.id()
                        << ", proc: " << pt_proc
                        << ") -- > element id: " << bulkData.identifier(element)
                        << " (proc: " << box_proc << ")" << std::endl;
                }
            }
        }
        messager::barrier();
    }

    if (messager::master())
    {
        std::cout << "\nFinished dumping search results for interface side "
                  << name_ << std::endl;
    }
}

void interfaceSideInfo::dumpFaceToFaceResults()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();

    for (label proci = 0; proci < messager::nProcs(); proci++)
    {
        if (proci == messager::myProcNo())
        {
            std::vector<std::vector<dgInfo*>>::iterator ii;
            for (ii = dgInfoVec_.begin(); ii != dgInfoVec_.end(); ++ii)
            {
                std::vector<dgInfo*>& theVec = (*ii);
                for (size_t k = 0; k < theVec.size(); ++k)
                {
                    dgInfo* dgInfo = theVec[k];

                    size_t cid = bulkData.identifier(dgInfo->currentElement_);
                    size_t oid = bulkData.identifier(dgInfo->opposingElement_);

                    std::cout
                        << "[" << messager::myProcNo() << "]\t" << label(cid)
                        << " " << label(oid)
                        << (dgInfo->gaussPointExposed_ ? " (exposed)" : "")
                        << std::endl;
                }
            }
        }
        messager::barrier();
    }

    if (messager::master())
    {
        std::cout
            << "\nFinished dumping face-to-face results for interface side "
            << name_ << std::endl;
    }
}

size_t interfaceSideInfo::errorCheck()
{
    // check for coincident nodes via intersection of parts provided
    std::vector<stk::mesh::EntityId> coincidentNodesVec;

    const stk::mesh::MetaData& meta_data = interfPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulk_data = interfPtr_->meshRef().bulkDataRef();

    stk::mesh::Selector s_locally_owned_intersected =
        meta_data.locally_owned_part() &
        stk::mesh::selectUnion(currentPartVec_) &
        stk::mesh::selectUnion(opposingPartVec_);

    stk::mesh::BucketVector const& node_buckets =
        interfPtr_->meshRef().bulkDataRef().get_buckets(
            stk::topology::NODE_RANK, s_locally_owned_intersected);

    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
         ib != node_buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            coincidentNodesVec.push_back(bulk_data.identifier(b[k]));
        }
    }

    // report the data if problem nodes were found
    if (coincidentNodesVec.size() > 0)
    {
        std::cout << std::endl;
        std::cout << "Non Conformal Alg (P" << messager::myProcNo()
                  << ") error found on surface: " << name_ << std::endl;
        std::cout << "========================================= " << std::endl;
        for (size_t k = 0; k < coincidentNodesVec.size(); ++k)
        {
            std::cout << "coincident nodeId found: " << coincidentNodesVec[k]
                      << std::endl;
        }
    }

    return coincidentNodesVec.size();
}

void interfaceSideInfo::transformCoordinateList(std::vector<scalar>& coordsList,
                                                label npts) const
{
    assert(coordsList.size() / SPATIAL_DIM == npts);

    for (label ni = 0; ni < npts; ++ni)
    {
        const utils::vectorViewC coords(&coordsList[ni * SPATIAL_DIM]);
        utils::vector new_coords =
            utils::transformVector(rotationMatrix_, coords);
        new_coords += translationVector_;

        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            coordsList[ni * SPATIAL_DIM + j] = new_coords(j);
        }
    }
}

template <>
void interfaceSideInfo::rotateVectorList<1>(std::vector<scalar>& vectorList,
                                            label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVector<1>(std::vector<scalar>& vector) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        utils::vectorView vector(&vectorList[ni * SPATIAL_DIM]);
        utils::transformVectorInPlace(rotationMatrix_, vector);
    }
}

template <>
void interfaceSideInfo::rotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        // clone
        utils::vector newVec;
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            newVec(j) = vectorList[offset];
        }

        // rotate
        utils::transformVectorInPlace(rotationMatrix_, newVec);

        // copy back
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            vectorList[offset] = newVec(j);
        }
    }
}

template <>
void interfaceSideInfo::rotateVector<SPATIAL_DIM>(
    std::vector<scalar>& vector) const
{
    assert(vector.size() == SPATIAL_DIM);
    utils::vectorView vec(vector.data());
    utils::transformVectorInPlace(rotationMatrix_, vec);
}

template <>
void interfaceSideInfo::reverseRotateVectorList<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVector<1>(
    std::vector<scalar>& vector) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        utils::vectorView vector(&vectorList[ni * SPATIAL_DIM]);
        utils::transformVectorInPlace(rotationMatrix_.transpose(), vector);
    }
}

template <>
void interfaceSideInfo::reverseRotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        // clone
        utils::vector newVec;
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            newVec(j) = vectorList[offset];
        }

        // rotate
        utils::transformVectorInPlace(rotationMatrix_.transpose(), newVec);

        // copy back
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            vectorList[offset] = newVec(j);
        }
    }
}

template <>
void interfaceSideInfo::reverseRotateVector<SPATIAL_DIM>(
    std::vector<scalar>& vector) const
{
    assert(vector.size() == SPATIAL_DIM);
    utils::vectorView vec(vector.data());
    utils::transformVectorInPlace(rotationMatrix_.transpose(), vec);
}

} // namespace accel

#endif /* HAS_INTERFACE */
