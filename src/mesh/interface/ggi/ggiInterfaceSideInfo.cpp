// File       : ggiInterfaceSideInfo.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: GGI (General Grid Interface) implementation of interfaceSideInfo
//              using polygon clipping (scan-line rasterization or exact
//              geometric intersection).
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifdef HAS_INTERFACE

// code
#include "ggiInterfaceSideInfo.h"
#include "interface.h"
#include "ipInfo.h"
#include "macros.h"
#include "mesh.h"
#include "messager.h"
#include "zone.h"

#ifdef USE_CLIPPER2
#include <clipper2/clipper.h>
#else
#include <RastClipper/RastClipper.h>
#endif

namespace accel
{

// ---------------------------------------------------------------------------
// Local geometry helpers
// ---------------------------------------------------------------------------

namespace
{

struct LocalFrame
{
    std::array<scalar, 3> centroid;
    std::array<scalar, 3> e1;
    std::array<scalar, 3> e2;
};

// Build an orthonormal 2-D coordinate frame from the node coordinates of a
// face.  nodeCoords is stored as [n0_x, n0_y, n0_z,  n1_x, ...].
// Requires numNodes >= 3 (all valid surface faces satisfy this).
LocalFrame buildLocalFrame(const scalar* nodeCoords, label numNodes)
{
    LocalFrame f;

    // Face centroid
    f.centroid = {0.0, 0.0, 0.0};
    for (label ni = 0; ni < numNodes; ++ni)
    {
        f.centroid[0] += nodeCoords[ni * SPATIAL_DIM];
        f.centroid[1] += nodeCoords[ni * SPATIAL_DIM + 1];
        f.centroid[2] += nodeCoords[ni * SPATIAL_DIM + 2];
    }
    const scalar invN = scalar(1) / scalar(numNodes);
    f.centroid[0] *= invN;
    f.centroid[1] *= invN;
    f.centroid[2] *= invN;

    // e1 = normalize( node[1] - node[0] )
    scalar ax = nodeCoords[SPATIAL_DIM] - nodeCoords[0];
    scalar ay = nodeCoords[SPATIAL_DIM + 1] - nodeCoords[1];
    scalar az = nodeCoords[SPATIAL_DIM + 2] - nodeCoords[2];
    const scalar la = std::sqrt(ax * ax + ay * ay + az * az);
    if (la > 1.0e-30)
    {
        ax /= la;
        ay /= la;
        az /= la;
    }
    f.e1 = {ax, ay, az};

    // face normal = normalize( cross(e1,  node[2] - node[0]) )
    const scalar bx = nodeCoords[2 * SPATIAL_DIM] - nodeCoords[0];
    const scalar by = nodeCoords[2 * SPATIAL_DIM + 1] - nodeCoords[1];
    const scalar bz = nodeCoords[2 * SPATIAL_DIM + 2] - nodeCoords[2];
    scalar nx = ay * bz - az * by;
    scalar ny = az * bx - ax * bz;
    scalar nz = ax * by - ay * bx;
    const scalar ln = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (ln > 1.0e-30)
    {
        nx /= ln;
        ny /= ln;
        nz /= ln;
    }

    // e2 = cross(normal, e1)  → forms a right-handed frame with e1 and normal
    f.e2 = {ny * az - nz * ay, nz * ax - nx * az, nx * ay - ny * ax};

    return f;
}

// Project the 3-D point (px, py, pz) onto the local 2-D plane.
inline std::array<scalar, 2>
projectTo2D(scalar px, scalar py, scalar pz, const LocalFrame& f)
{
    const scalar dx = px - f.centroid[0];
    const scalar dy = py - f.centroid[1];
    const scalar dz = pz - f.centroid[2];
    return {dx * f.e1[0] + dy * f.e1[1] + dz * f.e1[2],
            dx * f.e2[0] + dy * f.e2[1] + dz * f.e2[2]};
}

// Un-project a 2-D local-frame point back to 3-D physical space.
inline std::array<scalar, 3>
unprojectTo3D(scalar u, scalar v, const LocalFrame& f)
{
    return {f.centroid[0] + u * f.e1[0] + v * f.e2[0],
            f.centroid[1] + u * f.e1[1] + v * f.e2[1],
            f.centroid[2] + u * f.e1[2] + v * f.e2[2]};
}

// Build a LocalFrame from a 4-vertex SCS polygon stored as
// array<array<scalar,3>,4>.  Delegates to buildLocalFrame() via a flat copy.
LocalFrame buildScsFrame(const std::array<std::array<scalar, 3>, 4>& v)
{
    std::array<scalar, 4 * 3> coords;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < SPATIAL_DIM; ++j)
            coords[i * SPATIAL_DIM + j] = v[i][j];
    return buildLocalFrame(coords.data(), 4);
}

// Project all 4 vertices of a SCS polygon to 2D in the given frame.
std::array<std::array<scalar, 2>, 4>
projectScsTo2D(const std::array<std::array<scalar, 3>, 4>& v,
               const LocalFrame& frame)
{
    std::array<std::array<scalar, 2>, 4> out;
    for (int i = 0; i < 4; ++i)
        out[i] = projectTo2D(v[i][0], v[i][1], v[i][2], frame);
    return out;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Sort predicate for searchKeyPair_ (ascending first-key id)
// ---------------------------------------------------------------------------

struct ggiSortSearchKeysPredicate
{
    bool operator()(const std::pair<theKey, theKey>& p1,
                    const std::pair<theKey, theKey>& p2) const
    {
        return p1.first.id() < p2.first.id();
    }
};

// ---------------------------------------------------------------------------
// Constructor / destructor
// ---------------------------------------------------------------------------

ggiInterfaceSideInfo::ggiInterfaceSideInfo(
    interface* interfPtr,
    bool isMasterSide,
    const stk::mesh::PartVector currentPartVec,
    const stk::mesh::PartVector opposingPartVec,
    interfaceModelOption option,
    label bitmapResolution,
    scalar discernibleFraction,
    bool verbose,
    bool writeDiagnosisFile,
    const std::string name)
    : interfaceSideInfo(interfPtr,
                        isMasterSide,
                        currentPartVec,
                        opposingPartVec,
                        option,
                        name),
      bitmapResolution_(bitmapResolution),
      discernibleFraction_(discernibleFraction), verbose_(verbose),
      writeDiagnosisFile_(writeDiagnosisFile)
{
    if (messager::master())
    {
        std::cout << "\t[GGI Init] Side: " << name_
                  << " | Resolution: " << bitmapResolution_
                  << " | Discernible: " << discernibleFraction_ << std::endl;
    }
}

ggiInterfaceSideInfo::~ggiInterfaceSideInfo()
{
    deleteGgiInfo_();
}

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::setup()
{
}

void ggiInterfaceSideInfo::initialize()
{
    collectFaceData_();
    constructBoundingPoints_();
    constructBoundingBoxes_();
    doSearch_();
    determineElemsToGhost();
}

void ggiInterfaceSideInfo::update()
{
    reset_();
    collectFaceData_();
    constructBoundingPoints_();
    constructBoundingBoxes_();
    doSearch_();
    determineElemsToGhost();
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::deleteGgiInfo_()
{
    // Storage is owned by the base interfaceSideInfo.
    clearIpInfo_();
}

void ggiInterfaceSideInfo::reset_()
{
    faceDataVec_.clear();
    boundingSphereVec_.clear();
    boundingFaceElementBoxVec_.clear();
    searchKeyPair_.clear();
    deleteGgiInfo_();
    hasNonoverlap_ = false;
}

// ---------------------------------------------------------------------------
// collectFaceData_
//
// Iterates over all faces in currentPartVec_ and fills faceDataVec_ with
// per-face topology, master-element pointers, node coordinates, and
// IP centroid coordinates (interpolated via shape functions).
// Also pre-sizes ipInfoVec_ to match the number of collected faces.
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::collectFaceData_()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    const auto* coordsField = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    std::vector<stk::topology> parentTopo;
    std::vector<scalar> ws_shpFcn;

    stk::mesh::Selector selSides =
        (metaData.locally_owned_part() | metaData.aura_part()) &
        stk::mesh::selectUnion(currentPartVec_);

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selSides);

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        const stk::topology elemTopo = parentTopo[0];

        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(elemTopo);
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(b.topology());

        const label numIps = meFC->numIntPoints_;
        const label numNodes = meFC->nodesPerElement_;

        ws_shpFcn.resize(numIps * numNodes);
        meFC->shape_fcn(ws_shpFcn.data());

        std::array<scalar, 3> fcCentroid;

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            stk::mesh::Entity side = b[k];

            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            FaceData_ fd;
            fd.face = side;
            fd.element = element;
            fd.faceOrdinal = faceOrdinal;
            fd.globalFaceId = bulkData.identifier(side);
            fd.faceTopo = b.topology();
            fd.elemTopo = elemTopo;
            fd.meFC = meFC;
            fd.meSCS = meSCS;
            fd.numIps = numIps;
            fd.numNodes = numNodes;
            fd.nodeCoords.resize(numNodes * SPATIAL_DIM);
            fd.ipCoords.assign(numIps * SPATIAL_DIM, 0.0);

            // Gather node coordinates
            const stk::mesh::Entity* nodeRels = bulkData.begin_nodes(side);
            STK_ThrowAssert(bulkData.num_nodes(side) == numNodes);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                const scalar* c =
                    stk::mesh::field_data(*coordsField, nodeRels[ni]);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    fd.nodeCoords[ni * SPATIAL_DIM + j] = c[j];
            }

            // Interpolate IP centroid coordinates via shape functions
            for (label ip = 0; ip < numIps; ++ip)
            {
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    const scalar shp = ws_shpFcn[ip * numNodes + ni];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                        fd.ipCoords[ip * SPATIAL_DIM + j] +=
                            shp * fd.nodeCoords[ni * SPATIAL_DIM + j];
                }
            }

            // Compute SCS polygon vertices for each IP.
            // For any n-node polygon face numIps == numNodes (one SCS per
            // node). Each SCS is a quadrilateral whose vertices are, in
            // draw-order (node, fwd_edge, centroid, bwd_edge):
            //   [0] = face node ip
            //   [1] = midpoint(node_ip, node_{(ip+1)%numNodes})   (forward
            //   edge) [2] = face centroid [3] =
            //   midpoint(node_{(ip-1+numNodes)%numNodes}, node_ip) (backward
            //   edge)
            fd.scsVerts.resize(numIps);

            // Face centroid
            fcCentroid = {0.0, 0.0, 0.0};
            for (label ni = 0; ni < numNodes; ++ni)
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    fcCentroid[j] += fd.nodeCoords[ni * SPATIAL_DIM + j];
            for (label j = 0; j < SPATIAL_DIM; ++j)
                fcCentroid[j] /= static_cast<scalar>(numNodes);

            for (label ip = 0; ip < numIps; ++ip)
            {
                const label ipNext = (ip + 1) % numNodes;
                const label ipPrev = (ip - 1 + numNodes) % numNodes;

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // V0: face node ip
                    fd.scsVerts[ip][0][j] = fd.nodeCoords[ip * SPATIAL_DIM + j];
                    // V1: forward edge midpoint
                    fd.scsVerts[ip][1][j] =
                        0.5 * (fd.nodeCoords[ip * SPATIAL_DIM + j] +
                               fd.nodeCoords[ipNext * SPATIAL_DIM + j]);
                    // V2: face centroid
                    fd.scsVerts[ip][2][j] = fcCentroid[j];
                    // V3: backward edge midpoint
                    fd.scsVerts[ip][3][j] =
                        0.5 * (fd.nodeCoords[ipPrev * SPATIAL_DIM + j] +
                               fd.nodeCoords[ip * SPATIAL_DIM + j]);
                }
            }

            faceDataVec_.push_back(std::move(fd));
        }
    }

    // Pre-size ipInfoVec_ (base storage) so completeSearch() can index
    // directly
    ipInfoVec_.resize(faceDataVec_.size());
}

// ---------------------------------------------------------------------------
// constructBoundingPoints_
//
// Creates one circumscribed bounding sphere per current face, keyed by
// globalFaceId.  Using face-level spheres (rather than per-IP spheres)
// ensures that all opposing faces that could potentially overlap with any
// part of the current face are captured in the coarse search.
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::constructBoundingPoints_()
{
    for (const FaceData_& fd : faceDataVec_)
    {
        // Face centroid
        std::array<scalar, 3> centroid = {0.0, 0.0, 0.0};
        for (label ni = 0; ni < fd.numNodes; ++ni)
        {
            centroid[0] += fd.nodeCoords[ni * SPATIAL_DIM];
            centroid[1] += fd.nodeCoords[ni * SPATIAL_DIM + 1];
            centroid[2] += fd.nodeCoords[ni * SPATIAL_DIM + 2];
        }
        const scalar invN = scalar(1) / scalar(fd.numNodes);
        centroid[0] *= invN;
        centroid[1] *= invN;
        centroid[2] *= invN;

        // Circumscribed sphere radius = max node distance from centroid
        scalar maxDist = 0.0;
        for (label ni = 0; ni < fd.numNodes; ++ni)
        {
            const scalar dx = fd.nodeCoords[ni * SPATIAL_DIM] - centroid[0];
            const scalar dy = fd.nodeCoords[ni * SPATIAL_DIM + 1] - centroid[1];
            const scalar dz = fd.nodeCoords[ni * SPATIAL_DIM + 2] - centroid[2];
            maxDist = std::max(maxDist, std::sqrt(dx * dx + dy * dy + dz * dz));
        }

        Point centPt;
        centPt[0] = centroid[0];
        centPt[1] = centroid[1];
        centPt[2] = centroid[2];

        stk::search::IdentProc<uint64_t, label> theIdent(fd.globalFaceId,
                                                         messager::myProcNo());

        boundingSphereVec_.emplace_back(Sphere(centPt, maxDist + 1.0e-10),
                                        theIdent);
    }
}

// ---------------------------------------------------------------------------
// constructBoundingBoxes_
//
// Creates one axis-aligned bounding box per locally-owned opposing face.
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::constructBoundingBoxes_()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    const auto* coordsField = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    const scalar expansion = 1.0e-3;

    stk::mesh::Selector selOwned = metaData.locally_owned_part() &
                                   stk::mesh::selectUnion(opposingPartVec_);

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwned);

    Point minCorner, maxCorner;

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;
        const stk::mesh::Bucket::size_type nSides = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSides; ++iSide)
        {
            stk::mesh::Entity side = sideBucket[iSide];

            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                minCorner[j] = +1.0e16;
                maxCorner[j] = -1.0e16;
            }

            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(side);
            const label numNodes = bulkData.num_nodes(side);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                const utils::vectorViewC coords(
                    stk::mesh::field_data(*coordsField, nodeRels[ni]));

                switch (interfaceModelOption_)
                {
                    case interfaceModelOption::rotationalPeriodicity:
                        {
                            utils::vector newCoords =
                                utils::transformVector(rotationMatrix_, coords);
                            newCoords += translationVector_;

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

            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar span = maxCorner[j] - minCorner[j];
                minCorner[j] -= expansion * span + 1.0e-10;
                maxCorner[j] += expansion * span + 1.0e-10;
            }

            stk::search::IdentProc<uint64_t, label> theIdent(
                bulkData.identifier(side), messager::myProcNo());

            boundingFaceElementBoxVec_.emplace_back(Box(minCorner, maxCorner),
                                                    theIdent);
        }
    }
}

// ---------------------------------------------------------------------------
// doSearch_
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::doSearch_()
{
    stk::search::coarse_search(boundingSphereVec_,
                               boundingFaceElementBoxVec_,
                               searchMethod_,
                               messager::comm(),
                               searchKeyPair_);

    std::sort(searchKeyPair_.begin(),
              searchKeyPair_.end(),
              ggiSortSearchKeysPredicate());
}

// ---------------------------------------------------------------------------
// determineElemsToGhost
//
// For each (sphere, box) pair in the coarse-search result: if the opposing
// face box is owned by the current rank but the current-face sphere belongs
// to a different rank, ghost the opposing element to that rank so that
// completeSearch() can access its coordinates during rasterization.
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::determineElemsToGhost()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    auto& elementsToGhost = interfPtr_->elemsToGhost_;

    const unsigned myRank = messager::myProcNo();

    for (const auto& kp : searchKeyPair_)
    {
        const unsigned pt_proc = kp.first.proc();
        const unsigned box_proc = kp.second.proc();

        if ((box_proc == myRank) && (pt_proc != myRank))
        {
            stk::mesh::Entity side =
                bulkData.get_entity(metaData.side_rank(), kp.second.id());

            if (bulkData.is_valid(side))
            {
                STK_ThrowAssert(bulkData.num_elements(side) == 1);
                stk::mesh::Entity element = bulkData.begin_elements(side)[0];
                elementsToGhost.emplace_back(element,
                                             static_cast<int>(pt_proc));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// completeSearch
//
// Scan-line rasterization for area-fraction computation.
// After ghosting, all candidate opposing faces are locally available.
//
// For each current SCS (one per IP per face):
//   1. Build a local 2-D coordinate frame from the SCS polygon normal.
//   2. Init depth-buffer: BITCOL=0 (background), BITDEP=-1.0.
//   3. Draw current SCS at depth=1.0 unconditionally (ZFGEQU).
//   4. For each candidate opposing face and each of its SCS polygons, draw
//      at depth=0 only where current SCS was already drawn (ZFLESS: write
//      iff 0.0 < BITDEP).  Each opposing SCS gets a unique colour index.
//   5. Count interior pixels (4-neighbor check) per colour → area fractions.
//   6. Create one ggiInfo per (current SCS, opposing SCS) pair whose
//      area fraction exceeds discernibleFraction_.
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::completeSearch()
{
    const auto& bulkData = interfPtr_->meshRef().bulkDataRef();
    const auto& metaData = interfPtr_->meshRef().metaDataRef();

    const auto* coordsField = metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        interfPtr_->meshRef().getCoordinateFieldName());

    const int myRank = static_cast<int>(messager::myProcNo());

    // Build face-level candidate map: globalFaceId(current) → [opposing global
    // IDs]
    std::unordered_map<uint64_t, std::vector<uint64_t>> faceToCandidates;
    for (const auto& kp : searchKeyPair_)
    {
        if (static_cast<int>(kp.first.proc()) != myRank)
            continue;
        faceToCandidates[kp.first.id()].push_back(kp.second.id());
    }

    const label NBIT = bitmapResolution_;

    label exposedCount = 0;

    for (label fIdx = 0; fIdx < static_cast<label>(faceDataVec_.size()); ++fIdx)
    {
        const FaceData_& fd = faceDataVec_[fIdx];

        // ----- Collect candidate opposing faces -----

        const auto candIt = faceToCandidates.find(fd.globalFaceId);
        if (candIt == faceToCandidates.end())
        {
            // No candidates — mark all IPs as exposed
            for (label ip = 0; ip < fd.numIps; ++ip)
            {
                auto* ggiInfo = new accel::ggiInfo(myRank,
                                                   fd.globalFaceId,
                                                   ip,
                                                   fd.face,
                                                   fd.element,
                                                   fd.faceOrdinal,
                                                   fd.meFC,
                                                   fd.meSCS,
                                                   fd.elemTopo,
                                                   stk::mesh::Entity(),
                                                   stk::mesh::Entity(),
                                                   -1,
                                                   -1,
                                                   stk::topology(),
                                                   nullptr,
                                                   nullptr,
                                                   0.0,
                                                   0);
                ggiInfo->isExposed_ = true;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    ggiInfo->currentGaussPointCoords_[j] =
                        fd.ipCoords[ip * SPATIAL_DIM + j];
                const label faceDim = SPATIAL_DIM - 1;
                const scalar convFac = fd.meFC->scaleToStandardIsoFac_;
                for (label j = 0; j < faceDim; ++j)
                    ggiInfo->currentIsoParCoords_[j] =
                        convFac * fd.meFC->intgLoc_[ip * faceDim + j];
                ipInfoVec_[fIdx].push_back(ggiInfo);
                ++exposedCount;
            }
            hasNonoverlap_ = true;
            continue;
        }

        const std::vector<uint64_t>& candIds = candIt->second;
        const label numCands = static_cast<label>(candIds.size());

        // ----- Build opposing-face node data -----

        struct OppFaceData
        {
            stk::mesh::Entity face;
            stk::mesh::Entity element;
            label faceOrdinal;
            stk::topology faceTopo;
            stk::topology elemTopo;
            MasterElement* meFC;
            MasterElement* meSCS;
            label isGhosted;
            label numNodes;
            label numIps;
            std::vector<scalar> nodeCoords; // [numNodes × SPATIAL_DIM]
            bool valid;
        };

        std::vector<OppFaceData> oppData(numCands);

        for (label oi = 0; oi < numCands; ++oi)
        {
            stk::mesh::Entity oppFace =
                bulkData.get_entity(metaData.side_rank(), candIds[oi]);
            oppData[oi].valid = bulkData.is_valid(oppFace);
            if (!oppData[oi].valid)
                continue;

            stk::mesh::Entity oppElem = bulkData.begin_elements(oppFace)[0];

            oppData[oi].face = oppFace;
            oppData[oi].element = oppElem;
            oppData[oi].faceOrdinal =
                bulkData.begin_element_ordinals(oppFace)[0];
            oppData[oi].faceTopo = bulkData.bucket(oppFace).topology();
            oppData[oi].elemTopo = bulkData.bucket(oppElem).topology();
            oppData[oi].meFC = MasterElementRepo::get_surface_master_element(
                oppData[oi].faceTopo);
            oppData[oi].meSCS = MasterElementRepo::get_surface_master_element(
                oppData[oi].elemTopo);
            oppData[oi].isGhosted = bulkData.bucket(oppFace).owned() ? 0 : 1;

            oppData[oi].numNodes = bulkData.num_nodes(oppFace);
            oppData[oi].numIps = oppData[oi].meFC->numIntPoints_;
            oppData[oi].nodeCoords.resize(oppData[oi].numNodes * SPATIAL_DIM);

            const stk::mesh::Entity* oppNodes = bulkData.begin_nodes(oppFace);
            for (label ni = 0; ni < oppData[oi].numNodes; ++ni)
            {
                const scalar* c =
                    stk::mesh::field_data(*coordsField, oppNodes[ni]);

                switch (interfaceModelOption_)
                {
                    case interfaceModelOption::rotationalPeriodicity:
                        {
                            const utils::vectorViewC rawCoords(c);
                            utils::vector newCoords = utils::transformVector(
                                rotationMatrix_, rawCoords);
                            newCoords += translationVector_;

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                                oppData[oi].nodeCoords[ni * SPATIAL_DIM + j] =
                                    newCoords[j];
                        }
                        break;

                    case interfaceModelOption::translationalPeriodicity:
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                                oppData[oi].nodeCoords[ni * SPATIAL_DIM + j] =
                                    c[j] + translationVector_(j);
                        }
                        break;

                    case interfaceModelOption::generalConnection:
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                                oppData[oi].nodeCoords[ni * SPATIAL_DIM + j] =
                                    c[j];
                        }
                        break;
                }
            }
        }

        // ----- Build IP-index offset table for colour encoding -----
        // oppIpOffset[oi] = sum of numIps for all valid opposing faces j < oi.
        // oppIpOffset[numCands] = total number of opposing-SCS colour slots.
        // Using a prefix-sum avoids any assumption about a fixed
        // max-IPs-per-face and therefore works for arbitrary polygon faces.
        std::vector<label> oppIpOffset(numCands + 1, 0);
        for (label oi = 0; oi < numCands; ++oi)
            oppIpOffset[oi + 1] =
                oppIpOffset[oi] + (oppData[oi].valid ? oppData[oi].numIps : 0);

        // ----- Pre-compute opposing SCS 3-D vertices (IP-independent) -----
        // These are the same for every current-side IP; only the 2-D
        // projection changes.  Stored in a flat array indexed by
        // oppIpOffset[oi] + oip.

        const label totalOppScs = oppIpOffset[numCands];
        std::vector<std::array<std::array<scalar, 3>, 4>> oppScsVerts3D(
            totalOppScs);

        for (label oi = 0; oi < numCands; ++oi)
        {
            if (!oppData[oi].valid)
                continue;

            const label numOppNodes = oppData[oi].numNodes;
            const label numOppIps = oppData[oi].numIps;

            // Opposing face centroid
            std::array<scalar, 3> oppCentroid = {0.0, 0.0, 0.0};
            for (label ni = 0; ni < numOppNodes; ++ni)
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    oppCentroid[j] +=
                        oppData[oi].nodeCoords[ni * SPATIAL_DIM + j];
            for (label j = 0; j < SPATIAL_DIM; ++j)
                oppCentroid[j] /= static_cast<scalar>(numOppNodes);

            for (label oip = 0; oip < numOppIps; ++oip)
            {
                const label cidx = oppIpOffset[oi] + oip;
                const label oipNext = (oip + 1) % numOppNodes;
                const label oipPrev = (oip - 1 + numOppNodes) % numOppNodes;

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    oppScsVerts3D[cidx][0][j] =
                        oppData[oi].nodeCoords[oip * SPATIAL_DIM + j];
                    oppScsVerts3D[cidx][1][j] =
                        0.5 *
                        (oppData[oi].nodeCoords[oip * SPATIAL_DIM + j] +
                         oppData[oi].nodeCoords[oipNext * SPATIAL_DIM + j]);
                    oppScsVerts3D[cidx][2][j] = oppCentroid[j];
                    oppScsVerts3D[cidx][3][j] =
                        0.5 *
                        (oppData[oi].nodeCoords[oipPrev * SPATIAL_DIM + j] +
                         oppData[oi].nodeCoords[oip * SPATIAL_DIM + j]);
                }
            }
        }

        // ----- Per-IP rasterization via adapter -----

        for (label ip = 0; ip < fd.numIps; ++ip)
        {
            // Build local 2D frame from current SCS polygon
            const LocalFrame frame = buildScsFrame(fd.scsVerts[ip]);

            // Project current SCS vertices to 2D
            const auto curScs2D = projectScsTo2D(fd.scsVerts[ip], frame);

            // Build subject and clip polygons, dispatch to backend
#ifdef USE_CLIPPER2
            Clipper2Lib::PathD subjectPoly;
            subjectPoly.reserve(4);
            for (int vi = 0; vi < 4; ++vi)
                subjectPoly.emplace_back(curScs2D[vi][0], curScs2D[vi][1]);

            const double subjectArea = std::abs(Clipper2Lib::Area(subjectPoly));

            std::vector<double> fractions(totalOppScs, 0.0);
            std::vector<std::array<scalar, 2>> csCentroids2D(totalOppScs,
                                                             {0.0, 0.0});
            if (subjectArea > 1.0e-30)
            {
                for (label cidx = 0; cidx < totalOppScs; ++cidx)
                {
                    const auto oppScs2D =
                        projectScsTo2D(oppScsVerts3D[cidx], frame);
                    Clipper2Lib::PathD clip;
                    clip.reserve(4);
                    for (int vi = 0; vi < 4; ++vi)
                        clip.emplace_back(oppScs2D[vi][0], oppScs2D[vi][1]);

                    Clipper2Lib::PathsD subjects{subjectPoly};
                    Clipper2Lib::PathsD clips{clip};
                    Clipper2Lib::PathsD result = Clipper2Lib::Intersect(
                        subjects, clips, Clipper2Lib::FillRule::NonZero, 8);

                    double overlapArea = 0.0;
                    double cx = 0.0, cy = 0.0;
                    for (const auto& path : result)
                    {
                        const double pArea = std::abs(Clipper2Lib::Area(path));
                        overlapArea += pArea;

                        // Polygon centroid (shoelace-based)
                        const int np = static_cast<int>(path.size());
                        double pcx = 0.0, pcy = 0.0;
                        for (int pi = 0; pi < np; ++pi)
                        {
                            const int pj = (pi + 1) % np;
                            const double cross = path[pi].x * path[pj].y -
                                                 path[pj].x * path[pi].y;
                            pcx += (path[pi].x + path[pj].x) * cross;
                            pcy += (path[pi].y + path[pj].y) * cross;
                        }
                        if (pArea > 1.0e-30)
                        {
                            const double inv6A =
                                1.0 / (6.0 * Clipper2Lib::Area(path));
                            cx += pArea * pcx * inv6A;
                            cy += pArea * pcy * inv6A;
                        }
                    }
                    fractions[cidx] = overlapArea / subjectArea;
                    if (overlapArea > 1.0e-30)
                    {
                        csCentroids2D[cidx][0] = cx / overlapArea;
                        csCentroids2D[cidx][1] = cy / overlapArea;
                    }
                }
            }
#else
            RastClipper::PathD subjectPoly;
            subjectPoly.reserve(4);
            for (int vi = 0; vi < 4; ++vi)
                subjectPoly.push_back({curScs2D[vi][0], curScs2D[vi][1]});

            RastClipper::PathsD clipPolys(totalOppScs);
            for (label cidx = 0; cidx < totalOppScs; ++cidx)
            {
                const auto oppScs2D =
                    projectScsTo2D(oppScsVerts3D[cidx], frame);
                clipPolys[cidx].reserve(4);
                for (int vi = 0; vi < 4; ++vi)
                    clipPolys[cidx].push_back(
                        {oppScs2D[vi][0], oppScs2D[vi][1]});
            }

            const auto clipResult = RastClipper::IntersectionFractions(
                subjectPoly, clipPolys, NBIT);
            const auto& fractions = clipResult.fractions;
            std::vector<std::array<scalar, 2>> csCentroids2D(totalOppScs);
            for (label cidx = 0; cidx < totalOppScs; ++cidx)
                csCentroids2D[cidx] = {clipResult.centroids[cidx].x,
                                       clipResult.centroids[cidx].y};
#endif

            // ----- Build ggiInfo from computed fractions -----

            bool foundConnection = false;

            for (label oi = 0; oi < numCands; ++oi)
            {
                if (!oppData[oi].valid)
                    continue;
                const label numOppIps = oppData[oi].numIps;

                for (label oip = 0; oip < numOppIps; ++oip)
                {
                    const label cidx = oppIpOffset[oi] + oip;
                    const scalar frac = static_cast<scalar>(fractions[cidx]);
                    if (frac <= discernibleFraction_)
                        continue;

                    const OppFaceData& ofd = oppData[oi];
                    auto* ggiInfo = new accel::ggiInfo(myRank,
                                                       fd.globalFaceId,
                                                       ip,
                                                       fd.face,
                                                       fd.element,
                                                       fd.faceOrdinal,
                                                       fd.meFC,
                                                       fd.meSCS,
                                                       fd.elemTopo,
                                                       ofd.face,
                                                       ofd.element,
                                                       ofd.faceOrdinal,
                                                       oip,
                                                       ofd.elemTopo,
                                                       ofd.meFC,
                                                       ofd.meSCS,
                                                       frac,
                                                       ofd.isGhosted);

                    // CS centroid: un-project from 2D to 3D
                    {
                        const auto cen3D = unprojectTo3D(csCentroids2D[cidx][0],
                                                         csCentroids2D[cidx][1],
                                                         frame);
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                            ggiInfo->currentGaussPointCoords_[j] = cen3D[j];
                    }

                    // store iso-parametric coordinates
                    // (convert from CVFEM to standard -1:1 range)
                    const label faceDim = SPATIAL_DIM - 1;
                    const scalar cConvFac = fd.meFC->scaleToStandardIsoFac_;
                    for (label j = 0; j < faceDim; ++j)
                        ggiInfo->currentIsoParCoords_[j] =
                            cConvFac * fd.meFC->intgLoc_[ip * faceDim + j];
                    const scalar oConvFac = ofd.meFC->scaleToStandardIsoFac_;
                    for (label j = 0; j < faceDim; ++j)
                        ggiInfo->opposingIsoParCoords_[j] =
                            oConvFac * ofd.meFC->intgLoc_[oip * faceDim + j];

                    ipInfoVec_[fIdx].push_back(ggiInfo);
                    foundConnection = true;
                }
            }

            if (!foundConnection)
            {
                auto* ggiInfo = new accel::ggiInfo(myRank,
                                                   fd.globalFaceId,
                                                   ip,
                                                   fd.face,
                                                   fd.element,
                                                   fd.faceOrdinal,
                                                   fd.meFC,
                                                   fd.meSCS,
                                                   fd.elemTopo,
                                                   stk::mesh::Entity(),
                                                   stk::mesh::Entity(),
                                                   -1,
                                                   -1,
                                                   stk::topology(),
                                                   nullptr,
                                                   nullptr,
                                                   0.0,
                                                   0);
                ggiInfo->isExposed_ = true;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    ggiInfo->currentGaussPointCoords_[j] =
                        fd.ipCoords[ip * SPATIAL_DIM + j];
                const label faceDim = SPATIAL_DIM - 1;
                const scalar convFac = fd.meFC->scaleToStandardIsoFac_;
                for (label j = 0; j < faceDim; ++j)
                    ggiInfo->currentIsoParCoords_[j] =
                        convFac * fd.meFC->intgLoc_[ip * faceDim + j];
                ipInfoVec_[fIdx].push_back(ggiInfo);
                ++exposedCount;
                hasNonoverlap_ = true;
            }
        } // end IP loop
    } // end face loop

    // -----------------------------------------------------------------------
    // Post-loop statistics (gated on verbose_)
    // -----------------------------------------------------------------------

    if (!verbose_)
        return;

    label localTotalFaces = static_cast<label>(faceDataVec_.size());
    label localTotalIPs = 0;
    label localCoveredIPs = 0;
    label localConnections = 0;
    double localMinFrac = 1.0;
    double localMaxFrac = 0.0;
    double localSumFrac = 0.0;
    double localMaxConservErr = 0.0;
    double localSumConservErr = 0.0;

    for (const auto& faceVec : ipInfoVec_)
    {
        std::map<label, std::vector<const ggiInfo*>> ipMap;
        for (const ipInfo* ip : faceVec)
        {
            const ggiInfo* ggiInfo = static_cast<const accel::ggiInfo*>(ip);
            ipMap[ggiInfo->currentGaussPointId_].push_back(ggiInfo);
        }

        for (const auto& kv : ipMap)
        {
            ++localTotalIPs;
            const bool exposed = kv.second.front()->isExposed_;
            if (!exposed)
            {
                ++localCoveredIPs;
                double ipFracSum = 0.0;
                for (const ggiInfo* ggiInfo : kv.second)
                {
                    ++localConnections;
                    const double f =
                        static_cast<double>(ggiInfo->areaFraction_);
                    localMinFrac = std::min(localMinFrac, f);
                    localMaxFrac = std::max(localMaxFrac, f);
                    localSumFrac += f;
                    ipFracSum += f;
                }
                const double err = std::abs(1.0 - ipFracSum);
                localMaxConservErr = std::max(localMaxConservErr, err);
                localSumConservErr += err;
            }
        }
    }

    label g_hasNonoverlap = hasNonoverlap_ ? 1 : 0;
    label g_exposedCount = exposedCount;
    label g_totalFaces = localTotalFaces;
    label g_totalIPs = localTotalIPs;
    label g_coveredIPs = localCoveredIPs;
    label g_connections = localConnections;
    double g_minFrac = localMinFrac;
    double g_maxFrac = localMaxFrac;
    double g_sumFrac = localSumFrac;
    double g_maxConservErr = localMaxConservErr;
    double g_sumConservErr = localSumConservErr;

    messager::maxReduce(g_hasNonoverlap);
    messager::sumReduce(g_exposedCount);
    messager::sumReduce(g_totalFaces);
    messager::sumReduce(g_totalIPs);
    messager::sumReduce(g_coveredIPs);
    messager::sumReduce(g_connections);
    messager::minReduce(g_minFrac);
    messager::maxReduce(g_maxFrac);
    messager::sumReduce(g_sumFrac);
    messager::maxReduce(g_maxConservErr);
    messager::sumReduce(g_sumConservErr);

    if (messager::master())
    {
        const double avgFrac =
            g_connections > 0 ? g_sumFrac / static_cast<double>(g_connections)
                              : 0.0;
        const double avgConservErr =
            g_coveredIPs > 0
                ? g_sumConservErr / static_cast<double>(g_coveredIPs)
                : 0.0;

        std::cout << "\n  GGI completeSearch: " << name_ << "\n"
                  << "    faces processed    : " << g_totalFaces << "\n"
                  << "    total IPs (SCS)    : " << g_totalIPs << "\n"
                  << "    covered IPs        : " << g_coveredIPs << "\n"
                  << "    exposed IPs        : " << g_exposedCount << "\n"
                  << "    CS connections     : " << g_connections << "\n";

        if (g_connections > 0)
        {
            std::cout << std::fixed << std::setprecision(6)
                      << "    areaFrac min       : " << g_minFrac << "\n"
                      << "    areaFrac max       : " << g_maxFrac << "\n"
                      << "    areaFrac avg       : " << avgFrac << "\n"
                      << "    |1-sumFrac| max    : " << g_maxConservErr << "\n"
                      << "    |1-sumFrac| avg    : " << avgConservErr << "\n";
        }

        if (g_hasNonoverlap)
        {
            std::cout << "    *** WARNING: " << g_exposedCount
                      << " exposed (unmatched) IPs detected\n";
        }
        std::cout << std::endl;
    }

    // Detailed per-rank diagnostic files (independent of verbose_)
    if (writeDiagnosisFile_)
        provideDiagnosis();
}

// ---------------------------------------------------------------------------
// provideDiagnosis
// ---------------------------------------------------------------------------

void ggiInterfaceSideInfo::provideDiagnosis()
{
    const int rank = static_cast<int>(messager::myProcNo());

    // Per-rank diagnostic file
    const std::string fname =
        "ggi_diag_" + name_ + "_rank" + std::to_string(rank) + ".txt";
    std::ofstream fout(fname);

    auto& out = fout; // write to file; mirror key stats to stdout

    out << "GGI Diagnosis  surface: " << name_ << "  rank: " << rank << "\n";
    out << std::string(60, '=') << "\n\n";

    label totalIPs = 0;
    label exposedIPs = 0;
    scalar maxConservError = 0.0;
    scalar sumConservError = 0.0;

    for (label fi = 0; fi < static_cast<label>(ipInfoVec_.size()); ++fi)
    {
        const auto& faceVec = ipInfoVec_[fi];
        if (faceVec.empty())
            continue;

        const uint64_t faceId = faceVec.front()->globalFaceId_;
        out << "Face " << faceId << "  (" << faceVec.size()
            << " CS connections)\n";

        // Group by currentGaussPointId to compute per-IP sum of areaFractions
        std::map<label, std::vector<const ggiInfo*>> ipMap;
        for (const ipInfo* ip : faceVec)
        {
            const ggiInfo* ggiInfo = static_cast<const accel::ggiInfo*>(ip);
            ipMap[ggiInfo->currentGaussPointId_].push_back(ggiInfo);
        }

        for (const auto& kv : ipMap)
        {
            const label ipId = kv.first;
            const auto& ipEntries = kv.second;
            const ggiInfo* gi0 = ipEntries.front();

            const bool exposed = gi0->isExposed_;
            scalar fracSum = 0.0;
            for (const ggiInfo* ggiInfo : ipEntries)
                fracSum += ggiInfo->areaFraction_;

            const scalar conservErr = exposed ? 0.0 : std::abs(1.0 - fracSum);

            out << "  IP " << ipId << "  coords=(" << std::fixed
                << std::setprecision(6) << gi0->currentGaussPointCoords_[0]
                << ", " << gi0->currentGaussPointCoords_[1] << ", "
                << gi0->currentGaussPointCoords_[2] << ")"
                << "  nConnections=" << ipEntries.size()
                << "  sumFrac=" << fracSum;
            if (exposed)
                out << "  [EXPOSED]";
            else if (conservErr > 1.0e-4)
                out << "  *** CONSERVATION ERROR=" << conservErr;
            out << "\n";

            // detail each CS connection
            for (const ggiInfo* ggiInfo : ipEntries)
            {
                out << "    -> oppFace="
                    << interfPtr_->meshRef().bulkDataRef().identifier(
                           ggiInfo->opposingFace_)
                    << "  frac=" << ggiInfo->areaFraction_
                    << (ggiInfo->opposingFaceIsGhosted_ ? "  [ghosted]" : "")
                    << "\n";
            }

            // accumulate stats
            ++totalIPs;
            if (exposed)
                ++exposedIPs;
            else
            {
                sumConservError += conservErr;
                maxConservError = std::max(maxConservError, conservErr);
            }
        }
        out << "\n";
    }

    const label coveredIPs = totalIPs - exposedIPs;
    out << std::string(60, '-') << "\n";
    out << "Summary\n";
    out << "  Total IPs       : " << totalIPs << "\n";
    out << "  Covered IPs     : " << coveredIPs << "\n";
    out << "  Exposed IPs     : " << exposedIPs << "\n";
    if (coveredIPs > 0)
    {
        out << "  Max |1-sumFrac| : " << maxConservError << "\n";
        out << "  Avg |1-sumFrac| : "
            << sumConservError / static_cast<scalar>(coveredIPs) << "\n";
    }

    fout.close();

    // Mirror summary to stdout on each rank
    std::cout << "\nGGI Diagnosis surface: " << name_ << "  rank=" << rank
              << "  totalIPs=" << totalIPs << "  exposedIPs=" << exposedIPs
              << "  maxConservErr=" << maxConservError << "  (details in "
              << fname << ")\n";
}

// ---------------------------------------------------------------------------
// errorCheck
// ---------------------------------------------------------------------------

size_t ggiInterfaceSideInfo::errorCheck()
{
    std::vector<stk::mesh::EntityId> coincidentNodesVec;

    const stk::mesh::MetaData& metaData = interfPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData = interfPtr_->meshRef().bulkDataRef();

    stk::mesh::Selector s_coinc = metaData.locally_owned_part() &
                                  stk::mesh::selectUnion(currentPartVec_) &
                                  stk::mesh::selectUnion(opposingPartVec_);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, s_coinc);

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        for (stk::mesh::Bucket::size_type k = 0; k < b.size(); ++k)
            coincidentNodesVec.push_back(bulkData.identifier(b[k]));
    }

    if (!coincidentNodesVec.empty())
    {
        std::cout << "\nGGI (P" << messager::myProcNo()
                  << ") error found on surface: " << name_ << std::endl;
        std::cout << "========================================= " << std::endl;
        for (const auto& nid : coincidentNodesVec)
            std::cout << "coincident nodeId found: " << nid << std::endl;
    }

    return coincidentNodesVec.size();
}

} // namespace accel

#endif /* HAS_INTERFACE */
