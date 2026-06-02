// File       : surfaceComparator.h
// Created    : Tue Sep 16 2025 14:15:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Comparing/checking diagnostics between two surfaces
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SURFACECOMPARATOR_H
#define SURFACECOMPARATOR_H

// code
#include "vectorUtils.h"

namespace accel
{

namespace utils
{

// Geometric relation between two surfaces. Separation vector, separation
// angle, overlap and conformality are all available in 2D and 3D. The
// 3D-only helpers below back the moment-based separation-angle estimate.
class surfaceComparator
{
private:
    stk::mesh::BulkData& bulkData_;

    stk::mesh::MetaData& metaData_;

    stk::mesh::Field<scalar>* coordsSTKFieldPtr_;

    stk::mesh::PartVector surface1Parts_;

    stk::mesh::PartVector surface2Parts_;

    // area-weighted centroid of a surface (dimension generic)
    void calcCentroid_(const stk::mesh::PartVector& surfaceParts,
                       vector& centroid);

#if SPATIAL_DIM == 3

    // Separation angle

    label maxK_ = 3;

    // sample = projected triangle centroid in plane ⟂ axis
    struct sample2D
    {
        scalar x, y;   // in (e1,e2) plane
        scalar r, phi; // radius & angle
        scalar w;      // weight (area)
    };

    // triangle fan centroid/area for a generic polygon face
    struct faceTriCentroid
    {
        vector c; // centroid of triangle
        scalar A; // area of triangle
    };

    std::vector<sample2D> collectSamples_(const stk::mesh::PartVector& parts,
                                          const basis& B,
                                          const vector& p);

    void readNodeXYZ_(stk::mesh::Entity node, vector& out);

    void decomposeFaceToTries_(stk::mesh::Entity face,
                               std::vector<faceTriCentroid>& out);

    void allreduceComplexSum_(std::vector<std::complex<scalar>>& v);

#endif /* SPATIAL_DIM == 3 */

    // Overlap / conformality (dimension generic)

    typedef stk::search::IdentProc<stk::mesh::EntityKey, int> theEntityKey;
    typedef stk::search::Point<double> Point;
    typedef stk::search::Sphere<double> Sphere;
    typedef std::pair<Sphere, theEntityKey> sphereBoundingBox;

    typedef std::pair<stk::mesh::Selector, stk::mesh::Selector> selectorPair;
    typedef std::vector<std::pair<theEntityKey, theEntityKey>> searchKeyVector;

    std::string conformalityCheckSearchMethodName_ = "stk_kdtree";

    std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>
        nodePairCommunicator_;

    // collect the boundary nodes of a surface: in 3D these are the nodes on
    // edges shared by a single face; in 2D they are the endpoint nodes of the
    // edge-chain (nodes belonging to a single side).
    void getBoundaryNodes_(const stk::mesh::PartVector& surfaceParts,
                           std::set<stk::mesh::Entity>& boundaryNodes);

    void addRangeNodesToSharersOfDomainNodes_(
        const searchKeyVector& searchKeyPair,
        std::vector<stk::mesh::EntityProc>& sendNodes);

public:
    // stable matched-pair record (global id + owner); ghosted-side handles go
    // stale once checkConformality's search ghosting is destroyed
    struct matchedPair
    {
        stk::mesh::EntityId id1, id2; // surface1 (master), surface2 (slave)
        int owner1, owner2;
    };

    surfaceComparator(stk::mesh::PartVector surface1Parts,
                      stk::mesh::PartVector surface2Parts);

    bool checkOverlap(scalar overlapCheckSearchTolerance = 1e-3,
                      vector sepVec = vector::Zero(),
                      matrix rotMat = matrix::Identity());

    bool checkConformality(
        std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
            matchingNodePairVector,
        scalar conformalityCheckSearchTolerance = 1e-8,
        vector sepVec = vector::Zero(),
        matrix rotMat = matrix::Identity(),
        bool allowPartialOverlap = false);

    // matched pairs as stable ids+owners (populated by checkConformality)
    const std::vector<matchedPair>& matchedPairs() const
    {
        return matchedPairs_;
    }

private:
    std::vector<matchedPair> matchedPairs_;

public:
    vector determineSeparationVector();

    scalar determineSeparationAngle(vector rotationAxis, vector axisLocation);
};

} // namespace utils

} // namespace accel

#endif // SURFACECOMPARATOR_H
