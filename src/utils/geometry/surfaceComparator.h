// File : surfaceComparator.h
// Created : Tue Sep 16 2025 14:15:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Comparing/checking diagnostics between two surfaces
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SURFACECOMPARATOR_H
#define SURFACECOMPARATOR_H

// code
#include "vectorUtils.h"

namespace accel
{

#if SPATIAL_DIM == 3

namespace utils
{

class surfaceComparator
{
private:
    stk::mesh::BulkData& bulkData_;

    stk::mesh::MetaData& metaData_;

    stk::mesh::Field<scalar>* coordsSTKFieldPtr_;

    stk::mesh::PartVector surface1Parts_;

    stk::mesh::PartVector surface2Parts_;

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

    void calcCentroid_(const stk::mesh::PartVector& surfaceParts,
                       vector& centroid);

    // Conformality

    typedef stk::search::IdentProc<stk::mesh::EntityKey, int> theEntityKey;
    typedef stk::search::Point<double> Point;
    typedef stk::search::Sphere<double> Sphere;
    typedef std::pair<Sphere, theEntityKey> sphereBoundingBox;

    typedef std::pair<stk::mesh::Selector, stk::mesh::Selector> selectorPair;
    typedef std::vector<std::pair<theEntityKey, theEntityKey>> searchKeyVector;

    std::string conformalityCheckSearchMethodName_ = "stk_kdtree";

    std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>
        nodePairCommunicator_;

    void addRangeNodesToSharersOfDomainNodes_(
        const searchKeyVector& searchKeyPair,
        std::vector<stk::mesh::EntityProc>& sendNodes);

public:
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

    vector determineSeparationVector();

    scalar determineSeparationAngle(vector rotationAxis, vector axisLocation);
};

} // namespace utils

#endif /* SPATIAL_DIM == 3 */

} // namespace accel

#endif // SURFACECOMPARATOR_H
