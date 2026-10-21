// File       : pointInPolyhedron.h
// Created    : Fri Mar 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Point-in-polyhedron containment tests using half-space and ray
//              casting
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef POINTINPOLYHEDRON_H
#define POINTINPOLYHEDRON_H

// code
#include "types.h"

namespace accel
{

#if SPATIAL_DIM == 3

typedef stk::search::Point<scalar> Point;
typedef stk::search::IdentProc<label, label> Id;

namespace utils
{

enum class convexPolySearchAlgorithm
{
    halfSpace
};

enum class concavePolySearchAlgorithm
{
    rayCasting
};

class pointInPolyhedron
{
private:
    stk::mesh::BulkData& bulkData_;

    stk::mesh::MetaData& metaData_;

    stk::mesh::ConstPartVector& polyhedron_;

    stk::mesh::ConstPartVector& envelope_;

    stk::mesh::Field<scalar>* coordsSTKFieldPtr_;

    convexPolySearchAlgorithm convexSearchMethod_ =
        convexPolySearchAlgorithm::halfSpace;

    concavePolySearchAlgorithm concaveSearchMethod_ =
        concavePolySearchAlgorithm::rayCasting;

    scalar scale_ = 1.0;

    scalar shift_[3] = {0, 0, 0};

    scalar referencePoint_[3] = {0, 0, 0};

    scalar minCorner_[3] = {0, 0, 0};

    scalar maxCorner_[3] = {0, 0, 0};

    // absolute tolerance derived from the envelope size
    scalar geomTol_ = 0.0;

    label referenceRank_ = -1;

    bool convexPolyhedron_ = true;

    void determineEnvelopeBounds_();

    void checkIfConvexEnvelope_();

    void determineReferencePointOnEnvelope_();

    void filterConvexPolyhedron_(const std::vector<double>& scatter,
                                 std::vector<label>& inliers);

    void filterConvexPolyhedronHalfSpace_(const std::vector<double>& scatter,
                                          std::vector<label>& inliers);

    void filterConcavePolyhedron_(const std::vector<double>& scatter,
                                  std::vector<label>& inliers);

    void filterConcavePolyhedronRayCasting_(const std::vector<double>& scatter,
                                            std::vector<label>& inliers);

public:
    pointInPolyhedron(stk::mesh::ConstPartVector& polyhedron,
                      stk::mesh::ConstPartVector& envelope);

    void filter(const std::vector<double>& scatter,
                std::vector<label>& inliers);
};

scalar signedTetraVolume(const scalar* P,
                         const scalar* A,
                         const scalar* B,
                         const scalar* C);

void cross(const scalar* vec1, const scalar* vec2, scalar* result);

scalar dot(const scalar* vec1, const scalar* vec2);

void subtract(const scalar* vec1, const scalar* vec2, scalar* result);

scalar magnitude(const scalar* vec);

void unit(const scalar* vec, scalar* e);

void normalize(scalar* vec);

void direction(const scalar* A, const scalar* B, scalar* dir);

bool lineTriangleIntersection(const scalar* O,
                              const scalar* D,
                              const scalar* A,
                              const scalar* B,
                              const scalar* C,
                              const scalar tol);

bool lineQuadIntersection(const scalar* O,
                          const scalar* D,
                          const scalar* A,
                          const scalar* B,
                          const scalar* C,
                          const scalar* D_quad,
                          const scalar tol);

} // namespace utils

#else /* SPATIAL_DIM == 2 */

namespace utils
{

// 2D point-in-polygon test; points on the boundary count as inside.
class pointInPolygon
{
private:
    stk::mesh::BulkData& bulkData_;

    stk::mesh::MetaData& metaData_;

    stk::mesh::ConstPartVector& polygon_;

    stk::mesh::ConstPartVector& envelope_;

    stk::mesh::Field<scalar>* coordsSTKFieldPtr_;

    // envelope edges as x0,y0,x1,y1 per segment, replicated on every rank
    std::vector<scalar> segments_;

    // per-segment axis-aligned bounds, minx,miny,maxx,maxy, grown by geomTol_
    std::vector<scalar> segBounds_;

    // bounds of the whole envelope, grown by geomTol_
    scalar envMin_[2] = {0.0, 0.0};

    scalar envMax_[2] = {0.0, 0.0};

    scalar geomTol_ = 0.0;

    void gatherEnvelopeSegments_();

    bool onEnvelope_(const scalar* p) const;

public:
    pointInPolygon(stk::mesh::ConstPartVector& polygon,
                   stk::mesh::ConstPartVector& envelope);

    void filter(const std::vector<scalar>& scatter,
                std::vector<label>& inliers);
};

// the overset code refers to the containment test by one name in both
// dimensions
using pointInPolyhedron = pointInPolygon;

} // namespace utils

#endif /* SPATIAL_DIM == 3 */

} // namespace accel

#endif // POINTINPOLYHEDRON_H
