// File : pointInPolyhedron.cpp
// Created : Fri Mar 14 2025 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "pointInPolyhedron.h"
#include "messager.h"

namespace accel
{

#if SPATIAL_DIM == 3

namespace utils
{

pointInPolyhedron::pointInPolyhedron(stk::mesh::ConstPartVector& polyhedron,
                                     stk::mesh::ConstPartVector& envelope)
    : bulkData_(envelope.back()->mesh_bulk_data()),
      metaData_(envelope.back()->mesh_meta_data()), polyhedron_(polyhedron),
      envelope_(envelope), coordsSTKFieldPtr_(metaData_.get_field<scalar>(
                               stk::topology::NODE_RANK,
                               metaData_.coordinate_field_name()))
{
    // Determine envelope properties
    determineReferencePointOnEnvelope_();
    determineEnvelopeBounds_();
    checkIfConvexEnvelope_();
}

void pointInPolyhedron::filter(const std::vector<scalar>& scatter,
                               std::vector<label>& inliers)
{
    if (messager::master())
    {
        std::cout << "Started point-in-polyhedron detection process ..."
                  << std::endl;
    }

    if (convexPolyhedron_)
    {
        filterConvexPolyhedron_(scatter, inliers);
    }
    else
    {
        filterConcavePolyhedron_(scatter, inliers);
    }

    if (messager::master())
    {
        std::cout << "Done point-in-polyhedron detection process" << std::endl;
    }
}

void pointInPolyhedron::filterConvexPolyhedron_(
    const std::vector<scalar>& scatter,
    std::vector<label>& inliers)
{
    switch (convexSearchMethod_)
    {
        case convexPolySearchAlgorithm::halfSpace:
            filterConvexPolyhedronHalfSpace_(scatter, inliers);
            break;
    }
}

void pointInPolyhedron::filterConvexPolyhedronHalfSpace_(
    const std::vector<scalar>& scatter,
    std::vector<label>& inliers)
{
    // define some common selectors; select sides
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_edge1(SPATIAL_DIM), ws_edge2(SPATIAL_DIM),
        ws_normal(SPATIAL_DIM), ws_nn(SPATIAL_DIM);

    // define some common selectors; select owned nodes
    stk::mesh::Selector selSides =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(envelope_);

    for (label iPoint = 0; iPoint < inliers.size(); iPoint++)
    {
        const scalar* pointCoords = &scatter[3 * iPoint];

        scalar scaledPointCoords[SPATIAL_DIM];
        for (label i = 0; i < SPATIAL_DIM; i++)
        {
            scaledPointCoords[i] = (pointCoords[i] + shift_[i]) * scale_;
        }

        label flag = 1;

        // loop over all sides of the envelope and calculate the polyhedra
        // summed sign volume
        stk::mesh::BucketVector const& sideBuckets =
            bulkData_.get_buckets(metaData_.side_rank(), selSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            const stk::topology theTopo = sideBucket.topology();

            // extract master element
            MasterElement* meFC =
                MasterElementRepo::get_surface_master_element(theTopo);

            // extract master element specifics
            const label numNodesPerSide = meFC->nodesPerElement_;

            // set sizes
            ws_coords.resize(numNodesPerSide * SPATIAL_DIM);

            // get pointers
            scalar* p_coords = &ws_coords[0];

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                const auto& side = sideBucket[iSide];

                stk::mesh::Entity const* nodeRels = bulkData_.begin_nodes(side);

                // fill with nodal values
                for (label iNode = 0; iNode < numNodesPerSide; iNode++)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        p_coords[iNode * SPATIAL_DIM + i] =
                            (coords[i] + shift_[i]) * scale_;
                    }
                }

                ws_edge1[0] = p_coords[3 * 1 + 0] - p_coords[3 * 0 + 0];
                ws_edge1[1] = p_coords[3 * 1 + 1] - p_coords[3 * 0 + 1];
                ws_edge1[2] = p_coords[3 * 1 + 2] - p_coords[3 * 0 + 2];

                ws_edge2[0] = p_coords[3 * 2 + 0] - p_coords[3 * 0 + 0];
                ws_edge2[1] = p_coords[3 * 2 + 1] - p_coords[3 * 0 + 1];
                ws_edge2[2] = p_coords[3 * 2 + 2] - p_coords[3 * 0 + 2];

                ws_nn[0] = p_coords[3 * 0 + 0] - scaledPointCoords[0];
                ws_nn[1] = p_coords[3 * 0 + 1] - scaledPointCoords[1];
                ws_nn[2] = p_coords[3 * 0 + 2] - scaledPointCoords[2];

                // calculate normal of side
                cross(&ws_edge1[0], &ws_edge2[0], &ws_normal[0]);

                // calculate
                scalar res = dot(&ws_nn[0], &ws_normal[0]);

                if (res < -SMALL)
                {
                    flag = 0;
                    break;
                }
            }

            if (flag == 0)
                break;
        }

        stk::all_reduce_min(bulkData_.parallel(), &flag, &flag, 1);

        if (flag == 1)
        {
            inliers[iPoint] = 1;
        }
    }
}

void pointInPolyhedron::filterConcavePolyhedron_(
    const std::vector<scalar>& scatter,
    std::vector<label>& inliers)
{
    switch (concaveSearchMethod_)
    {
        case concavePolySearchAlgorithm::rayCasting:
            filterConcavePolyhedronRayCasting_(scatter, inliers);
            break;
    }
}

void pointInPolyhedron::filterConcavePolyhedronRayCasting_(
    const std::vector<scalar>& scatter,
    std::vector<label>& inliers)
{
    // define some common selectors; select sides
    std::vector<scalar> ws_coords;

    // define some common selectors; select owned nodes
    stk::mesh::Selector selSides =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(envelope_);

    for (label iPoint = 0; iPoint < inliers.size(); iPoint++)
    {
        const scalar* pointCoords = &scatter[3 * iPoint];

        scalar scaledPointCoords[SPATIAL_DIM];
        for (label i = 0; i < SPATIAL_DIM; i++)
        {
            scaledPointCoords[i] = (pointCoords[i] + shift_[i]) * scale_;
        }

        scalar dir[SPATIAL_DIM] = {0, 0, 0};

        label numIntersections = 0;

        // loop over all sides of the envelope and calculate the polyhedra
        // summed sign volume
        stk::mesh::BucketVector const& sideBuckets =
            bulkData_.get_buckets(metaData_.side_rank(), selSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            const stk::topology theTopo = sideBucket.topology();

            // extract master element
            MasterElement* meFC =
                MasterElementRepo::get_surface_master_element(theTopo);

            // extract master element specifics
            const label numNodesPerSide = meFC->nodesPerElement_;

            // set sizes
            ws_coords.resize(numNodesPerSide * SPATIAL_DIM);

            // get pointers
            scalar* p_coords = &ws_coords[0];

            if (numNodesPerSide == SPATIAL_DIM)
            {
                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    const auto& side = sideBucket[iSide];

                    stk::mesh::Entity const* nodeRels =
                        bulkData_.begin_nodes(side);

                    // fill with nodal values
                    for (label iNode = 0; iNode < numNodesPerSide; iNode++)
                    {
                        stk::mesh::Entity node = nodeRels[iNode];

                        const scalar* coords =
                            stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            p_coords[iNode * SPATIAL_DIM + i] =
                                (coords[i] + shift_[i]) * scale_;
                        }
                    }

                    // get direction from test point to origin
                    dir[0] = referencePoint_[0] - scaledPointCoords[0];
                    dir[1] = referencePoint_[1] - scaledPointCoords[1];
                    dir[2] = referencePoint_[2] - scaledPointCoords[2];

                    normalize(dir);

                    // check if test point intersects the side
                    bool intersects =
                        lineTriangleIntersection(scaledPointCoords,
                                                 dir,
                                                 &ws_coords[0],
                                                 &ws_coords[3],
                                                 &ws_coords[6]);

                    if (intersects)
                    {
                        numIntersections++;
                    }
                }
            }
            else if (numNodesPerSide == 4)
            {
                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    const auto& side = sideBucket[iSide];

                    stk::mesh::Entity const* nodeRels =
                        bulkData_.begin_nodes(side);

                    // fill with nodal values
                    for (label iNode = 0; iNode < numNodesPerSide; iNode++)
                    {
                        stk::mesh::Entity node = nodeRels[iNode];

                        const scalar* coords =
                            stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            p_coords[iNode * SPATIAL_DIM + i] =
                                (coords[i] + shift_[i]) * scale_;
                        }
                    }

                    // get direction from test point to origin
                    dir[0] = referencePoint_[0] - scaledPointCoords[0];
                    dir[1] = referencePoint_[1] - scaledPointCoords[1];
                    dir[2] = referencePoint_[2] - scaledPointCoords[2];

                    normalize(dir);

                    // check if test point intersects the side
                    bool intersects = lineQuadIntersection(scaledPointCoords,
                                                           dir,
                                                           &ws_coords[0],
                                                           &ws_coords[3],
                                                           &ws_coords[6],
                                                           &ws_coords[9]);

                    if (intersects)
                    {
                        numIntersections++;
                    }
                }
            }
        }

        stk::all_reduce_sum(
            bulkData_.parallel(), &numIntersections, &numIntersections, 1);

        // check if odd
        if (numIntersections % 2 != 0)
        {
            inliers[iPoint] = 1;
        }
    }
}

void pointInPolyhedron::determineEnvelopeBounds_()
{
    scalar maxCorner[3] = {-BIG, -BIG, -BIG};
    scalar minCorner[3] = {BIG, BIG, BIG};

    scalar gMinCorner[3], gMaxCorner[3];

    // loop over all nodes of the polyhedron
    stk::mesh::Selector selNodes =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(polyhedron_);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData_.get_buckets(stk::topology::NODE_RANK, selNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunk in bucket
        const scalar* coordsb =
            stk::mesh::field_data(*coordsSTKFieldPtr_, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            minCorner[0] = std::min(
                minCorner[0], (coordsb[3 * iNode + 0] + shift_[0]) * scale_);
            minCorner[1] = std::min(
                minCorner[0], (coordsb[3 * iNode + 1] + shift_[1]) * scale_);
            minCorner[2] = std::min(
                minCorner[0], (coordsb[3 * iNode + 2] + shift_[2]) * scale_);

            maxCorner[0] = std::max(
                maxCorner[0], (coordsb[3 * iNode + 0] + shift_[0]) * scale_);
            maxCorner[1] = std::max(
                maxCorner[0], (coordsb[3 * iNode + 1] + shift_[1]) * scale_);
            maxCorner[2] = std::max(
                maxCorner[0], (coordsb[3 * iNode + 2] + shift_[2]) * scale_);
        }
    }

    stk::all_reduce_min(bulkData_.parallel(), &minCorner[0], &gMinCorner[0], 3);
    stk::all_reduce_max(bulkData_.parallel(), &maxCorner[0], &gMaxCorner[0], 3);

    if (messager::master())
    {
        std::cout << "envelope properties: " << std::endl;
        std::cout << "\tenvelope bounds: min(" << gMinCorner[0] << " "
                  << gMinCorner[1] << " " << gMinCorner[2] << "), max("
                  << gMaxCorner[0] << " " << gMaxCorner[1] << " "
                  << gMaxCorner[2] << ")" << std::endl;
    }
}

void pointInPolyhedron::checkIfConvexEnvelope_()
{
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_conn(SPATIAL_DIM), ws_center(SPATIAL_DIM),
        ws_edge1(SPATIAL_DIM), ws_edge2(SPATIAL_DIM), ws_normal(SPATIAL_DIM);

    stk::mesh::Selector selSides =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(envelope_);

    label foundConcavity = 0;

    // loop over all sides of the envelope and calculate the polyhedra summed
    // sign volume
    stk::mesh::BucketVector const& sideBuckets =
        bulkData_.get_buckets(metaData_.side_rank(), selSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        const stk::topology theTopo = sideBucket.topology();

        // extract master element
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(theTopo);

        // extract master element specifics
        const label numNodesPerSide = meFC->nodesPerElement_;

        // set sizes
        ws_coords.resize(numNodesPerSide * SPATIAL_DIM);

        // get pointers
        scalar* p_coords = &ws_coords[0];

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            const auto& side = sideBucket[iSide];

            stk::mesh::Entity const* nodeRels = bulkData_.begin_nodes(side);

            // zero center
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                ws_center[i] = 0;
            }

            // fill with nodal values
            for (label iNode = 0; iNode < numNodesPerSide; iNode++)
            {
                stk::mesh::Entity node = nodeRels[iNode];

                const scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    p_coords[iNode * SPATIAL_DIM + i] =
                        (coords[i] + shift_[i]) * scale_;
                    ws_center[i] += p_coords[iNode * SPATIAL_DIM + i];
                }
            }

            ws_center[0] /= static_cast<scalar>(numNodesPerSide);
            ws_center[1] /= static_cast<scalar>(numNodesPerSide);
            ws_center[2] /= static_cast<scalar>(numNodesPerSide);

            ws_edge1[0] =
                p_coords[SPATIAL_DIM * 1 + 0] - p_coords[SPATIAL_DIM * 0 + 0];
            ws_edge1[1] =
                p_coords[SPATIAL_DIM * 1 + 1] - p_coords[SPATIAL_DIM * 0 + 1];
            ws_edge1[2] =
                p_coords[SPATIAL_DIM * 1 + 2] - p_coords[SPATIAL_DIM * 0 + 2];

            ws_edge2[0] =
                p_coords[SPATIAL_DIM * 2 + 0] - p_coords[SPATIAL_DIM * 0 + 0];
            ws_edge2[1] =
                p_coords[SPATIAL_DIM * 2 + 1] - p_coords[SPATIAL_DIM * 0 + 1];
            ws_edge2[2] =
                p_coords[SPATIAL_DIM * 2 + 2] - p_coords[SPATIAL_DIM * 0 + 2];

            // calculate normal of side
            cross(&ws_edge1[0], &ws_edge2[0], &ws_normal[0]);
            normalize(&ws_normal[0]);

            ws_conn[0] = ws_center[0] - referencePoint_[0];
            ws_conn[1] = ws_center[1] - referencePoint_[1];
            ws_conn[2] = ws_center[2] - referencePoint_[2];

            if (dot(&ws_conn[0], &ws_normal[0]) < 0.0)
            {
                foundConcavity = 1;
                break;
            }
        }

        if (foundConcavity == 1)
        {
            break;
        }
    }

    stk::all_reduce_max(
        bulkData_.parallel(), &foundConcavity, &foundConcavity, 1);

    if (foundConcavity == 1)
    {
        convexPolyhedron_ = false;
    }
    else
    {
        convexPolyhedron_ = true;
    }

    if (messager::master())
    {
        if (convexPolyhedron_)
        {
            std::cout
                << "\t**polyhedron is convex, utilizing half-space detection "
                   "method"
                << std::endl;
        }
        else
        {
            std::cout << "\t**polyhedron is concave, utilizing ray-casting "
                         "detection method"
                      << std::endl;
        }
    }
}

void pointInPolyhedron::determineReferencePointOnEnvelope_()
{
    // get reference point, the center of the first element in the background
    // mesh
    stk::mesh::Selector selElements =
        metaData_.locally_owned_part() & stk::mesh::selectUnion(polyhedron_);

    stk::mesh::BucketVector const& elemBuckets =
        bulkData_.get_buckets(stk::topology::ELEMENT_RANK, selElements);
    for (stk::mesh::BucketVector::const_iterator ib = elemBuckets.begin();
         ib != elemBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elemBucket = **ib;

        const stk::mesh::Bucket::size_type nElemsPerBucket = elemBucket.size();

        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elemBucket.topology());

        // extract master element specifics
        const label numNodesPerElement = meSCS->nodesPerElement_;

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElemsPerBucket;
             ++iElement)
        {
            const auto& element = elemBucket[iElement];

            stk::mesh::Entity const* nodeRels = bulkData_.begin_nodes(element);

            // fill with nodal values
            for (label iNode = 0; iNode < numNodesPerElement; iNode++)
            {
                stk::mesh::Entity node = nodeRels[iNode];

                const scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr_, node);

                referencePoint_[0] += (coords[0] + shift_[0]) * scale_;
                referencePoint_[1] += (coords[1] + shift_[1]) * scale_;
                referencePoint_[2] += (coords[2] + shift_[2]) * scale_;
            }

            referencePoint_[0] /= static_cast<scalar>(numNodesPerElement);
            referencePoint_[1] /= static_cast<scalar>(numNodesPerElement);
            referencePoint_[2] /= static_cast<scalar>(numNodesPerElement);

            referenceRank_ = bulkData_.parallel_rank();

            break;
        }

        break;
    }

    // broadcast the reference point from the proc of the highest rank to all
    // other procs
    stk::all_reduce_max(
        bulkData_.parallel(), &referenceRank_, &referenceRank_, 1);
    MPI_Bcast(&referencePoint_[0],
              SPATIAL_DIM,
              MPI_DOUBLE,
              referenceRank_,
              MPI_COMM_WORLD);

    if (messager::master())
    {
        std::cout << "\treference point in polyhedron: " << referencePoint_[0]
                  << " " << referencePoint_[1] << " " << referencePoint_[2]
                  << "; owner processor -> " << referenceRank_ << std::endl;
    }
}

scalar signedTetraVolume(const scalar* P,
                         const scalar* A,
                         const scalar* B,
                         const scalar* C)
{
    // Compute determinant of the 3x3 matrix
    scalar det =
        (A[0] - P[0]) *
            ((B[1] - P[1]) * (C[2] - P[2]) - (B[2] - P[2]) * (C[1] - P[1])) -
        (A[1] - P[1]) *
            ((B[0] - P[0]) * (C[2] - P[2]) - (B[2] - P[2]) * (C[0] - P[0])) +
        (A[2] - P[2]) *
            ((B[0] - P[0]) * (C[1] - P[1]) - (B[1] - P[1]) * (C[0] - P[0]));

    return det / 6.0;
}

void cross(const scalar* vec1, const scalar* vec2, scalar* result)
{
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

scalar dot(const scalar* vec1, const scalar* vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void subtract(const scalar* vec1, const scalar* vec2, scalar* result)
{
    result[0] = vec1[0] - vec2[0];
    result[1] = vec1[1] - vec2[1];
    result[2] = vec1[2] - vec2[2];
}

scalar magnitude(const scalar* vec)
{
    scalar mag = 0.0;
    for (label i = 0; i < SPATIAL_DIM; i++)
    {
        mag += vec[i] * vec[i];
    }
    return std::sqrt(mag);
}

void normalize(scalar* vec)
{
    scalar mag = magnitude(vec);
    vec[0] = vec[0] / mag;
    vec[1] = vec[1] / mag;
    vec[2] = vec[2] / mag;
}

scalar distance(const scalar* A, const scalar* B)
{
    return std::sqrt(std::pow(B[0] - A[0], 2.0) + std::pow(B[1] - A[1], 2.0) +
                     std::pow(B[2] - A[2], 2.0));
}

void direction(const scalar* A, const scalar* B, scalar* dir)
{
    subtract(B, A, dir);
    normalize(dir);
}

bool lineTriangleIntersection(const scalar* O,
                              const scalar* D,
                              const scalar* A,
                              const scalar* B,
                              const scalar* C)
{
    // Compute the normal of the triangle
    scalar AB[3], AC[3], N[3];
    subtract(B, A, AB);
    subtract(C, A, AC);
    cross(AB, AC, N);

    // Check if the line is parallel to the plane
    scalar denom = dot(N, D);
    if (std::abs(denom) < SMALL)
    {
        return false; // Line is parallel to the plane
    }

    // Compute the intersection parameter t
    scalar AO[3];
    subtract(A, O, AO);
    scalar t = dot(N, AO) / denom;

    // Check if the intersection is in the direction of D (t >= 0)
    if (t < 0)
    {
        return false; // Intersection is in the opposite direction of D
    }

    // Compute the intersection point
    scalar P[3] = {O[0] + t * D[0], O[1] + t * D[1], O[2] + t * D[2]};

    // if intersection point and point O are very close return false
    if (distance(P, O) < SMALL)
    {
        return false;
    }

    // Compute vectors for barycentric coordinates
    scalar v0[3], v1[3], v2[3];
    subtract(B, A, v0);
    subtract(C, A, v1);
    subtract(P, A, v2);

    // Compute dot products
    scalar d00 = dot(v0, v0);
    scalar d01 = dot(v0, v1);
    scalar d11 = dot(v1, v1);
    scalar d20 = dot(v2, v0);
    scalar d21 = dot(v2, v1);

    // Compute barycentric coordinates
    scalar denom_bary = d00 * d11 - d01 * d01;
    scalar u = (d11 * d20 - d01 * d21) / denom_bary;
    scalar v = (d00 * d21 - d01 * d20) / denom_bary;

    // Check if the point is inside the triangle
    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

bool lineQuadIntersection(const scalar* O,
                          const scalar* D,
                          const scalar* A,
                          const scalar* B,
                          const scalar* C,
                          const scalar* D_quad)
{
    // Split the quad into two triangles: ABC and ACD
    // Check intersection with the first triangle (ABC)
    if (lineTriangleIntersection(O, D, A, B, C))
    {
        return true;
    }

    // Check intersection with the second triangle (ACD)
    if (lineTriangleIntersection(O, D, A, C, D_quad))
    {
        return true;
    }

    // No intersection
    return false;
}

} // namespace utils

#endif /* SPATIAL_DIM = 3 */

} // namespace accel
