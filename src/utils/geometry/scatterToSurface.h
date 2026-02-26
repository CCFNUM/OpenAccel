// File       : scatterToSurface.h
// Created    : Thu Sep 25 2025 14:15:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Interpolate data from scatter of points
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SCATTERTOSURFACE_H
#define SCATTERTOSURFACE_H

// code
#include "vectorUtils.h"

// external
#include "nanoflann.hpp"

namespace accel
{
namespace utils
{

class scatterToSurface
{
private:
    stk::mesh::BulkData& bulkData_;

    stk::mesh::MetaData& metaData_;

    stk::mesh::Field<scalar>* coordsSTKFieldPtr_;

    stk::mesh::PartVector surfaceParts_;

    // Inverse-Distancec-Weighting parameters

    // common value is 2
    scalar distPow_ = 2.0;

    // stencil of interpolation: a value of 1 means pick the closest
    label donorCount_ = 4;

    // Point cloud adaptor
    class PointCloud
    {
    public:
        PointCloud() = default;

        explicit PointCloud(label nPoints) : pts(nPoints)
        {
        }

        struct Point
        {
            scalar x, y, z;
        };

        std::vector<Point> pts;

        inline size_t kdtree_get_point_count() const
        {
            return pts.size();
        }

        inline scalar kdtree_get_pt(const size_t idx, label dim) const
        {
            if (dim == 0)
                return pts[idx].x;
            if (dim == 1)
                return pts[idx].y;
            return pts[idx].z;
        }

        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const
        {
            return false;
        }
    };

    // KD-tree (keep DistanceType=scalar, SPATIAL_DIM, and expose IndexType)
    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<scalar, PointCloud>,
        PointCloud,
        SPATIAL_DIM>;

    using IndexT =
        typename my_kd_tree_t::IndexType; // robust wrt nanoflann build

    std::unique_ptr<PointCloud> cloudPtr_ = nullptr;
    std::unique_ptr<my_kd_tree_t> indexPtr_ = nullptr;

    // Per-query buffers
    std::vector<IndexT> retIndexes_;
    std::vector<scalar> outDistsSqr_;

    // ResultSet matches DistanceType=scalar, IndexType=IndexT
    std::unique_ptr<nanoflann::KNNResultSet<scalar, IndexT>> resultSetPtr_;

public:
    scatterToSurface(std::vector<scalar>& scatterPoints,
                     stk::mesh::PartVector surfaceParts,
                     scalar distancePowerParameter = 2.0,
                     label donorPointsCount = 4);

    void interpolateToField(std::vector<scalar>& scatterValues,
                            stk::mesh::Field<scalar>* fieldPtr);
};

} // namespace utils
} // namespace accel

#endif // SCATTERTOSURFACE_H
