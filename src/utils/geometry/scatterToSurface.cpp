// File       : scatterToSurface.cpp
// Created    : Thu Sep 25 2025 14:15:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "scatterToSurface.h"

namespace accel
{
namespace utils
{

scatterToSurface::scatterToSurface(std::vector<scalar>& scatterPoints,
                                   stk::mesh::PartVector surfaceParts,
                                   scalar distancePowerParameter,
                                   label donorPointsCount)
    : bulkData_(surfaceParts.back()->mesh_bulk_data()),
      metaData_(surfaceParts.back()->mesh_meta_data()),
      coordsSTKFieldPtr_(
          metaData_.get_field<scalar>(stk::topology::NODE_RANK,
                                      metaData_.coordinate_field_name())),
      surfaceParts_(std::move(surfaceParts)), distPow_(distancePowerParameter),
      donorCount_(donorPointsCount)
{
    // ensure donorCount_ is greater than 1
    if (donorCount_ < 1)
    {
        errorMsg("donor_count must be at least 1");
    }

    const label nPoints =
        static_cast<label>(scatterPoints.size() / SPATIAL_DIM);

    // Build point cloud from flat array
    cloudPtr_ = std::make_unique<PointCloud>(nPoints);
    for (label i = 0; i < nPoints; ++i)
    {
        cloudPtr_->pts[i].x = scatterPoints[0 * nPoints + i];
        cloudPtr_->pts[i].y = scatterPoints[1 * nPoints + i];
#if SPATIAL_DIM == 3
        cloudPtr_->pts[i].z = scatterPoints[2 * nPoints + i];
#elif SPATIAL_DIM == 2
        cloudPtr_->pts[i].z = scalar(0);
#endif
    }

    // KD-tree
    indexPtr_ = std::make_unique<my_kd_tree_t>(
        SPATIAL_DIM,
        *(cloudPtr_.get()),
        nanoflann::KDTreeSingleIndexAdaptorParams(10));
    indexPtr_->buildIndex();

    // Clamp k to available points
    donorCount_ = std::max<label>(1, std::min<label>(donorCount_, nPoints));

    // Buffers + ResultSet (NOTE: do NOT init() once forever; re-init per query)
    retIndexes_.assign(donorCount_, IndexT(0));
    outDistsSqr_.assign(donorCount_, scalar(0));

    resultSetPtr_ =
        std::make_unique<nanoflann::KNNResultSet<scalar, IndexT>>(donorCount_);
}

void scatterToSurface::interpolateToField(std::vector<scalar>& scatterValues,
                                          stk::mesh::Field<scalar>* fieldPtr)
{
    // -- Sanity checks
    if (fieldPtr->entity_rank() != stk::topology::NODE_RANK)
        errorMsg("Field must be of node rank");

    const auto& nodeBuckets = bulkData_.get_buckets(
        stk::topology::NODE_RANK,
        metaData_.universal_part() & stk::mesh::selectUnion(surfaceParts_));

    // Components per node (flat layout)
    const label N = fieldPtr->max_size();

    // Number of scatter points
    const label nPoints = static_cast<label>(cloudPtr_->pts.size());

    // scatterValues is SoA: [comp0: nPoints][comp1: nPoints]...
    if (static_cast<label>(scatterValues.size()) != N * nPoints)
        errorMsg("scatterValues size must be N * nPoints (SoA layout)");

    // IDW parameters
    const scalar eps = static_cast<scalar>(1e-12); // distance floor
    const scalar p = distPow_;                     // typical: 2.0

    for (const stk::mesh::Bucket* pb : nodeBuckets)
    {
        const stk::mesh::Bucket& b = *pb;
        const auto length = b.size();

        const scalar* coords = stk::mesh::field_data(*coordsSTKFieldPtr_, b);
        scalar* values = stk::mesh::field_data(*fieldPtr, b);

        for (stk::mesh::Bucket::size_type kNode = 0; kNode < length; ++kNode)
        {
            // Query k nearest donors around this node
            resultSetPtr_->init(retIndexes_.data(), outDistsSqr_.data());
            indexPtr_->findNeighbors(
                *resultSetPtr_,
                &coords[kNode * SPATIAL_DIM],
                nanoflann::SearchParameters(/*checks*/ 10));

            const label k = static_cast<label>(donorCount_);

            // Check exact hit among returned neighbors (d^2 ~ 0)
            label exactIdx = -1;
            for (label r = 0; r < k; ++r)
            {
                if (outDistsSqr_[r] <= eps * eps)
                {
                    exactIdx = static_cast<label>(retIndexes_[r]);
                    break;
                }
            }

            for (label comp = 0; comp < N; ++comp)
            {
                if (exactIdx >= 0)
                {
                    // Perfect (or near-perfect) coincidence -> copy value
                    values[kNode * N + comp] =
                        scatterValues[comp * nPoints + exactIdx];
                    continue;
                }

                // Plain IDW over the k neighbors
                scalar num = 0.0;
                scalar den = 0.0;

                for (label r = 0; r < k; ++r)
                {
                    const IndexT idx = retIndexes_[r];
                    const scalar d2 = outDistsSqr_[r];
                    // Guard against zero distance
                    const scalar d = std::sqrt(std::max(d2, eps * eps));
                    const scalar w = static_cast<scalar>(1.0) / std::pow(d, p);

                    num += w * scatterValues[comp * nPoints + idx];
                    den += w;
                }

                if (den > 0)
                {
                    values[kNode * N + comp] = num / den;
                }
                else
                {
                    // Ultra-degenerate fallback: simple average of neighbors
                    scalar acc = 0.0;
                    for (label r = 0; r < k; ++r)
                    {
                        const IndexT idx = retIndexes_[r];
                        acc += scatterValues[comp * nPoints + idx];
                    }
                    values[kNode * N + comp] = (k > 0) ? acc / k : 0.0;
                }
            }
        }
    }
}

#if 0
void scatterToSurface::interpolateToField(std::vector<scalar>& scatterValues,
                                          stk::mesh::Field<scalar>* fieldPtr)
{
    // -- Sanity checks
    if (fieldPtr->entity_rank() != stk::topology::NODE_RANK)
    {
        errorMsg("Field must be of node rank");
    }

    const auto& nodeBuckets = bulkData_.get_buckets(
        stk::topology::NODE_RANK,
        metaData_.locally_owned_part() & stk::mesh::selectUnion(surfaceParts_));

    // Number of components per node (flat layout)
    const label N = fieldPtr->max_size();

    // Number of scatter points
    const label nPoints = static_cast<label>(cloudPtr_->pts.size());

    // scatterValues is assumed SoA: [comp0: nPoints][comp1: nPoints]...
    if (static_cast<label>(scatterValues.size()) != N * nPoints)
    {
        errorMsg("scatterValues size must be N * nPoints (SoA layout)");
    }

    // Build a local tangent frame using ONLY the k nearest centroids
    // themselves. We pick two neighbor vectors with the largest |cross| area to
    // define a stable normal.
    auto tangentFrameFromKNN_ = [&](const scalar* nodeXYZ,
                                    const std::vector<IndexT>& idxs,
                                    vector& t1,
                                    vector& t2,
                                    vector& n)
    {
    // Node position
#if SPATIAL_DIM == 3
        vector Xn({nodeXYZ[0], nodeXYZ[1], nodeXYZ[2]});
#elif SPATIAL_DIM == 2
        vector Xn({nodeXYZ[0], nodeXYZ[1]});
#endif

        const size_t k = idxs.size();
        std::vector<vector> d;
        d.reserve(k);
        for (size_t r = 0; r < k; ++r)
        {
            const auto& P = cloudPtr_->pts[idxs[r]];
#if SPATIAL_DIM == 3
            vector Xc({P.x, P.y, P.z});
#elif SPATIAL_DIM == 2
            vector Xc({P.x, P.y});
#endif
            d.push_back(Xc - Xn);
        }

#if SPATIAL_DIM == 3
        // Pick pair with max |cross| to define robust normal
        scalar bestA = 0.0;
        vector bestN = vector::Zero();
        label iBest = -1, jBest = -1;
        for (size_t i = 0; i < k; ++i)
            for (size_t j = i + 1; j < k; ++j)
            {
                vector c = cross(d[i], d[j]);
                scalar A = norm(c);
                if (A > bestA)
                {
                    bestA = A;
                    bestN = c;
                    iBest = label(i);
                    jBest = label(j);
                }
            }

        if (bestA > 0)
        {
            n = unit(bestN);
            vector base =
                (norm(d[iBest]) >= norm(d[jBest])) ? d[iBest] : d[jBest];
            base = base - n * dot(n, base); // remove normal component
            const scalar bl = norm(base);
            if (bl > 0)
                t1 = base / bl;
            else
            {
                vector tmp =
                    (std::fabs(n[2]) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
                t1 = unit(cross(tmp, n));
            }
            t2 = cross(n, t1);
        }
        else
        {
            // nearly collinear; synthesize a frame
            vector dir = (k > 0 ? unit(d[0]) : vector({1, 0, 0}));
            vector tmp =
                (std::fabs(dir[2]) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
            n = unit(cross(dir, tmp));
            t1 = unit(cross(n, dir));
            t2 = cross(n, t1);
        }
#elif SPATIAL_DIM == 2
        // 2D: canonical frame
        n = vector({0, 0});
        t1 = vector({1, 0});
        t2 = vector({0, 1});
#endif
    };

    auto affinePlaneAtOrigin_ = [](const std::vector<scalar>& xi,
                                   const std::vector<scalar>& eta,
                                   const std::vector<scalar>& u,
                                   scalar& u_hat) -> bool
    {
        const label k = static_cast<label>(xi.size());
        label I = -1, J = -1, K = -1;
        scalar bestA2 = 0;

        auto area2 = [&](label i, label j, label m) -> scalar
        {
            return std::abs((xi[j] - xi[i]) * (eta[m] - eta[i]) -
                            (xi[m] - xi[i]) * (eta[j] - eta[i]));
        };

        for (label i = 0; i < k; i++)
            for (label j = i + 1; j < k; j++)
                for (label m = j + 1; m < k; m++)
                {
                    scalar A2 = area2(i, j, m);
                    if (A2 > bestA2)
                    {
                        bestA2 = A2;
                        I = i;
                        J = j;
                        K = m;
                    }
                }

        if (bestA2 <= 1e-14)
            return false; // nearly collinear

        auto det3 = [](scalar M[3][3]) -> scalar
        {
            return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
                   M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                   M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        };

        scalar M[3][3] = {
            {1, xi[I], eta[I]}, {1, xi[J], eta[J]}, {1, xi[K], eta[K]}};
        scalar rhs[3] = {u[I], u[J], u[K]};
        scalar D = det3(M);
        if (std::abs(D) <= 1e-18)
            return false;

        // Only need 'a' (value at origin), so replace col 0 with rhs:
        scalar Ma[3][3] = {{rhs[0], M[0][1], M[0][2]},
                           {rhs[1], M[1][1], M[1][2]},
                           {rhs[2], M[2][1], M[2][2]}};
        u_hat = det3(Ma) / D; // a
        return true;
    };

    auto maxAngularGap_ = [](const std::vector<scalar>& xi,
                             const std::vector<scalar>& eta) -> scalar
    {
        const scalar PI = std::acos(scalar(-1));
        const label k = static_cast<label>(xi.size());
        if (k < 2)
            return 2 * PI;

        std::vector<scalar> ang(k);
        for (label i = 0; i < k; i++)
            ang[i] = std::atan2(eta[i], xi[i]);
        std::sort(ang.begin(), ang.end());
        scalar maxGap = 0;
        for (label i = 0; i < k - 1; i++)
            maxGap = std::max(maxGap, ang[i + 1] - ang[i]);
        maxGap = std::max(maxGap, (ang.front() + 2 * PI) - ang.back());
        return maxGap;
    };

    for (const stk::mesh::Bucket* pb : nodeBuckets)
    {
        const stk::mesh::Bucket& b = *pb;
        const auto length = b.size();

        const scalar* coords = stk::mesh::field_data(*coordsSTKFieldPtr_, b);
        scalar* values = stk::mesh::field_data(*fieldPtr, b);

        for (stk::mesh::Bucket::size_type kNode = 0; kNode < length; ++kNode)
        {
            // Query exactly donorCount_ neighbors for this node
            resultSetPtr_->init(retIndexes_.data(), outDistsSqr_.data());
            indexPtr_->findNeighbors(*resultSetPtr_,
                                     &coords[kNode * SPATIAL_DIM],
                                     nanoflann::SearchParameters(10));

            // Build tangent frame from THESE k centroids (no faces)
            vector t1, t2, n;
            tangentFrameFromKNN_(
                &coords[kNode * SPATIAL_DIM], retIndexes_, t1, t2, n);

            // Node position
#if SPATIAL_DIM == 3
            vector Xn({coords[kNode * 3 + 0],
                       coords[kNode * 3 + 1],
                       coords[kNode * 3 + 2]});
#elif SPATIAL_DIM == 2
            vector Xn({coords[kNode * 2 + 0], coords[kNode * 2 + 1]});
#endif

            // Precompute planar coordinates and check exact-hit
            const label k = static_cast<label>(donorCount_);
            std::vector<scalar> xi(k), eta(k);
            std::vector<IndexT> idxOrder(k);
            label hitIdx = -1;
            scalar minRp = std::numeric_limits<scalar>::max();

            for (label r = 0; r < k; ++r)
            {
                const IndexT idx = retIndexes_[r];
                idxOrder[r] = idx;

                // Exact (or extremely close) 3D coincidence?
                if (outDistsSqr_[r] <= scalar(1e-24))
                {
                    hitIdx = static_cast<label>(idx);
                    // still fill xi/eta for completeness
                }

#if SPATIAL_DIM == 3
                vector Xc({cloudPtr_->pts[idx].x,
                           cloudPtr_->pts[idx].y,
                           cloudPtr_->pts[idx].z});
#elif SPATIAL_DIM == 2
                vector Xc({cloudPtr_->pts[idx].x, cloudPtr_->pts[idx].y});
#endif
                auto d = Xc - Xn;
                xi[r] = dot(d, t1);
                eta[r] = dot(d, t2);

                const scalar rp = std::sqrt(xi[r] * xi[r] + eta[r] * eta[r]);
                if (rp < minRp)
                    minRp = rp;
            }

            // Geometric stencil quality (one-sided if big angular gap)
            const scalar gap = maxAngularGap_(xi, eta);

            // For each component, compute value
            for (label comp = 0; comp < N; ++comp)
            {
                // If exact hit, just copy that centroid's component
                if (hitIdx >= 0)
                {
                    values[kNode * N + comp] =
                        scatterValues[comp * nPoints + hitIdx];
                    continue;
                }

                // Collect component values for our k donors
                std::vector<scalar> uu(k);
                for (label r = 0; r < k; ++r)
                {
                    const IndexT idx = idxOrder[r];
                    uu[r] = scatterValues[comp * nPoints + idx];
                }

                // Try affine (proper extrapolation) at origin
                scalar u_aff = 0.0;
                const bool okAff = affinePlaneAtOrigin_(xi, eta, uu, u_aff);

                // Planar Shepard (convex averaging)
                const scalar cap = std::max<scalar>(1e-12, 0.25 * minRp);
                scalar num = 0.0, den = 0.0;
                for (label r = 0; r < k; ++r)
                {
                    scalar rp = std::sqrt(xi[r] * xi[r] + eta[r] * eta[r]);
                    rp = std::max(rp, cap);
                    const scalar w = 1.0 / std::pow(rp + SMALL, distPow_);
                    num += w * uu[r];
                    den += w;
                }

                scalar u_idw = (den > 0) ? (num / (den + SMALL)) : 0.0;

                // Selection rule:
                // - If affine succeeded AND the stencil is clearly one-sided
                // (gap > 180°),
                //   use affine (true extrapolation).
                // - Otherwise prefer IDW (interior/interpolation).
                const scalar PI = std::acos(scalar(-1));
                if (okAff && gap > PI)
                {
                    values[kNode * N + comp] = u_aff;
                }
                else if (den > 0)
                {
                    values[kNode * N + comp] = u_idw;
                }
                else
                {
                    // ultra-degenerate: 3D Shepard as last resort
                    scalar num3 = 0.0, den3 = 0.0;
                    for (label r = 0; r < k; ++r)
                    {
                        const IndexT idx = idxOrder[r];
                        const scalar d2 = outDistsSqr_[r];
                        const scalar w =
                            1.0 / std::pow(std::sqrt(d2) + SMALL, distPow_);
                        num3 += w * scatterValues[comp * nPoints + idx];
                        den3 += w;
                    }
                    values[kNode * N + comp] =
                        (den3 > 0) ? num3 / (den3 + SMALL) : 0.0;
                }
            }
        }
    }
}
#endif

} // namespace utils
} // namespace accel
