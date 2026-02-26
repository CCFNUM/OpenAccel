// File       : vectorUtils.h
// Created    : Sat Dec 21 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Vector and tensor transformation utilities with rotation matrix
// support
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef VECTORUTILS_H
#define VECTORUTILS_H

// code
#include "types.h"

namespace accel
{

namespace utils
{
using vector = rvectorD;
using matrix = rtensorD;
using vectorView = rvectorDView;
using vectorViewC = rvectorDViewC;
using tensorView = rtensorDView;
using tensorViewC = rtensorDViewC;

inline vector transformVector(const matrix& Q, const vector& v)
{
    return Q * v;
}

inline matrix transformTensor(const matrix& Q, const matrix& tensor)
{
    return Q * tensor * Q.transpose();
}

inline void transformVectorInPlace(const matrix& Q, Eigen::Ref<vector> v)
{
    v = Q * v;
}

inline void transformTensorInPlace(const matrix& Q, Eigen::Ref<matrix> tensor)
{
    tensor = Q * tensor * Q.transpose();
}
#if SPATIAL_DIM == 2
inline rtensor2 getRotationMatrix(scalar angle, const vector& axis)
{
    return Eigen::Rotation2D<scalar>(angle).toRotationMatrix();
}
#else
inline rtensor3 getRotationMatrix(scalar angle, const vector& axis)
{
    return matrix(Eigen::AngleAxis<scalar>(angle, axis));
}
#endif

// Orthonormal basis around axis
class basis
{
private:
    vector a_;  // unit axis
    vector e1_; // unit, perpendicular to a
    vector e2_; // e2 = a x e1

public:
    basis() = default;

    basis(const vector& axis)
    {
        a_ = axis.normalized();
        vector tmp =
            (std::fabs(a_[0]) < 0.9) ? vector(1, 0, 0) : vector(0, 1, 0);
        e1_ = (a_.cross(tmp)).normalized();
        e2_ = a_.cross(e1_);
    }

    inline void projectToPlane(const vector& x3,
                               const vector& p,
                               scalar& X,
                               scalar& Y) const
    {
        // remove axis component relative to point on axis (p)
        vector r = x3 - p;
        vector rperp = r - a_ * a_.dot(r);
        X = rperp.dot(e1_);
        Y = rperp.dot(e2_);
    }
};

} // namespace utils

} // namespace accel

#endif // VECTORUTILS_H
