// File       : rigidBody.h
// Created    : Thu Feb 13 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Rigid body dynamics with linear and angular motion
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RIGIDBODY_H
#define RIGIDBODY_H

// code
#include "mesh.h"
#include "vectorUtils.h"

namespace accel
{

class realm;

class rigidBody
{
private:
    realm* realmPtr_;

    stk::mesh::PartVector parts_;

    // ID for the rigid body
    std::string name_;

    // mass
    scalar m_ = 0.0;

    // mass momentum of intertia
    utils::matrix I_;

    // linear acceleration
    utils::vector a_;

    // linear velocity
    utils::vector U_;

    // angular acceleration
    utils::vector alpha_;

    // angular velocity
    utils::vector omega_;

    // External force applied on body
    utils::vector Fext_;

    // External torque applied on body
    utils::vector Mext_;

    // angle measure from previous position
    utils::vector theta_;

    // centroid displacement
    utils::vector dx_;

    // enabled directions in linear motion
    ivector3 linearDOFs_ = ivector3::Zero();

    // enabled directions in angular motion
    ivector3 angularDOFs_ = ivector3::Zero();

public:
    // Constructors
    rigidBody(realm* realmPtr);

    // IO
    void read(const YAML::Node& input);

    // Methods

    void update(const utils::vector F, const utils::vector M, scalar dt);

    // Access

    std::string name() const
    {
        return name_;
    }

    utils::vector U() const
    {
        return U_;
    }

    utils::vector omega() const
    {
        return omega_;
    }

    utils::vector theta() const
    {
        return theta_;
    }

    utils::vector dx() const
    {
        return dx_;
    }

    // Static functions

    static std::vector<scalar>
    getQuatFromEulerxyz(const utils::vector& bodyAngle);

    static std::vector<scalar>
    getEulerxyzFromQuat(const std::vector<scalar>& quat);

    static utils::matrix getRotmatFromQuat(std::vector<scalar>& quat,
                                           bool to_frame = true);

    static utils::vector convertVectToOrigFrame(const utils::vector& vect,
                                                const utils::vector& bodyAngle,
                                                bool toLabFrame);

    static scalar assessRigidEulerRHS(scalar momentRHS,
                                      scalar omega1RHS,
                                      scalar omega2RHS,
                                      scalar i0LHS,
                                      scalar i1RHS,
                                      scalar i2RHS);

    static void evaluateQuatRHS(std::vector<scalar>& quat,
                                utils::vector& omega,
                                std::vector<scalar>& qrhs);
};

} // namespace accel

#endif // RIGIDBODY_H
