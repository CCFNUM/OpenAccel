// File : rigidBody.cpp
// Created : Thu Feb 13 2025 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "rigidBody.h"
#include "realm.h"

namespace accel
{

rigidBody::rigidBody(realm* realmPtr) : realmPtr_(realmPtr)
{
}

void rigidBody::read(const YAML::Node& input)
{
    // set name
    name_ = input["name"].template as<std::string>();

    // collect parts
    if (input["basic_settings"])
    {
        auto names = input["basic_settings"]["location"]
                         .template as<std::vector<std::string>>();
        for (auto name : names)
        {
            // get part
            const auto& part =
                realmPtr_->meshPtr()->metaDataPtr()->get_part(name);
            parts_.push_back(part);
        }
    }
    else
    {
        errorMsg("location block is not provided in rigid body block basic "
                 "settings");
    }

    // basic settings
    if (input["basic_settings"])
    {
        const YAML::Node& basicSettingsBlock = input["basic_settings"];

        // set mass
        m_ = basicSettingsBlock["mass"].template as<scalar>();

        // set moment of inertia
        I_(0, 0) = basicSettingsBlock["mass_moment_of_inertia"]["xx"]
                       .template as<scalar>();
        I_(1, 1) = basicSettingsBlock["mass_moment_of_inertia"]["yy"]
                       .template as<scalar>();
        I_(2, 2) = basicSettingsBlock["mass_moment_of_inertia"]["zz"]
                       .template as<scalar>();
        I_(0, 1) = basicSettingsBlock["mass_moment_of_inertia"]["xy"]
                       .template as<scalar>();
        I_(0, 2) = basicSettingsBlock["mass_moment_of_inertia"]["xz"]
                       .template as<scalar>();
        I_(1, 2) = basicSettingsBlock["mass_moment_of_inertia"]["yz"]
                       .template as<scalar>();
        I_(1, 0) = I_(0, 1);
        I_(2, 0) = I_(0, 2);
        I_(2, 1) = I_(1, 2);
    }
    else
    {
        errorMsg("basic settngs block is not provided in rigid body block");
    }

    // dynamics
    if (input["dynamics"])
    {
        const YAML::Node& dynamicsBlock = input["dynamics"];

        if (dynamicsBlock["external_force_definitions"])
        {
            for (const auto& FextDef :
                 dynamicsBlock["external_force_definitions"])
            {
                Fext_ += utils::vector(
                    FextDef["value"].template as<std::vector<scalar>>().data());
            }
        }

        if (dynamicsBlock["external_moment_definitions"])
        {
            for (const auto& MextDef :
                 dynamicsBlock["external_moment_definitions"])
            {
                Mext_ += utils::vector(
                    MextDef["value"].template as<std::vector<scalar>>().data());
            }
        }

        if (dynamicsBlock["degrees_of_freedom"])
        {
            const auto& DOFBlock = dynamicsBlock["degrees_of_freedom"];

            if (DOFBlock["translational_degrees_of_freedom"])
            {
                const auto& transDOFBlock =
                    DOFBlock["translational_degrees_of_freedom"];
                std::string option =
                    transDOFBlock["option"].template as<std::string>();
                if (option == "none")
                {
                }
                else if (option == "x_axis")
                {
                    linearDOFs_(0) = 1;
                }
                else if (option == "y_axis")
                {
                    linearDOFs_(1) = 1;
                }
                else if (option == "z_axis")
                {
                    linearDOFs_(2) = 1;
                }
                else if (option == "x_and_y_axes")
                {
                    linearDOFs_(0) = 1;
                    linearDOFs_(1) = 1;
                }
                else if (option == "y_and_z_axes")
                {
                    linearDOFs_(1) = 1;
                    linearDOFs_(2) = 1;
                }
                else if (option == "x_and_z_axes")
                {
                    linearDOFs_(0) = 1;
                    linearDOFs_(2) = 1;
                }
                else if (option == "x_y_and_z_axes")
                {
                    linearDOFs_(0) = 1;
                    linearDOFs_(1) = 1;
                    linearDOFs_(2) = 1;
                }
            }

            if (DOFBlock["rotational_degrees_of_freedom"])
            {
                const auto& rotDOFBlock =
                    DOFBlock["rotational_degrees_of_freedom"];
                std::string option =
                    rotDOFBlock["option"].template as<std::string>();
                if (option == "none")
                {
                }
                else if (option == "x_axis")
                {
                    angularDOFs_(0) = 1;
                }
                else if (option == "y_axis")
                {
                    angularDOFs_(1) = 1;
                }
                else if (option == "z_axis")
                {
                    angularDOFs_(2) = 1;
                }
                else if (option == "x_and_y_axes")
                {
                    angularDOFs_(0) = 1;
                    angularDOFs_(1) = 1;
                }
                else if (option == "y_and_z_axes")
                {
                    angularDOFs_(1) = 1;
                    angularDOFs_(2) = 1;
                }
                else if (option == "x_and_z_axes")
                {
                    angularDOFs_(0) = 1;
                    angularDOFs_(2) = 1;
                }
                else if (option == "x_y_and_z_axes")
                {
                    angularDOFs_(0) = 1;
                    angularDOFs_(1) = 1;
                    angularDOFs_(2) = 1;
                }
            }
        }
    }
    else
    {
        errorMsg("dynamics block is not provided in rigid body block");
    }

    // initial conditions
    if (input["initial_conditions"])
    {
        const YAML::Node& initialConditionsBlock = input["initial_conditions"];

        if (initialConditionsBlock["linear_velocity"])
        {
            U_ = utils::vector(initialConditionsBlock["linear_velocity"]
                                   .template as<std::vector<scalar>>()
                                   .data());
        }

        if (initialConditionsBlock["linear_acceleration"])
        {
            a_ = utils::vector(initialConditionsBlock["linear_acceleration"]
                                   .template as<std::vector<scalar>>()
                                   .data());
        }

        if (initialConditionsBlock["angular_velocity"])
        {
            omega_ = utils::vector(initialConditionsBlock["angular_velocity"]
                                       .template as<std::vector<scalar>>()
                                       .data());
        }

        if (initialConditionsBlock["angular_acceleration"])
        {
            alpha_ =
                utils::vector(initialConditionsBlock["angular_acceleration"]
                                  .template as<std::vector<scalar>>()
                                  .data());
        }
    }
    else
    {
        errorMsg(
            "initial conditions block is not provided in rigid body block");
    }
}

void rigidBody::update(const utils::vector F, const utils::vector M, scalar dt)
{
    // Get acceleration and update velocity from updated force
    if (m_ > 0.0)
    {
        a_ = a_ * 0.5 + 0.5 * (F + Fext_) / m_;
        U_ = U_ + dt * a_;
    }

    // enable/disable dof's
    U_.array() *= linearDOFs_.array().cast<scalar>();

    // Convert back into body frame
    utils::vector bodyFrameOme = convertVectToOrigFrame(omega_, theta_, false);
    utils::vector bodyFrameMom =
        convertVectToOrigFrame(M + Mext_, theta_, false);

    // advance omega
    {
        // Advance using RK4
        utils::vector omegam(bodyFrameOme);
        std::vector<scalar> rk_fact(4, 0.0);
        scalar k[3][4];

        rk_fact[0] = 0.5 * dt;
        rk_fact[1] = 0.5 * dt;
        rk_fact[2] = dt;

        for (size_t i = 0; i < 4; ++i)
        {
            k[0][i] = 0.0;
            k[1][i] = 0.0;
            k[2][i] = 0.0;
        }

        for (size_t i = 0; i < 4; ++i)
        {
            k[0][i] = assessRigidEulerRHS(bodyFrameMom(0),
                                          omegam(2),
                                          omegam(1),
                                          I_(0, 0),
                                          I_(2, 2),
                                          I_(1, 1));
            k[1][i] = assessRigidEulerRHS(bodyFrameMom(1),
                                          omegam(0),
                                          omegam(2),
                                          I_(1, 1),
                                          I_(0, 0),
                                          I_(2, 2));
            k[2][i] = assessRigidEulerRHS(bodyFrameMom(2),
                                          omegam(1),
                                          omegam(0),
                                          I_(2, 2),
                                          I_(1, 1),
                                          I_(0, 0));
            for (size_t j = 0; j < 3; ++j)
            {
                omegam[j] = bodyFrameOme[j] + rk_fact[i] * dt;
            }
        }

        bodyFrameOme[0] +=
            (dt / 6.0) * (k[0][0] + 2.0 * k[0][1] + 2.0 * k[0][2] + k[0][3]);
        bodyFrameOme[1] +=
            (dt / 6.0) * (k[1][0] + 2.0 * k[1][1] + 2.0 * k[1][2] + k[1][3]);
        bodyFrameOme[2] +=
            (dt / 6.0) * (k[2][0] + 2.0 * k[2][1] + 2.0 * k[2][2] + k[2][3]);
    }

    // Convert back to lab frame
    omega_ = convertVectToOrigFrame(bodyFrameOme, theta_, true);

    // enable/disable dof's
    omega_.array() *= angularDOFs_.array().cast<scalar>();

    // update body centroid displacement: total displacement
    dx_ += U_ * dt;

    // angle -> quaternion -> RK4 using Omega -> angle
    std::vector<scalar> quat = getQuatFromEulerxyz(theta_);
    std::vector<scalar> quatm(quat);
    std::vector<scalar> qrhs(4, 0.0);
    std::vector<scalar> rk_fact(4, 0.0);
    scalar k[4][4];

    rk_fact[0] = 0.5 * dt;
    rk_fact[1] = 0.5 * dt;
    rk_fact[2] = dt;

    for (size_t i = 0; i < 4; ++i)
    {
        evaluateQuatRHS(quatm, omega_, qrhs);
        std::copy(qrhs.begin(), qrhs.end(), k[i]);
        for (size_t j = 0; j < 4; ++j)
        {
            quatm[j] = quat[j] + dt * rk_fact[i] * qrhs[j];
        }
    }
    quat[0] += dt / 6.0 * (k[0][0] + 2.0 * k[1][0] + 2.0 * k[2][0] + k[3][0]);
    quat[1] += dt / 6.0 * (k[0][1] + 2.0 * k[1][1] + 2.0 * k[2][1] + k[3][1]);
    quat[2] += dt / 6.0 * (k[0][2] + 2.0 * k[1][2] + 2.0 * k[2][2] + k[3][2]);
    quat[3] += dt / 6.0 * (k[0][3] + 2.0 * k[1][3] + 2.0 * k[2][3] + k[3][3]);

    // Extract angle from Quat
    theta_ = rvector3(getEulerxyzFromQuat(quat).data());
}

std::vector<scalar>
rigidBody::getQuatFromEulerxyz(const utils::vector& bodyAngle)
{
    std::vector<scalar> quat(4, 0.0);
    std::vector<scalar> c(3);
    std::vector<scalar> s(3);

    c[0] = std::cos(bodyAngle(0) * 0.5);
    c[1] = std::cos(bodyAngle(1) * 0.5);
    c[2] = std::cos(bodyAngle(2) * 0.5);
    s[0] = std::sin(bodyAngle(0) * 0.5);
    s[1] = std::sin(bodyAngle(1) * 0.5);
    s[2] = std::sin(bodyAngle(2) * 0.5);

    quat[0] = c[0] * c[1] * c[2] + s[0] * s[1] * s[2];
    quat[1] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
    quat[2] = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
    quat[3] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];

    scalar qmag = std::sqrt(quat[0] * quat[0] + quat[1] * quat[1] +
                            quat[2] * quat[2] + quat[3] * quat[3]);

    for (size_t i = 0; i < 4; ++i)
        quat[i] /= qmag;

    return quat;
}

std::vector<scalar>
rigidBody::getEulerxyzFromQuat(const std::vector<scalar>& quat)
{
    std::vector<scalar> angle(3, 0.0);
    angle[0] = std::atan2(2.0 * (quat[0] * quat[1] + quat[2] * quat[3]),
                          quat[0] * quat[0] - quat[1] * quat[1] -
                              quat[2] * quat[2] + quat[3] * quat[3]);
    angle[1] = std::asin(-2.0 * (quat[1] * quat[3] - quat[0] * quat[2]));
    angle[2] = std::atan2(2.0 * (quat[2] * quat[1] + quat[0] * quat[3]),
                          quat[0] * quat[0] + quat[1] * quat[1] -
                              quat[2] * quat[2] - quat[3] * quat[3]);

    return angle;
}

utils::matrix rigidBody::getRotmatFromQuat(std::vector<scalar>& quat,
                                           bool to_frame)
{
    utils::matrix rotMat;

    // assuming rotMat was originally row ordered
    rotMat(0, 0) = 2.0 * (quat[0] * quat[0] + quat[1] * quat[1] - 0.5);
    rotMat(0, 1) = 2.0 * (quat[1] * quat[2] - quat[0] * quat[3]);
    rotMat(0, 2) = 2.0 * (quat[1] * quat[3] + quat[0] * quat[2]);
    rotMat(1, 0) = 2.0 * (quat[1] * quat[2] + quat[0] * quat[3]);
    rotMat(1, 1) = 2.0 * (quat[0] * quat[0] + quat[2] * quat[2] - 0.5);
    rotMat(1, 2) = 2.0 * (quat[2] * quat[3] - quat[0] * quat[1]);
    rotMat(2, 0) = 2.0 * (quat[1] * quat[3] - quat[0] * quat[2]);
    rotMat(2, 1) = 2.0 * (quat[2] * quat[3] + quat[0] * quat[1]);
    rotMat(2, 2) = 2.0 * (quat[0] * quat[0] + quat[3] * quat[3] - 0.5);

    if (!to_frame)
    {
        // rotation matrices are always orthogonal
        return rotMat.transpose();

        // CODE BELOW WON'T WORK ANYMORE
        // // Calculate inverse of rotation matrix
        // scalar det =
        // p_rotMat[0] *
        // (p_rotMat[4] * p_rotMat[8] - p_rotMat[7] * p_rotMat[5]) -
        // p_rotMat[1] *
        // (p_rotMat[3] * p_rotMat[8] - p_rotMat[6] * p_rotMat[5]) +
        // p_rotMat[2] *
        // (p_rotMat[3] * p_rotMat[7] - p_rotMat[4] * p_rotMat[6]);

        // if (det < FLT_MIN)
        // {
        // rotMat.fill(0);
        // p_rotMat[0] = 1.0;
        // p_rotMat[4] = 1.0;
        // p_rotMat[8] = 1.0;
        // }
        // else
        // {
        // utils::matrix rotMati;
        // scalar* p_rotMati = rotMati.data();

        // det = 1.0 / det;
        // p_rotMati[0] =
        // (p_rotMat[4] * p_rotMat[8] - p_rotMat[5] * p_rotMat[7]) *
        // det;
        // p_rotMati[1] =
        // -(p_rotMat[1] * p_rotMat[8] - p_rotMat[2] * p_rotMat[7]) *
        // det;
        // p_rotMati[2] =
        // (p_rotMat[1] * p_rotMat[5] - p_rotMat[2] * p_rotMat[4]) *
        // det;
        // p_rotMati[3] =
        // -(p_rotMat[3] * p_rotMat[8] - p_rotMat[5] * p_rotMat[6]) *
        // det;
        // p_rotMati[4] =
        // (p_rotMat[0] * p_rotMat[8] - p_rotMat[2] * p_rotMat[6]) *
        // det;
        // p_rotMati[5] =
        // -(p_rotMat[0] * p_rotMat[5] - p_rotMat[3] * p_rotMat[2]) *
        // det;
        // p_rotMati[6] =
        // (p_rotMat[3] * p_rotMat[7] - p_rotMat[4] * p_rotMat[6]) *
        // det;
        // p_rotMati[7] =
        // -(p_rotMat[0] * p_rotMat[7] - p_rotMat[1] * p_rotMat[6]) *
        // det;
        // p_rotMati[8] =
        // (p_rotMat[0] * p_rotMat[4] - p_rotMat[1] * p_rotMat[3]) *
        // det;

        // return rotMati;
        // }
    }
    return rotMat;
}

utils::vector rigidBody::convertVectToOrigFrame(const utils::vector& vect,
                                                const utils::vector& bodyAngle,
                                                bool toLabFrame)
{
    // Get quaternion from angles
    std::vector<scalar> quat = getQuatFromEulerxyz(bodyAngle);

    // Get rotation matrix
    utils::matrix rotMat = getRotmatFromQuat(quat, toLabFrame);

    utils::vector vectReturn = rotMat * vect;

    return vectReturn;
}

scalar rigidBody::assessRigidEulerRHS(scalar momentRHS,
                                      scalar omega1RHS,
                                      scalar omega2RHS,
                                      scalar i0LHS,
                                      scalar i1RHS,
                                      scalar i2RHS)
{
    scalar rhs = 0.0;
    if (i0LHS > FLT_MIN)
    {
        rhs = (momentRHS - (i1RHS - i2RHS) * omega1RHS * omega2RHS) / i0LHS;
    }

    return rhs;
}

void rigidBody::evaluateQuatRHS(std::vector<scalar>& quat,
                                utils::vector& omega,
                                std::vector<scalar>& qrhs)
{
    qrhs[0] =
        (-quat[1] * omega(0) - quat[2] * omega(1) - quat[3] * omega(2)) * 0.5;
    qrhs[1] =
        (quat[0] * omega(0) + quat[3] * omega(1) - quat[2] * omega(2)) * 0.5;
    qrhs[2] =
        (-quat[3] * omega(0) + quat[0] * omega(1) + quat[1] * omega(2)) * 0.5;
    qrhs[3] =
        (quat[2] * omega(0) - quat[1] * omega(1) + quat[0] * omega(2)) * 0.5;
}

} // namespace accel
