// File       : zoneTransformation.cpp
// Created    : Sat Jan 11 2025 20:15:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "zoneTransformation.h"
#include "zone.h"

namespace accel
{

// Constructors

zoneTransformation::zoneTransformation(zone* zonePtr) : zonePtr_(zonePtr)
{
}

// Methods

void zoneTransformation::read(const YAML::Node& inputNode)
{
    type_ = convertMeshMotionTypeFromString(
        inputNode["option"].template as<std::string>());

    switch (type_)
    {
        case meshMotionType::stationary:
            break;

        case meshMotionType::translating:
            {
                if (inputNode["velocity"])
                {
                    translation_.v_ =
                        utils::vector(inputNode["velocity"]
                                          .template as<std::vector<scalar>>()
                                          .data());
                }
                else
                {
                    errorMsg("velocity in domain " + zonePtr_->name() +
                             " not provided");
                }
            }
            break;

        case meshMotionType::rotating:
            {
                // read and store origin of rotation
                if (inputNode["origin"])
                {
                    const utils::vector originOfRotation =
                        utils::vector(inputNode["origin"]
                                          .template as<std::vector<scalar>>()
                                          .data());
                    assert(originOfRotation.size() == SPATIAL_DIM);
                    for (int i = 0; i < SPATIAL_DIM; i++)
                    {
                        rotation_.origin_[i] = originOfRotation[i];
                    }
                }
                else
                {
                    errorMsg("origin in domain " + zonePtr_->name() +
                             " not provided");
                }

                // read angular velocity magnitude
                if (inputNode["angular_velocity"])
                {
                    rotation_.omega_ =
                        inputNode["angular_velocity"].template as<scalar>();
                }
                else
                {
                    errorMsg("angular_velocity in domain " + zonePtr_->name() +
                             " not provided");
                }

#if SPATIAL_DIM == 2
                // axis of rotation always in k-dir (any user input will be
                // ignored)
                rotation_.axis_ = {0, 0, 1};

                // fill coriolis matrix
                rotation_.coriolisMatrix_ << 0, -rotation_.omega_,
                    rotation_.omega_, 0;
#elif SPATIAL_DIM == 3
                // read and store axis of rotation
                if (inputNode["axis"])
                {
                    rotation_.axis_ =
                        utils::vector(inputNode["axis"]
                                          .template as<std::vector<scalar>>()
                                          .data());

                    // normalize to ensure unit vector
                    rotation_.axis_.normalize();
                }
                else
                {
                    errorMsg("axis in domain " + zonePtr_->name() +
                             " not provided");
                }

                // fill coriolis matrix
                rotation_.coriolisMatrix_ << 0,
                    -rotation_.omega_ * rotation_.axis_[2],
                    rotation_.omega_ * rotation_.axis_[1],
                    rotation_.omega_ * rotation_.axis_[2], 0,
                    -rotation_.omega_ * rotation_.axis_[0],
                    -rotation_.omega_ * rotation_.axis_[1],
                    rotation_.omega_ * rotation_.axis_[0], 0;
#endif
            }
            break;
    }
}

void zoneTransformation::update()
{
}

} // namespace accel
