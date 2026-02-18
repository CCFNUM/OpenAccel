// File : zoneTransformation.h
// Created : Sat Jan 11 2025 20:15:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Rigid body translation and rotation motion for mesh zones
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ZONETRANSFORMATION_H
#define ZONETRANSFORMATION_H

// code
#include "vectorUtils.h"

namespace accel
{

struct translationDictionary
{
    utils::vector v_;
};

struct rotationDictionary
{
    scalar omega_;

    utils::vector origin_;

    utils::vector axis_;

    utils::matrix coriolisMatrix_;
};

// forward declaration
class zone;

class zoneTransformation
{
private:
    zone* zonePtr_;

    meshMotionType type_ = meshMotionType::stationary;

    translationDictionary translation_;

    rotationDictionary rotation_;

public:
    // Constructors
    zoneTransformation(zone* zonePtr);

    // Methods

    void read(const YAML::Node& inputNode);

    void update();

    // Access

    meshMotionType type() const
    {
        return type_;
    }

    translationDictionary& translation()
    {
        return translation_;
    }

    const translationDictionary& translation() const
    {
        return translation_;
    }

    rotationDictionary& rotation()
    {
        return rotation_;
    }

    const rotationDictionary& rotation() const
    {
        return rotation_;
    }
};

} // namespace accel

#endif // ZONETRANSFORMATION_H
