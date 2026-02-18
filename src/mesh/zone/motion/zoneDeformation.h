// File : zoneDeformation.h
// Created : Sat Jan 11 2025 20:15:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Mesh deformation specification and displacement diffusion
// parameters
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef ZONEDEFORMATION_H
#define ZONEDEFORMATION_H

// code
#include "vectorUtils.h"

namespace accel
{

struct displacementDiffusionDictionary
{
    meshStiffnessSpecificationType meshStiffnessSpecification_ =
        meshStiffnessSpecificationType::value;

    scalar meshStiffnessValue_ = 1.0;

    scalar modelExponent_ = 2.0;

    scalar referenceLengthScale_ = 1.0;

    // Blended distance and small volumes model parameters
    scalar blendedVolumeWeight_ = 0.5;
    scalar blendedDistanceWeight_ = 0.5;
    scalar blendedVolumeExponent_ = 2.0;
    scalar blendedDistanceExponent_ = 2.0;
};

// forward declaration
class zone;

class zoneDeformation
{
private:
    zone* zonePtr_;

    meshDeformationSpecificationType specification_ =
        meshDeformationSpecificationType::none;

    meshDeformationModel model_ = meshDeformationModel::displacementDiffusion;

    bool displacementRelativeToPreviousMesh_ = true;

    displacementDiffusionDictionary displacementDiffusion_;

public:
    // Constructors
    zoneDeformation(zone* zonePtr);

    // Methods

    void read(const YAML::Node& inputNode);

    void update();

    // Getters and Setters

    meshDeformationSpecificationType specification() const
    {
        return specification_;
    }

    void setSpecification(meshDeformationSpecificationType specification)
    {
        specification_ = specification;
    }

    meshDeformationModel model() const
    {
        return model_;
    }

    void setModel(meshDeformationModel model)
    {
        model_ = model;
    }

    bool displacementRelativeToPreviousMesh() const
    {
        return displacementRelativeToPreviousMesh_;
    }

    void setDisplacementRelativeToPreviousMesh(bool state)
    {
        displacementRelativeToPreviousMesh_ = state;
    }

    displacementDiffusionDictionary& displacementDiffusion()
    {
        return displacementDiffusion_;
    }

    const displacementDiffusionDictionary& displacementDiffusion() const
    {
        return displacementDiffusion_;
    }
};

} // namespace accel

#endif // ZONEDEFORMATION_H
