// File : zoneDeformation.cpp
// Created : Sat Jan 11 2025 20:15:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "zoneDeformation.h"
#include "zone.h"

namespace accel
{

// Constructors

zoneDeformation::zoneDeformation(zone* zonePtr) : zonePtr_(zonePtr)
{
}

// Methods

void zoneDeformation::read(const YAML::Node& inputNode)
{
    if (inputNode["option"])
    {
        specification_ = convertMeshDeformationSpecificationFromString(
            inputNode["option"].template as<std::string>());

        switch (specification_)
        {
            case meshDeformationSpecificationType::none:
                break;

            case meshDeformationSpecificationType::regionsOfMotionSpecified:
                {
                    if (inputNode["mesh_motion_model"])
                    {
                        const auto& meshMotionModelBlock =
                            inputNode["mesh_motion_model"];
                        model_ = convertMeshDeformationModelTypeFromString(
                            meshMotionModelBlock["option"]
                                .template as<std::string>());

                        if (meshMotionModelBlock["mesh_stiffness"])
                        {
                            const auto& meshStiffnessBlock =
                                meshMotionModelBlock["mesh_stiffness"];

                            if (meshStiffnessBlock["option"])
                            {
                                displacementDiffusion_
                                    .meshStiffnessSpecification_ =
                                    convertMeshStiffnessSpecificationTypeFromString(
                                        meshStiffnessBlock["option"]
                                            .template as<std::string>());

                                if (displacementDiffusion_
                                        .meshStiffnessSpecification_ ==
                                    meshStiffnessSpecificationType::value)
                                {
                                    if (meshStiffnessBlock["value"])
                                    {
                                        displacementDiffusion_
                                            .meshStiffnessValue_ =
                                            meshStiffnessBlock["value"]
                                                .template as<scalar>();
                                    }
                                    else
                                    {
                                        errorMsg("Mesh stiffness value not "
                                                 "provided");
                                    }
                                }
                                else if (displacementDiffusion_
                                             .meshStiffnessSpecification_ ==
                                         meshStiffnessSpecificationType::
                                             increaseNearSmallVolumes)
                                {
                                    if (meshStiffnessBlock["model_exponent"])
                                    {
                                        displacementDiffusion_.modelExponent_ =
                                            meshStiffnessBlock["model_exponent"]
                                                .template as<scalar>();
                                    }
                                    else
                                    {
                                        errorMsg(
                                            "Mesh stiffness model exponent not "
                                            "provided");
                                    }
                                }
                                else if (displacementDiffusion_
                                             .meshStiffnessSpecification_ ==
                                         meshStiffnessSpecificationType::
                                             increaseNearBoundaries)
                                {
                                    if (meshStiffnessBlock["model_exponent"])
                                    {
                                        displacementDiffusion_.modelExponent_ =
                                            meshStiffnessBlock["model_exponent"]
                                                .template as<scalar>();
                                    }
                                    else
                                    {
                                        errorMsg(
                                            "Mesh stiffness model exponent not "
                                            "provided");
                                    }

                                    if (meshStiffnessBlock
                                            ["reference_length_scale"])
                                    {
                                        displacementDiffusion_
                                            .referenceLengthScale_ =
                                            meshStiffnessBlock
                                                ["reference_length_scale"]
                                                    .template as<scalar>();
                                    }
                                }
                                else if (displacementDiffusion_
                                             .meshStiffnessSpecification_ ==
                                         meshStiffnessSpecificationType::
                                             blendedDistanceAndSmallVolumes)
                                {
                                    if (meshStiffnessBlock["volume_weight"])
                                    {
                                        displacementDiffusion_
                                            .blendedVolumeWeight_ =
                                            meshStiffnessBlock["volume_weight"]
                                                .template as<scalar>();
                                    }

                                    if (meshStiffnessBlock["distance_weight"])
                                    {
                                        displacementDiffusion_
                                            .blendedDistanceWeight_ =
                                            meshStiffnessBlock
                                                ["distance_weight"]
                                                    .template as<scalar>();
                                    }

                                    if (meshStiffnessBlock["volume_exponent"])
                                    {
                                        displacementDiffusion_
                                            .blendedVolumeExponent_ =
                                            meshStiffnessBlock
                                                ["volume_exponent"]
                                                    .template as<scalar>();
                                    }

                                    if (meshStiffnessBlock["distance_exponent"])
                                    {
                                        displacementDiffusion_
                                            .blendedDistanceExponent_ =
                                            meshStiffnessBlock
                                                ["distance_exponent"]
                                                    .template as<scalar>();
                                    }
                                }
                            }
                            else
                            {
                                errorMsg(
                                    "Option for mesh stiffness not provided");
                            }
                        }
                        else
                        {
                            errorMsg("Mesh stiffness not provided");
                        }
                    }
                    else
                    {
                        errorMsg("Mesh motion model not provided");
                    }

                    if (inputNode["displacement_relative_to"])
                    {
                        if (inputNode["displacement_relative_to"]
                                .template as<std::string>() == "previous_mesh")
                        {
                            displacementRelativeToPreviousMesh_ = true;
                        }
                        else if (inputNode["displacement_relative_to"]
                                     .template as<std::string>() ==
                                 "initial_mesh")
                        {
                            displacementRelativeToPreviousMesh_ = false;
                        }
                        else
                        {
                            errorMsg("Unrecognized value " +
                                     inputNode["displacement_relative_to"]
                                         .template as<std::string>() +
                                     " for displacement_relative_to in the "
                                     "mesh_deformation block");
                        }
                    }
                }
                break;

            case meshDeformationSpecificationType::inherent:
                break;
        }
    }
    else
    {
        errorMsg("option key not provided for mesh deformation");
    }
}

void zoneDeformation::update()
{
}

} // namespace accel
