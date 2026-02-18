// File : zoneIO.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "boundary.h"
#include "controls.h"
#include "zone.h"

namespace accel
{

// IO

void zone::read(const YAML::Node& domain)
{
    // zone interior parts
    if (domain["location"])
    {
        std::vector<std::string> location =
            domain["location"].template as<std::vector<std::string>>();
        for (const auto& locationName : location)
        {
            stk::mesh::Part* part =
                meshRef().metaDataRef().get_part(locationName);

            if (!part)
            {
                errorMsg("Part " + locationName + " not found in the mesh");
            }

            interiorParts_.push_back(part);
        }
    }
    else
    {
        errorMsg("zone: `location` block is missing for domain `" + name_ +
                 "`");
    }

    // domain models
    if (domain["domain_models"])
    {
        // is zone moving as rigid body?
        if (domain["domain_models"]["domain_motion"])
        {
            // create domain motion object
            transformationPtr_ = std::make_unique<zoneTransformation>(this);

            // read
            transformationPtr_->read(domain["domain_models"]["domain_motion"]);

            // set some flags
            switch (transformationPtr_->type())
            {
                case meshMotionType::stationary:
                    {
                        // delete the domain motion pointer as we do not need it
                        // throughout
                        transformationPtr_.reset();
                    }
                    break;

                case meshMotionType::rotating:
                    {
                        if (meshPtr_->controlsRef().isTransient())
                        {
                            meshTransforming_ = true;
                        }
                        else
                        {
                            frameRotating_ = true;
                        }
                    }
                    break;

                case meshMotionType::translating:
                    {
                        if (meshPtr_->controlsRef().isTransient())
                        {
                            meshTransforming_ = true;
                        }
                        else
                        {
                            errorMsg(
                                "translation and general mesh motion types are "
                                "only supported for transient simulation");
                        }
                    }
                    break;
            }

            // determine stationary parts of the zone
            if (transformationPtr_ &&
                domain["domain_models"]["domain_motion"]["stationary_parts"])
            {
                std::vector<std::string> stationaryPartNames =
                    domain["domain_models"]["domain_motion"]["stationary_parts"]
                        .template as<std::vector<std::string>>();

                for (auto stationaryPartName : stationaryPartNames)
                {
                    stk::mesh::Part* part =
                        meshRef().metaDataRef().get_part(stationaryPartName);

                    if (!part)
                    {
                        errorMsg("Part " + stationaryPartName +
                                 " not found in the mesh");
                    }

                    stationaryParts_.push_back(part);

                    // check if part is among interior parts
                    if (interiorParts_.end() ==
                        std::find(
                            interiorParts_.begin(), interiorParts_.end(), part))
                    {
                        errorMsg("stationary parts must be parts of the zone");
                    }
                }

                assert(interiorParts_ != stationaryParts_);
            }
        }

        // is zone mesh deforming?
        if (domain["domain_models"]["mesh_deformation"])
        {
            // create mesh deformation object
            deformationPtr_ = std::make_unique<zoneDeformation>(this);

            // read
            deformationPtr_->read(domain["domain_models"]["mesh_deformation"]);

            // set some flags
            if (deformationPtr_->specification() !=
                meshDeformationSpecificationType::none)
            {
                if (meshPtr_->controlsRef().isTransient())
                {
                    meshDeforming_ = true;
                }
                else
                {
                    errorMsg("mesh deformation is only supported for transient "
                             "simulations");
                }
            }
            else
            {
                // delete mesh deformation pointer as we do not need it
                // throughout
                deformationPtr_.reset();
            }

            // determine stationary parts of the zone
            if (deformationPtr_ &&
                domain["domain_models"]["mesh_deformation"]["stationary_parts"])
            {
                std::vector<std::string> stationaryPartNames =
                    domain["domain_models"]["domain_motion"]["stationary_parts"]
                        .template as<std::vector<std::string>>();

                for (auto stationaryPartName : stationaryPartNames)
                {
                    stk::mesh::Part* part =
                        meshRef().metaDataRef().get_part(stationaryPartName);

                    if (!part)
                    {
                        errorMsg("Part " + stationaryPartName +
                                 " not found in the mesh");
                    }

                    stationaryParts_.push_back(part);

                    // check if part is among interior parts
                    if (interiorParts_.end() ==
                        std::find(
                            interiorParts_.begin(), interiorParts_.end(), part))
                    {
                        errorMsg("stationary parts must be parts of the zone");
                    }
                }

                assert(interiorParts_ != stationaryParts_);
            }
        }
    }

    // zone boundary parts
    if (domain["boundaries"])
    {
        const auto& boundariesBlockArray = domain["boundaries"];
        for (label iBoundary = 0;
             iBoundary < static_cast<label>(boundariesBlockArray.size());
             iBoundary++)
        {
            const auto& boundaryBlock = boundariesBlockArray[iBoundary];

            // retrieve boundary name
            std::string name = boundaryBlock["name"].template as<std::string>();

            // get type
            auto type = convertBoundaryPhysicalTypeFromString(
                boundaryBlock["type"].template as<std::string>());

            std::unique_ptr<boundary> boundary_ptr;

            switch (type)
            {
                case boundaryPhysicalType::wall:
                    {
                        // create wall boundary and add to global boundaries
                        // array of the zone
                        boundary_ptr =
                            std::make_unique<wall>(this, name, iBoundary);
                    }
                    break;

                default:
                    {
                        // create boundary and add to global boundaries array of
                        // the zone
                        boundary_ptr =
                            std::make_unique<boundary>(this, name, iBoundary);
                    }
                    break;
            }

            // read
            boundary_ptr->read(boundaryBlock);

            // put into global boundary vector
            boundaryVector_.push_back(std::move(boundary_ptr));
        }
    }
    else
    {
        errorMsg("zone: `boundaries` block is missing for domain `" +
                 domain["name"].template as<std::string>() +
                 "`.\nboundary_conditions:\n -name: <boundary name>\n type: "
                 "<boundary type>\n location: <location list>");
    }
}

} // namespace accel
