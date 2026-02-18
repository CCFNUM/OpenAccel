// File       : interfaceIO.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

#include "interface.h"
#include "messager.h"
#include "zone.h"

namespace accel
{

// IO

void interface::read(const YAML::Node& inputNode)
{
    // Primary interface attributes: obligatory
    option_ = convertInterfaceModelOptionFromString(
        inputNode["option"].template as<std::string>());

    // In case of rotational periodicity, rotation axis and axis location
    // are obligatory
    if (option_ == interfaceModelOption::rotationalPeriodicity)
    {
        // read rotation axis and check dims
        rotationAxis_ = utils::vector(inputNode["rotation_axis"]
                                          .template as<std::vector<scalar>>()
                                          .data());

        // read axis location and check dims
        axisLocation_ = utils::vector(inputNode["axis_location"]
                                          .template as<std::vector<scalar>>()
                                          .data());
    }

    // fluid-fluid? fluid-solid? solid-solid? a further check of consistency
    // will be made later
    if (inputNode["type"])
    {
        type_ = convertInterfaceTypeFromString(
            inputNode["type"].template as<std::string>());
    }

    // non-overlap will be treated as a slip wall?
    if (inputNode["slip_non_overlap"])
    {
        isSlipNonOverlap_ = inputNode["slip_non_overlap"].template as<bool>();
    }

    // conformality check tolerance
    if (inputNode["conformality_check_tolerance"])
    {
        conformalityTolerance_ =
            inputNode["conformality_check_tolerance"].template as<scalar>();
    }

    // force non-conformal treatment for conformal interface: this enables DG
    // even if conformal, but it is always true for a fluid-solid interface
    if (type_ == interfaceType::fluid_solid)
    {
        isForceNonconformalTreatment_ = true;

        if (inputNode["force_nonconformal_treatment"])
        {
            if (messager::master())
            {
                warningMsg("force_nonconformal_treatment key will be ignored "
                           "for fluid-solid interface. DG is ALWAYS enabled "
                           "for fluid-solid interfaces");
            }
        }
    }
    else
    {
        if (inputNode["force_nonconformal_treatment"])
        {
            isForceNonconformalTreatment_ =
                inputNode["force_nonconformal_treatment"].template as<bool>();
        }
    }

    // overlap check tolerance
    if (inputNode["overlap_check_tolerance"])
    {
        overlapTolerance_ =
            inputNode["overlap_check_tolerance"].template as<scalar>();
    }

    // penalty factor is by default 1, a smaller value may produce smoother
    // coefficients
    if (inputNode["penalty_factor"])
    {
        penaltyFactor_ = inputNode["penalty_factor"].template as<scalar>();
    }

    std::string masterZoneName, slaveZoneName;
    std::vector<std::string> masterRegionNameList, slaveRegionNameList;

    // Master side
    if (inputNode["side1"])
    {
        masterZoneName =
            inputNode["side1"]["domain"].template as<std::string>();
        masterRegionNameList = inputNode["side1"]["region_list"]
                                   .template as<std::vector<std::string>>();
    }
    else
    {
        errorMsg("side1 block is not provided in the yaml input file");
    }

    // Slave side
    if (inputNode["side2"])
    {
        slaveZoneName = inputNode["side2"]["domain"].template as<std::string>();
        slaveRegionNameList = inputNode["side2"]["region_list"]
                                  .template as<std::vector<std::string>>();
    }
    else
    {
        errorMsg("side2 block is not provided in the yaml input file");
    }

    // find pair zone names (the zone to which each side belong)
    pairZoneNames_ = std::make_pair(masterZoneName, slaveZoneName);

    // determine the indices of the two zones
    label masterZoneIndex = -1;
    label slaveZoneIndex = -1;
    for (zone* zone : this->meshPtr()->zoneVector())
    {
        if (zone->name() == masterZoneName)
        {
            masterZoneIndex = zone->index();
        }

        if (zone->name() == slaveZoneName)
        {
            slaveZoneIndex = zone->index();
        }
    }

    assert(masterZoneIndex != -1);
    assert(slaveZoneIndex != -1);

    pairZoneIndices_ = std::make_pair(masterZoneIndex, slaveZoneIndex);

    // Parameters for interface sides

    // search-related parameters
    scalar expandBoxPercentage = 0.0;
    scalar searchTolerance = 0.01;
    bool activateDynamicSearchAlgorithm = false;
    bool useShifted = false;
    bool clipIsoParametricCoords = false;
    std::string searchMethodName = "stk_kdtree";

    if (inputNode["expand_box_percentage"])
    {
        expandBoxPercentage =
            inputNode["expand_box_percentage"].template as<scalar>();
    }

    if (inputNode["clip_isoparametric_coordinates"])
    {
        clipIsoParametricCoords =
            inputNode["clip_isoparametric_coordinates"].template as<bool>();
    }

    if (inputNode["search_method"])
    {
        searchMethodName =
            inputNode["search_method"].template as<std::string>();
    }

    if (inputNode["search_tolerance"])
    {
        searchTolerance = inputNode["search_tolerance"].template as<scalar>();
    }

    if (inputNode["activate_dynamic_search_algorithm"])
    {
        activateDynamicSearchAlgorithm =
            inputNode["activate_dynamic_search_algorithm"].template as<bool>();
    }

    if (inputNode["gauss_lobatto_quadrature"])
    {
        useShifted = inputNode["gauss_lobatto_quadrature"].template as<bool>();
    }

    // Find pair parts
    stk::mesh::PartVector masterParts, slaveParts;
    for (const auto& locationName : masterRegionNameList)
    {
        masterParts.push_back(meshRef().metaDataRef().get_part(locationName));
    }
    for (const auto& locationName : slaveRegionNameList)
    {
        slaveParts.push_back(meshRef().metaDataRef().get_part(locationName));
    }

    // Master side
    masterInfoPtr_ =
        std::make_unique<interfaceSideInfo>(this,
                                            true,
                                            masterParts,
                                            slaveParts,
                                            option_,
                                            expandBoxPercentage / 100.0,
                                            searchMethodName,
                                            clipIsoParametricCoords,
                                            searchTolerance,
                                            activateDynamicSearchAlgorithm,
                                            useShifted,
                                            name_ + "_master_side");

    // Slave side
    slaveInfoPtr_ =
        std::make_unique<interfaceSideInfo>(this,
                                            false,
                                            slaveParts,
                                            masterParts,
                                            option_,
                                            expandBoxPercentage / 100.0,
                                            searchMethodName,
                                            clipIsoParametricCoords,
                                            searchTolerance,
                                            activateDynamicSearchAlgorithm,
                                            useShifted,
                                            name_ + "_slave_side");
}

} // namespace accel

#endif /* HAS_INTERFACE */
