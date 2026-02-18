// File : fieldBroker.cpp
// Created : Tue Feb 20 2024 12:55:24 (+0100)
// Author : Fabian Wermelinger
// Description: Field broker implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "fieldBroker.h"
#include "dataHandler.h"
#include "domain.h"
#include "initialConditions.h"
#include "macros.h"
#include "realm.h"
#include "types.h"

namespace accel
{

fieldBroker::fieldBroker(realm* realm) : realmPtr_(realm)
{
    // Create transport field instances (see also TODO for protected transport
    // field API in fieldBroker.h)
    URef();
    rhoRef();
    mDotRef();

    // if any frame of a domain is moving, a new field "relative velocity"
    // will be calculated as a post-processed field
    if (this->meshRef().anyZoneFrameRotating() ||
        this->meshRef().anyZoneMeshTransforming() ||
        this->meshRef().anyZoneMeshDeforming())
    {
        UrRef();
    }
}

void fieldBroker::setupVelocity(const std::shared_ptr<domain> domain)
{
    if (URef().isZoneUnset(domain->index()))
    {
        URef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            URef(), realm::U_ID, domain);

#ifdef HAS_INTERFACE
        // 1) Register velocity side fields on fluid-solid interface
        // (only fluid)
        //    side, because we deal with this as a typical no-slip wall
        // 2) Register also on an interface side which has a non-overlap
        // portion
        //    because the non-overlap portion might be treated as as a no-slip
        //    subject to user's preference
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                if (interf->isSlipNonOverlap())
                    continue;

                // register side fields only if there are non-overlaps
                // while the current domain is a fluid domain
                URef().registerSideFieldsForInterfaceSide(
                    interf->index(), true, true);
                URef().registerSideFieldsForInterfaceSide(
                    interf->index(), false, true);
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // consider only the side of the interface that belong to
                    // the current fluid domain
                    URef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
                else
                {
                    if (interf->isSlipNonOverlap())
                        continue;

                    // register side fields only if there are non-overlaps while
                    // the current domain is a fluid domain
                    URef().registerSideFieldsForInterfaceSide(
                        interf->index(),
                        interf->isMasterZone(domain->index()),
                        true);
                }
            }
        }
#endif /* HAS_INTERFACE */

        // boundary conditions for this domain
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc = URef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option =
                            "subsonic"; // default: subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "velocity_components")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    if (massAndMomentumNode["velocity"])
                                    {
                                        bc.query<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            "velocity");
                                    }
                                    else
                                    {
#if SPATIAL_DIM == 2
                                        bc.queryMulti<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            {"u", "v"});
#elif SPATIAL_DIM == 3
                                        bc.queryMulti<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            {"u", "v", "w"});
#endif
                                    }
                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "static_pressure")
                                {
                                    bc.setType(boundaryConditionType::
                                                   specifiedDirection);

                                    // query flow direction
                                    if (boundaryDetailsNode["flow_direction"])
                                    {
                                        const auto& flowDirectionNode =
                                            boundaryDetailsNode
                                                ["flow_direction"];

                                        if (flowDirectionNode["option"])
                                        {
                                            auto flowDirectionOption =
                                                flowDirectionNode["option"]
                                                    .template as<std::string>();

                                            if (flowDirectionOption ==
                                                "normal_to_boundary_condition")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "normal_to_boundary_"
                                                    "condition");
                                            }
                                            else if (flowDirectionOption ==
                                                     "cartesian_components")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cartesian_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y", "z"});
#endif
                                                }
                                            }
                                            else if (flowDirectionOption ==
                                                     "cylindrical_components")
                                            {
                                                assert(SPATIAL_DIM == 3);

                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cylindrical_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta", "z"});
#endif
                                                }

                                                // query rotation axis
                                                bc.query<SPATIAL_DIM>(
                                                    flowDirectionNode,
                                                    "rotation_axis",
                                                    "rotation_axis");
                                            }
                                            else
                                            {
                                                errorMsg("Invalid option for "
                                                         "flow_direction");
                                            }
                                        }
                                        else
                                        {
                                            errorMsg("option not provided for "
                                                     "flow direction node");
                                        }
                                    }
                                    else
                                    {
                                        errorMsg("flow_direction node is not "
                                                 "provided for static pressure "
                                                 "inlet");
                                    }

                                    URef().registerSideFlowDirectionFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option == "normal_speed")
                                {
                                    bc.setType(
                                        boundaryConditionType::normalSpeed);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "normal_speed");
                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "total_pressure")
                                {
                                    bc.setType(boundaryConditionType::
                                                   specifiedDirection);

                                    // query flow direction
                                    if (boundaryDetailsNode["flow_direction"])
                                    {
                                        const auto& flowDirectionNode =
                                            boundaryDetailsNode
                                                ["flow_direction"];

                                        if (flowDirectionNode["option"])
                                        {
                                            auto flowDirectionOption =
                                                flowDirectionNode["option"]
                                                    .template as<std::string>();

                                            if (flowDirectionOption ==
                                                "normal_to_boundary_condition")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "normal_to_boundary_"
                                                    "condition");
                                            }
                                            else if (flowDirectionOption ==
                                                     "cartesian_components")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cartesian_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y", "z"});
#endif
                                                }
                                            }
                                            else if (flowDirectionOption ==
                                                     "cylindrical_components")
                                            {
                                                assert(SPATIAL_DIM == 3);

                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cylindrical_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta", "z"});
#endif
                                                }

                                                // query rotation axis
                                                bc.query<SPATIAL_DIM>(
                                                    flowDirectionNode,
                                                    "rotation_axis",
                                                    "rotation_axis");
                                            }
                                            else
                                            {
                                                errorMsg("Invalid option for "
                                                         "flow_direction");
                                            }
                                        }
                                        else
                                        {
                                            errorMsg("option not provided for "
                                                     "flow direction node");
                                        }
                                    }
                                    else
                                    {
                                        errorMsg("flow_direction node is not "
                                                 "provided for total pressure "
                                                 "inlet");
                                    }

                                    URef().registerSideFlowDirectionFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option == "mass_flow_rate")
                                {
                                    bc.setType(
                                        boundaryConditionType::massFlowRate);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "mass_flow_rate");

                                    // query flow direction
                                    if (boundaryDetailsNode["flow_direction"])
                                    {
                                        const auto& flowDirectionNode =
                                            boundaryDetailsNode
                                                ["flow_direction"];

                                        if (flowDirectionNode["option"])
                                        {
                                            auto flowDirectionOption =
                                                flowDirectionNode["option"]
                                                    .template as<std::string>();

                                            if (flowDirectionOption ==
                                                "normal_to_boundary_condition")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "normal_to_boundary_"
                                                    "condition");
                                            }
                                            else if (flowDirectionOption ==
                                                     "cartesian_components")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cartesian_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y", "z"});
#endif
                                                }
                                            }
                                            else if (flowDirectionOption ==
                                                     "cylindrical_components")
                                            {
                                                assert(SPATIAL_DIM == 3);

                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cylindrical_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"r", "theta", "z"});
#endif
                                                }

                                                // query rotation axis
                                                bc.query<SPATIAL_DIM>(
                                                    flowDirectionNode,
                                                    "rotation_axis",
                                                    "rotation_axis");
                                            }
                                            else
                                            {
                                                errorMsg("Invalid option for "
                                                         "flow_direction");
                                            }
                                        }
                                        else
                                        {
                                            errorMsg("option not provided for "
                                                     "flow direction node");
                                        }
                                    }
                                    else
                                    {
                                        errorMsg("flow_direction node is not "
                                                 "provided for mass flow rate "
                                                 "inlet");
                                    }

                                    URef().registerSideFlowDirectionFields(
                                        domain->index(), iBoundary);
                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        std::string("option for ") +
                                        "mass_and_momentum" +
                                        std::string(
                                            " not provided at subsonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum node not provided");
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option ==
                                    "velocity_components_and_static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    if (massAndMomentumNode["velocity"])
                                    {
                                        bc.query<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            "velocity");
                                    }
                                    else
                                    {
#if SPATIAL_DIM == 2
                                        bc.queryMulti<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            {"u", "v"});
#elif SPATIAL_DIM == 3
                                        bc.queryMulti<SPATIAL_DIM>(
                                            massAndMomentumNode,
                                            "value",
                                            {"u", "v", "w"});
#endif
                                    }
                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(std::string("option for ") +
                                             "mass_and_momentum" +
                                             std::string(" not provided at "
                                                         "supersonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum node not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option =
                            "subsonic"; // default: subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "average_static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "mass_flow_rate")
                                {
                                    bc.setType(
                                        boundaryConditionType::massFlowRate);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "mass_flow_rate");
                                }
                                else
                                {
                                    errorMsg(std::string("option for ") +
                                             "mass_and_momentum" +
                                             std::string(" not provided at "
                                                         "subsonic outlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum" +
                                         std::string(" not provided"));
                            }
                        }
                        else if (option == "supersonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for outlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option =
                            "subsonic"; // default: subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "opening_pressure")
                                {
                                    bc.setType(boundaryConditionType::
                                                   specifiedDirection);

                                    // query flow direction
                                    if (boundaryDetailsNode["flow_direction"])
                                    {
                                        const auto& flowDirectionNode =
                                            boundaryDetailsNode
                                                ["flow_direction"];

                                        if (flowDirectionNode["option"])
                                        {
                                            auto flowDirectionOption =
                                                flowDirectionNode["option"]
                                                    .template as<std::string>();

                                            if (flowDirectionOption ==
                                                "normal_to_boundary_condition")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "normal_to_boundary_"
                                                    "condition");
                                            }
                                            else if (flowDirectionOption ==
                                                     "cartesian_components")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cartesian_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y", "z"});
#endif
                                                }
                                            }
                                            else
                                            {
                                                errorMsg("Invalid option for "
                                                         "flow_direction");
                                            }
                                        }
                                        else
                                        {
                                            errorMsg("option not provided for "
                                                     "flow direction node");
                                        }
                                    }
                                    else
                                    {
                                        errorMsg("flow_direction node is not "
                                                 "provided for opening");
                                    }

                                    URef().registerSideFlowDirectionFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option == "static_pressure")
                                {
                                    bc.setType(boundaryConditionType::
                                                   specifiedDirection);

                                    // query flow direction
                                    if (boundaryDetailsNode["flow_direction"])
                                    {
                                        const auto& flowDirectionNode =
                                            boundaryDetailsNode
                                                ["flow_direction"];

                                        if (flowDirectionNode["option"])
                                        {
                                            auto flowDirectionOption =
                                                flowDirectionNode["option"]
                                                    .template as<std::string>();

                                            if (flowDirectionOption ==
                                                "normal_to_boundary_condition")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "normal_to_boundary_"
                                                    "condition");
                                            }
                                            else if (flowDirectionOption ==
                                                     "cartesian_components")
                                            {
                                                bc.addRawData(
                                                    "flow_direction_option",
                                                    "cartesian_components");

                                                if (flowDirectionNode
                                                        ["flow_direction"])
                                                {
                                                    bc.query<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        "flow_direction");
                                                }
                                                else
                                                {
#if SPATIAL_DIM == 2
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y"});
#elif SPATIAL_DIM == 3
                                                    bc.queryMulti<SPATIAL_DIM>(
                                                        flowDirectionNode,
                                                        "flow_direction",
                                                        {"x", "y", "z"});
#endif
                                                }
                                            }
                                            else
                                            {
                                                errorMsg("Invalid option for "
                                                         "flow_direction");
                                            }
                                        }
                                        else
                                        {
                                            errorMsg("option not provided for "
                                                     "flow direction node");
                                        }
                                    }
                                    else
                                    {
                                        errorMsg("flow_direction node is not "
                                                 "provided for opening");
                                    }

                                    URef().registerSideFlowDirectionFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg("option for mass_and_momentum not "
                                             "provided at subsonic opening");
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum node not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for opening not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["mass_and_momentum"])
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];

                            std::string option = "no_slip_wall"; // default
                            if (massAndMomentumNode["option"])
                            {
                                option = massAndMomentumNode["option"]
                                             .template as<std::string>();
                            }

                            if (option == "no_slip_wall")
                            {
                                if (massAndMomentumNode["wall_velocity"])
                                {
                                    bc.setType(boundaryConditionType::noSlip);

                                    if (massAndMomentumNode["wall_velocity"]
                                                           ["option"])
                                    {
                                        std::string wallVelocityOption =
                                            massAndMomentumNode
                                                ["wall_velocity"]["option"]
                                                    .template as<std::string>();
                                        if (wallVelocityOption ==
                                            "cartesian_components")
                                        {
                                            if (massAndMomentumNode
                                                    ["wall_velocity"]
                                                    ["wall_velocity"])
                                            {
                                                bc.query<SPATIAL_DIM>(
                                                    massAndMomentumNode
                                                        ["wall_velocity"],
                                                    "value",
                                                    "wall_velocity");
                                            }
                                            else
                                            {
#if SPATIAL_DIM == 2
                                                bc.queryMulti<SPATIAL_DIM>(
                                                    massAndMomentumNode
                                                        ["wall_velocity"],
                                                    "value",
                                                    {"wall_u", "wall_v"});
#elif SPATIAL_DIM == 3
                                                bc.queryMulti<SPATIAL_DIM>(
                                                    massAndMomentumNode
                                                        ["wall_velocity"],
                                                    "value",
                                                    {"wall_u",
                                                     "wall_v",
                                                     "wall_w"});
#endif
                                            }
                                        }
                                        else if (wallVelocityOption ==
                                                 "rotating_wall")
                                        {
                                            scalar omega = 0.0;
                                            std::vector<scalar> axis(
                                                SPATIAL_DIM, 0);
                                            std::vector<scalar> origin(
                                                SPATIAL_DIM, 0);

                                            // get angular velocity
                                            if (massAndMomentumNode
                                                    ["wall_velocity"]
                                                    ["angular_velocity"])
                                            {
                                                omega = massAndMomentumNode
                                                            ["wall_velocity"]
                                                            ["angular_velocity"]
                                                                .template as<
                                                                    scalar>();
                                            }
                                            else
                                            {
                                                errorMsg(
                                                    "angular_velocity is not "
                                                    "set for rotating_wall");
                                            }

                                            // get rotation axis
                                            if (massAndMomentumNode
                                                    ["wall_velocity"]
                                                    ["rotation_axis"])
                                            {
                                                axis =
                                                    massAndMomentumNode
                                                        ["wall_velocity"]
                                                        ["rotation_axis"]
                                                            .template as<
                                                                std::vector<
                                                                    scalar>>();
                                            }
                                            else
                                            {
                                                errorMsg(
                                                    "rotation_axis is not set "
                                                    "for rotating_wall");
                                            }

                                            // get rotation origin
                                            if (massAndMomentumNode
                                                    ["wall_velocity"]["origin"])
                                            {
                                                origin =
                                                    massAndMomentumNode
                                                        ["wall_velocity"]
                                                        ["origin"]
                                                            .template as<
                                                                std::vector<
                                                                    scalar>>();
                                            }
                                            else
                                            {
                                                errorMsg("origin is not set "
                                                         "for rotating_wall");
                                            }

                                            // construct the expression
                                            std::vector<std::string> expression(
                                                SPATIAL_DIM);
#if SPATIAL_DIM == 3
                                            expression[0] =
                                                std::to_string(omega) + " * (" +
                                                std::to_string(axis[1]) +
                                                " * (z - " +
                                                std::to_string(origin[2]) +
                                                ") - " +
                                                std::to_string(axis[2]) +
                                                " * (y - " +
                                                std::to_string(origin[1]) +
                                                "))";

                                            expression[1] =
                                                std::to_string(omega) + " * (" +
                                                std::to_string(axis[2]) +
                                                " * (x - " +
                                                std::to_string(origin[0]) +
                                                ") - " +
                                                std::to_string(axis[0]) +
                                                " * (z - " +
                                                std::to_string(origin[2]) +
                                                "))";

                                            expression[2] =
                                                std::to_string(omega) + " * (" +
                                                std::to_string(axis[0]) +
                                                " * (y - " +
                                                std::to_string(origin[1]) +
                                                ") - " +
                                                std::to_string(axis[1]) +
                                                " * (x - " +
                                                std::to_string(origin[0]) +
                                                "))";
#elif SPATIAL_DIM == 2
                                            expression[0] =
                                                std::to_string(omega) + " * (" +
                                                std::to_string(0.0) +
                                                " * (z - " +
                                                std::to_string(origin[2]) +
                                                ") - " + std::to_string(1.0) +
                                                " * (y - " +
                                                std::to_string(origin[1]) +
                                                "))";

                                            expression[1] =
                                                std::to_string(omega) + " * (" +
                                                std::to_string(1.0) +
                                                " * (x - " +
                                                std::to_string(origin[0]) +
                                                ") - " + std::to_string(0.0) +
                                                " * (z - " +
                                                std::to_string(origin[2]) +
                                                "))";
#endif
                                            bc.addExpression<SPATIAL_DIM>(
                                                "value", expression);
                                        }
                                        else if (wallVelocityOption ==
                                                 "counter_rotating_wall")
                                        {
                                            // this option must only be used
                                            // in the case of a rotating
                                            // domain
                                            assert(domain->zoneRef()
                                                       .frameRotating() ||
                                                   domain->zoneRef()
                                                       .meshTransforming());
                                            assert(domain->zoneRef()
                                                       .transformationRef()
                                                       .type() ==
                                                   meshMotionType::rotating);

                                            // easiest way to const_cast ..
                                            // otherwise will be a mess
                                            dynamic_cast<wall&>(
                                                const_cast<::accel::domain*>(
                                                    domain)
                                                    ->zonePtr()
                                                    ->boundaryRef(iBoundary))
                                                .setCounterRotating(true);
                                            bc.setConstantValue<SPATIAL_DIM>(
                                                "value",
                                                std::vector<scalar>(SPATIAL_DIM,
                                                                    0));
                                        }
                                        else
                                        {
                                            errorMsg("Invalid wall_velocity "
                                                     "option");
                                        }
                                    }
                                    else
                                    {
                                        // cartersian components by default
                                        bc.queryConstantValue<SPATIAL_DIM>(
                                            massAndMomentumNode
                                                ["wall_velocity"],
                                            "value",
                                            "value");
                                    }

                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    bc.setType(boundaryConditionType::noSlip);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "value",
                                        std::vector<scalar>(SPATIAL_DIM, 0));
                                    URef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                            }
                            else if (option == "free_slip_wall")
                            {
                                bc.setType(boundaryConditionType::slip);
                            }
                            else
                            {
                                bc.setType(boundaryConditionType::noSlip);
                                bc.setConstantValue<SPATIAL_DIM>(
                                    "value",
                                    std::vector<scalar>(SPATIAL_DIM, 0));
                                URef().registerSideFields(domain->index(),
                                                          iBoundary);
                            }

                            if (massAndMomentumNode
                                    ["wall_velocity_relative_to"])
                            {
                                if (massAndMomentumNode
                                        ["wall_velocity_relative_to"]
                                            .template as<std::string>() ==
                                    "mesh_motion")
                                {
                                    // easiest way to const_cast ..
                                    // otherwise will be a mess
                                    dynamic_cast<wall&>(
                                        const_cast<::accel::domain*>(domain)
                                            ->zonePtr()
                                            ->boundaryRef(iBoundary))
                                        .setWallVelocityRelativeToMeshMotion(
                                            true);
                                }
                                else
                                {
                                    errorMsg("Not yet implemented");
                                }
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::noSlip);
                            bc.setConstantValue<SPATIAL_DIM>(
                                "value", std::vector<scalar>(SPATIAL_DIM, 0));
                            URef().registerSideFields(domain->index(),
                                                      iBoundary);
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::symmetry);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupPressure(const std::shared_ptr<domain> domain)
{
    if (pRef().isZoneUnset(domain->index()))
    {
        pRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            pRef(), realm::p_ID, domain);

        // boundary conditions
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc = pRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option = "subsonic"; // default
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "velocity_components")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::staticPressure);
                                    // relative pressure here is the relative
                                    // gauge pressure
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "normal_speed")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "mass_flow_rate")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "total_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::totalPressure);
                                    // relative pressure here is the relative
                                    // total pressure
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        std::string("option for ") +
                                        "mass_and_momentum" +
                                        std::string(
                                            " not provided at subsonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum" +
                                         std::string(" not provided"));
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option ==
                                    "velocity_components_and_static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::staticPressure);
                                    // relative pressure here is the relative
                                    // gauge pressure
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        std::string("option for ") +
                                        "mass_and_momentum" +
                                        std::string(
                                            " not provided at subsonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum" +
                                         std::string(" not provided"));
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option = "subsonic"; // default
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::staticPressure);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "average_static_pressure")
                                {
                                    bc.setType(boundaryConditionType::
                                                   averageStaticPressure);
                                    bc.query<1>(massAndMomentumNode,
                                                "average_static_pressure",
                                                "relative_pressure");
                                    bc.query<1>(massAndMomentumNode,
                                                "blend_factor",
                                                "pressure_profile_blend");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "normal_speed")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else if (option == "mass_flow_rate")
                                {
                                    bc.setType(
                                        boundaryConditionType::massFlowRate);
                                }
                                else
                                {
                                    errorMsg(std::string("option for ") +
                                             "mass_and_momentum" +
                                             std::string(" not provided at "
                                                         "subsonic outlet"));
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum" +
                                         std::string(" not provided"));
                            }
                        }
                        else if (option == "supersonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option = "subsonic"; // default
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];
                            if (massAndMomentumNode)
                            {
                                std::string option =
                                    massAndMomentumNode["option"]
                                        .template as<std::string>();

                                if (option == "static_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::staticPressure);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "opening_pressure")
                                {
                                    bc.setType(
                                        boundaryConditionType::totalPressure);
                                    bc.query<1>(massAndMomentumNode,
                                                "value",
                                                "relative_pressure");
                                    pRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg("option for mass_and_momentum not "
                                             "provided at subsonic opening");
                                }
                            }
                            else
                            {
                                errorMsg("mass_and_momentum not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for opening not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupTemperature(const std::shared_ptr<domain> domain)
{
    if (TRef().isZoneUnset(domain->index()))
    {
        TRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            TRef(), realm::T_ID, domain);

        // boundary conditions for this domain
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc = TRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["heat_transfer"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "static_temperature")
                                {
                                    bc.setType(boundaryConditionType::
                                                   staticTemperature);
                                    bc.query<1>(
                                        node, "value", "static_temperature");
                                    TRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "total_temperature")
                                {
                                    bc.setType(boundaryConditionType::
                                                   totalTemperature);
                                    bc.query<1>(
                                        node, "value", "total_temperature");
                                    TRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        std::string("option for ") +
                                        "heat_transfer" +
                                        std::string(
                                            " not provided at subsonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("heat_transfer" +
                                         std::string(" not provided"));
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["heat_transfer"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "static_temperature")
                                {
                                    bc.setType(boundaryConditionType::
                                                   staticTemperature);
                                    bc.query<1>(
                                        node, "value", "static_temperature");
                                    TRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(std::string("option for ") +
                                             "heat_transfer" +
                                             std::string(" not provided at "
                                                         "supersonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("heat_transfer" +
                                         std::string(" not provided"));
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else if (option == "supersonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["heat_transfer"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "static_temperature")
                                {
                                    bc.setType(boundaryConditionType::
                                                   staticTemperature);
                                    bc.query<1>(
                                        node, "value", "static_temperature");
                                    TRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "opening_temperature")
                                {
                                    bc.setType(boundaryConditionType::
                                                   totalTemperature);
                                    bc.query<1>(
                                        node, "value", "opening_temperature");
                                    TRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        std::string("option for ") +
                                        "heat_transfer" +
                                        std::string(
                                            " not provided at subsonic inlet"));
                                }
                            }
                            else
                            {
                                errorMsg("heat_transfer" +
                                         std::string(" not provided"));
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["heat_transfer"])
                        {
                            const auto& node =
                                boundaryDetailsNode["heat_transfer"];
                            std::string option =
                                node["option"].template as<std::string>();

                            if (option == "adiabatic")
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                            }
                            else if (option == "temperature")
                            {
                                bc.setType(
                                    boundaryConditionType::specifiedValue);
                                bc.query<1>(node, "value", "fixed_temperature");
                                TRef().registerSideFields(domain->index(),
                                                          iBoundary);
                            }
                            else if (option == "heat_flux")
                            {
                                bc.setType(
                                    boundaryConditionType::specifiedFlux);
                                bc.query<1>(node, "flux", "heat_flux_in");
                                TRef().registerSideFluxField(domain->index(),
                                                             iBoundary);
                            }
                            else
                            {
                                errorMsg(std::string("option for ") +
                                         "heat_transfer" +
                                         std::string(" not provided at wall"));
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupSpecificEnthalpy(const std::shared_ptr<domain> domain)
{
    if (hRef().isZoneUnset(domain->index()))
    {
        hRef().setZone(domain->index());

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            auto& bc = h0Ref().boundaryConditionRef(domain->index(), iBoundary);

            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();
            boundaryConditionType bcType =
                TRef().boundaryConditionRef(domain->index(), iBoundary).type();

            switch (type)
            {
                case boundaryPhysicalType::wall:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::specifiedValue:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    hRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                case boundaryPhysicalType::inlet:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::staticTemperature:
                            case boundaryConditionType::totalTemperature:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    hRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::staticTemperature:
                            case boundaryConditionType::totalTemperature:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    hRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }
}

void fieldBroker::setupSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    if (h0Ref().isZoneUnset(domain->index()))
    {
        h0Ref().setZone(domain->index());

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            auto& bc = h0Ref().boundaryConditionRef(domain->index(), iBoundary);

            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();
            boundaryConditionType bcType =
                TRef().boundaryConditionRef(domain->index(), iBoundary).type();

            switch (type)
            {
                case boundaryPhysicalType::wall:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::specifiedValue:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    h0Ref().registerSideFields(domain->index(),
                                                               iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                case boundaryPhysicalType::inlet:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::staticTemperature:
                            case boundaryConditionType::totalTemperature:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    h0Ref().registerSideFields(domain->index(),
                                                               iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::staticTemperature:
                            case boundaryConditionType::totalTemperature:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    h0Ref().registerSideFields(domain->index(),
                                                               iBoundary);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }
}

void fieldBroker::setupWallScale(const std::shared_ptr<domain> domain)
{
    if (yScaleRef().isZoneUnset(domain->index()))
    {
        yScaleRef().setZone(domain->index());

#ifdef HAS_INTERFACE
        // 1) Register yscale side fields on fluid-solid interface (only fluid)
        //    side, because we deal with this as a typical no-slip wall
        // 2) Register also on an interface side which has a non-overlap portion
        //    because the non-overlap portion shall be treated as a no-slip wall
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // consider only the side of the interface that belong to
                    // the current domain
                    yScaleRef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                    interf->interfaceSideInfoPtr(domain->index())
                        ->dataHandlerRef()
                        .setConstantValue<1>(yScaleRef().name() + "_value",
                                             {0});
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            auto& bc =
                yScaleRef().boundaryConditionRef(domain->index(), iBoundary);

            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();
            boundaryConditionType bcType =
                URef().boundaryConditionRef(domain->index(), iBoundary).type();

            switch (type)
            {
                case boundaryPhysicalType::wall:
                    {
                        switch (bcType)
                        {
                            case boundaryConditionType::noSlip:
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.setConstantValue<1>("value", {0});
                                    yScaleRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                break;

                            default:
                                bc.setType(boundaryConditionType::zeroGradient);
                                break;
                        }
                    }
                    break;

                default:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;
            }
        }
    }
}

void fieldBroker::setupTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    if (kRef().isZoneUnset(domain->index()))
    {
        kRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            kRef(), turbRealm::k_ID, domain);

#ifdef HAS_INTERFACE
        // 1) Register k side fields on fluid-solid interface (only fluid)
        //    side, because we deal with this as a typical no-slip wall
        // 2) Register also on an interface side which has a non-overlap
        // portion
        //    because the non-overlap portion shall be treated as a no-slip
        //    wall
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // consider only the side of the interface that belong to
                    // the current domain
                    kRef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        // boundary conditions
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc = kRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega" ||
                                    option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "k");
                                    kRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option == "intensity_and_length_scale")
                                {
                                    bc.setType(boundaryConditionType::
                                                   intensityAndLengthScale);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(node, "l", "eddy_length_scale");
                                    kRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (option ==
                                         "intensity_and_eddy_viscosity_ratio")
                                {
                                    bc.setType(
                                        boundaryConditionType::
                                            intensityAndEddyViscosityRatio);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(
                                        node, "r", "eddy_viscosity_ratio");
                                    kRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega" ||
                                    option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "k");
                                    kRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega" ||
                                    option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "k");
                                    kRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["mass_and_momentum"])
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];

                            std::string option = "no_slip_wall"; // default
                            if (massAndMomentumNode["option"])
                            {
                                option = massAndMomentumNode["option"]
                                             .template as<std::string>();
                            }

                            if (option == "no_slip_wall")
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                                kRef().registerSideFields(
                                    domain->index(),
                                    iBoundary); // required for wall-function
                            }
                            else
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                            kRef().registerSideFields(
                                domain->index(),
                                iBoundary); // required for wall-function
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    if (omegaRef().isZoneUnset(domain->index()))
    {
        omegaRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            omegaRef(), turbRealm::omega_ID, domain);

#ifdef HAS_INTERFACE
        // 1) Register omega side fields on fluid-solid interface (only
        // fluid)
        //    side, because we deal with this as a typical no-slip wall
        // 2) Register also on an interface side which has a non-overlap
        // portion
        //    because the non-overlap portion shall be treated as a no-slip
        //    wall
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // consider only the side of the interface that belong to
                    // the current domain
                    omegaRef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        // boundary conditions for this domain
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc =
                omegaRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "omega");
                                    omegaRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option == "intensity_and_length_scale")
                                {
                                    bc.setType(boundaryConditionType::
                                                   intensityAndLengthScale);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(node, "l", "eddy_length_scale");
                                    omegaRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option ==
                                         "intensity_and_eddy_viscosity_ratio")
                                {
                                    bc.setType(
                                        boundaryConditionType::
                                            intensityAndEddyViscosityRatio);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(
                                        node, "r", "eddy_viscosity_ratio");
                                    omegaRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "omega");
                                    omegaRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_omega")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "omega");
                                    omegaRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["mass_and_momentum"])
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];

                            std::string option = "no_slip_wall"; // default
                            if (massAndMomentumNode["option"])
                            {
                                option = massAndMomentumNode["option"]
                                             .template as<std::string>();
                            }

                            if (option == "no_slip_wall")
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                                omegaRef().registerSideFields(
                                    domain->index(),
                                    iBoundary); // required for wall-function
                            }
                            else
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                            omegaRef().registerSideFields(
                                domain->index(),
                                iBoundary); // required for wall-function
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    if (epsilonRef().isZoneUnset(domain->index()))
    {
        epsilonRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            epsilonRef(), turbRealm::epsilon_ID, domain);

#ifdef HAS_INTERFACE
        // 1) Register epsilon side fields on fluid-solid interface (only fluid)
        //    side, because we deal with this as a typical no-slip wall
        // 2) Register also on an interface side which has a non-overlap
        // portion
        //    because the non-overlap portion shall be treated as a no-slip
        //    wall
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // consider only the side of the interface that belong to
                    // the current domain
                    epsilonRef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        // boundary conditions
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc =
                epsilonRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "epsilon");
                                    epsilonRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option == "intensity_and_length_scale")
                                {
                                    bc.setType(boundaryConditionType::
                                                   intensityAndLengthScale);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(node, "l", "eddy_length_scale");
                                    epsilonRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else if (option ==
                                         "intensity_and_eddy_viscosity_ratio")
                                {
                                    bc.setType(
                                        boundaryConditionType::
                                            intensityAndEddyViscosityRatio);
                                    bc.query<1>(
                                        node, "I", "fractional_intensity");
                                    bc.query<1>(
                                        node, "r", "eddy_viscosity_ratio");
                                    epsilonRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else if (option == "supersonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "epsilon");
                                    epsilonRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::outlet:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                        else
                        {
                            errorMsg("option for inlet not available");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        std::string option =
                            "subsonic"; // default subsonic/supersonic
                        if (boundaryDetailsNode["option"])
                        {
                            option = boundaryDetailsNode["option"]
                                         .template as<std::string>();
                        }

                        if (option == "subsonic")
                        {
                            const auto& node =
                                boundaryDetailsNode["turbulence"];
                            if (node)
                            {
                                std::string option =
                                    node["option"].template as<std::string>();

                                if (option == "k_and_epsilon")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<1>(node, "value", "epsilon");
                                    epsilonRef().registerSideFields(
                                        domain->index(), iBoundary);
                                }
                                else
                                {
                                    errorMsg(
                                        "option for turbulence not provided "
                                        "at subsonic inlet");
                                }
                            }
                            else
                            {
                                errorMsg("turbulence not provided");
                            }
                        }
                        else
                        {
                            errorMsg("option for inlet not implemented yet");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["mass_and_momentum"])
                        {
                            const auto& massAndMomentumNode =
                                boundaryDetailsNode["mass_and_momentum"];

                            std::string option = "no_slip_wall"; // default
                            if (massAndMomentumNode["option"])
                            {
                                option = massAndMomentumNode["option"]
                                             .template as<std::string>();
                            }

                            if (option == "no_slip_wall")
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                                epsilonRef().registerSideFields(
                                    domain->index(),
                                    iBoundary); // required for wall-function
                            }
                            else
                            {
                                bc.setType(boundaryConditionType::zeroGradient);
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                            epsilonRef().registerSideFields(
                                domain->index(),
                                iBoundary); // required for wall-function
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupTransitionOnsetReynoldsNumber(
    const std::shared_ptr<domain> domain)
{
    if (ReThetaRef().isZoneUnset(domain->index()))
    {
        ReThetaRef().setZone(domain->index());

        initialCondition::setupFieldInitializationOverDomainFromValues(
            ReThetaRef(), domain->index(), {1.0});

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            auto& bc =
                ReThetaRef().boundaryConditionRef(domain->index(), iBoundary);

            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::opening:
                    {
                        bc.setType(boundaryConditionType::specifiedValue);

                        // values of ReTheta at inlet/opening will always be
                        // deduced from velocity and turbulent kinetic energy.
                        // This is done in the transition SST model

                        ReThetaRef().registerSideFields(domain->index(),
                                                        iBoundary);
                    }
                    break;

                default:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;
            }
        }
    }
}

void fieldBroker::setupTurbulentIntermittency(
    const std::shared_ptr<domain> domain)
{
    if (gammaRef().isZoneUnset(domain->index()))
    {
        gammaRef().setZone(domain->index());

        initialCondition::setupFieldInitializationOverDomainFromValues(
            gammaRef(), domain->index(), {1.0});

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            auto& bc =
                gammaRef().boundaryConditionRef(domain->index(), iBoundary);

            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::opening:
                    {
                        bc.setType(boundaryConditionType::specifiedValue);
                        bc.setConstantValue<1>("value", {1.0});
                        gammaRef().registerSideFields(domain->index(),
                                                      iBoundary);
                    }
                    break;

                default:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;
            }
        }
    }
}

void fieldBroker::setupYoungModulus(const std::shared_ptr<domain> domain)
{
    if (ERef().isZoneUnset(domain->index()))
    {
        SET_YOUNG_MODULUS_PROPERTY(this->ERef(), domain->getYAMLMaterial());
    }
}

void fieldBroker::setupPoissonRatio(const std::shared_ptr<domain> domain)
{
    if (nuRef().isZoneUnset(domain->index()))
    {
        SET_POISSON_RATIO_PROPERTY(this->nuRef(), domain->getYAMLMaterial());
    }
}

void fieldBroker::setupDensity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        SET_DENSITY_PROPERTY(this->rhoRef(), domain->getYAMLMaterial());

        // strong side fields?
        // this is the case for supersonic inlets
        if (domain->type() == domainType::fluid &&
            domain->isMaterialCompressible())
        {
            for (label iBoundary = 0;
                 iBoundary < domain->zonePtr()->nBoundaries();
                 iBoundary++)
            {
                auto type = domain->zonePtr()->boundaryPtr(iBoundary)->type();
                auto UBCType =
                    URef()
                        .boundaryConditionRef(domain->index(), iBoundary)
                        .type();
                auto pBCType =
                    pRef()
                        .boundaryConditionRef(domain->index(), iBoundary)
                        .type();
                switch (type)
                {
                    case boundaryPhysicalType::inlet:
                        {
                            switch (UBCType)
                            {
                                case boundaryConditionType::specifiedValue:
                                    {
                                        switch (pBCType)
                                        {
                                            case boundaryConditionType::
                                                staticPressure:
                                                {
                                                    rhoRef().registerSideFields(
                                                        domain->index(),
                                                        iBoundary);
                                                }
                                                break;

                                            default:
                                                break;
                                        }
                                    }
                                    break;

                                default:
                                    break;
                            }
                        }
                        break;

                    default:
                        break;
                }
            }
        }
    }
}

void fieldBroker::setupDynamicViscosity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        SET_DYNAMIC_VISCOSITY_PROPERTY(this->muRef(),
                                       domain->getYAMLMaterial());
    }
}

void fieldBroker::setupSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        SET_SPECIFIC_HEAT_CAPACITY_PROPERTY(this->cpRef(),
                                            domain->getYAMLMaterial());
    }
}

void fieldBroker::setupThermalConductivity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        SET_THERMAL_CONDUCTIVITY_PROPERTY(this->lambdaRef(),
                                          domain->getYAMLMaterial());
    }
}

void fieldBroker::setupThermalExpansivity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->beta_)
    {
        SET_THERMAL_EXPANSIVITY_PROPERTY(this->betaRef(),
                                         domain->getYAMLMaterial());
    }
}

void fieldBroker::setupCompressibility(const std::shared_ptr<domain> domain)
{
    if (psiRef().isZoneUnset(domain->index()))
    {
        psiRef().setZone(domain->index());
    }
}

void fieldBroker::setupDisplacement(const std::shared_ptr<domain> domain)
{
    if (DRef().isZoneUnset(domain->index()))
    {
        DRef().setZone(domain->index());

        // set initial conditions for mesh displacement to 0
        initialCondition::setupFieldInitializationOverDomainFromValues(
            DRef(), domain->index(), std::vector<scalar>(SPATIAL_DIM, 0.0));

#ifdef HAS_INTERFACE
        // for a fluid-solid interface, the interface side in this domain will
        // act as a dirichlet boundary for the displacement diffusion equation,
        // which is being solved over this domain, thus, side fields on the
        // this side must be registered
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // this is a fluid-structure interface. Displacement
                    // side flux field must be defined on this side
                    DRef().registerSideFluxFieldForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        setupBoundaryConditions_(domain,
                                 // anonymous function to set boundary
                                 // conditions for this model
                                 [this](const ::accel::domain* domain,
                                        const label iBoundary,
                                        const boundaryPhysicalType bc_type,
                                        const YAML::Node boundaryDetailsNode)
        {
            auto& bc = DRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
 std::cout << "Setting boundary conditions:\n";
 std::cout << "\tdomain name: " << domain->name() << "\n";
 std::cout << "\tdomain index: " << domain->index() << "\n";
 std::cout << "\tpatch index: " << iBoundary << "\n";
 std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
 std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::symmetry);
                    break;

                default:
                    {
                        if (boundaryDetailsNode["mesh_motion"])
                        {
                            const auto& meshMotionNode =
                                boundaryDetailsNode["mesh_motion"];

                            std::string option =
                                meshMotionNode["option"]
                                    .template as<std::string>();

                            if (option == "specified_displacement")
                            {
                                const auto& displacementNode =
                                    meshMotionNode["displacement"];

                                if (displacementNode["option"])
                                {
                                    std::string displacementOption =
                                        displacementNode["option"]
                                            .template as<std::string>();

                                    if (displacementOption ==
                                        "cartesian_components")
                                    {
                                        bc.setType(boundaryConditionType::
                                                       specifiedValue);

                                        if (displacementNode
                                                ["cartesian_components"])
                                        {
                                            bc.query<SPATIAL_DIM>(
                                                displacementNode,
                                                "value",
                                                "cartesian_components");
                                        }
                                        else
                                        {
#if SPATIAL_DIM == 2
                                            bc.queryMulti<SPATIAL_DIM>(
                                                displacementNode,
                                                "value",
                                                {"x_component", "y_component"});
#elif SPATIAL_DIM == 3
                                            bc.queryMulti<SPATIAL_DIM>(
                                                displacementNode,
                                                "value",
                                                {"x_component",
                                                 "y_component",
                                                 "z_component"});
#endif
                                        }

                                        DRef().registerSideFields(
                                            domain->index(), iBoundary);
                                    }
                                    else
                                    {
                                        errorMsg("invalid option for "
                                                 "displacement node");
                                    }
                                }
                                else
                                {
                                    errorMsg("option for "
                                             "specified_displacment "
                                             "not provided "
                                             "for mesh motion");
                                }
                            }
                            else if (option == "periodic_displacement")
                            {
                                bc.setType(boundaryConditionType::
                                               periodicDisplacement);

                                const auto& displacementNode =
                                    meshMotionNode["displacement"];

                                if (displacementNode["frequency"])
                                {
                                    bc.query<1>(displacementNode,
                                                "frequency",
                                                "frequency");
                                }
                                else
                                {
                                    errorMsg("frequency not "
                                             "provided for "
                                             "displacement");
                                }

                                if (displacementNode["value"])
                                {
                                    bc.query<SPATIAL_DIM>(
                                        displacementNode, "value", "value");
                                }
                                else
                                {
                                    errorMsg("value not provided for "
                                             "displacement node");
                                }

                                DRef().registerSideFields(domain->index(),
                                                          iBoundary);
                            }
                            else if (option == "rigid_body_solution")
                            {
                                bc.setType(
                                    boundaryConditionType::rigidBodySolution);

                                bc.addRawData("rigid_body",
                                              meshMotionNode["rigid_body"]
                                                  .template as<std::string>());

                                DRef().registerSideFields(domain->index(),
                                                          iBoundary);
                            }
                            else
                            {
                                errorMsg("invalid option for "
                                         "mesh_motion node");
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::specifiedValue);
                            bc.setConstantValue<SPATIAL_DIM>(
                                "value", std::vector<scalar>(SPATIAL_DIM, 0));
                            DRef().registerSideFields(domain->index(),
                                                      iBoundary);
                        }
                    }
                    break;
            }
        });
    }
}

void fieldBroker::setupTotalDisplacement(const std::shared_ptr<domain> domain)
{
    if (DtRef().isZoneUnset(domain->index()))
    {
        DtRef().setZone(domain->index());

        // set initial conditions for total mesh displacement to 0
        initialCondition::setupFieldInitializationOverDomainFromValues(
            DtRef(), domain->index(), std::vector<scalar>(SPATIAL_DIM, 0.0));
    }
}

void fieldBroker::setupMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (UmRef().isZoneUnset(domain->index()))
    {
        UmRef().setZone(domain->index());

        // set initial conditions for total mesh displacement to 0
        initialCondition::setupFieldInitializationOverDomainFromValues(
            UmRef(), domain->index(), std::vector<scalar>(SPATIAL_DIM, 0.0));
    }
}

void fieldBroker::setupDivMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (divUmRef().isZoneUnset(domain->index()))
    {
        divUmRef().setZone(domain->index());

        // set initial conditions for total mesh displacement to 0
        initialCondition::setupFieldInitializationOverDomainFromValues(
            divUmRef(), domain->index(), {0.0});
    }
}

void fieldBroker::setupMassFlowRate(const std::shared_ptr<domain> domain)
{
    if (this->mDotRef().isZoneUnset(domain->index()))
    {
        this->mDotRef().setZone(domain->index());
        this->mDotRef().divRef().setZone(domain->index());

#ifdef HAS_INTERFACE
        // register mass flux side field for interfaces in fluid domain
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                this->mDotRef().registerSideFieldsForInterfaceSide(
                    interf->index(), true);
                this->mDotRef().registerSideFieldsForInterfaceSide(
                    interf->index(), false);
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // do nothing
                }
                else
                {
                    this->mDotRef().registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
            }
        }
#endif /* HAS_INTERFACE */

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::opening:
                case boundaryPhysicalType::outlet:
                    this->mDotRef().registerSideField(domain->index(),
                                                      iBoundary);
                    break;

                default:
                    break;
            }
        }
    }
}

void fieldBroker::setupHeatFlowRate(const std::shared_ptr<domain> domain)
{
    if (this->qDotRef().isZoneUnset(domain->index()))
    {
        this->qDotRef().setZone(domain->index());

#ifdef HAS_INTERFACE
        // register heat flux side field for interfaces
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                this->qDotRef().registerSideFieldsForInterfaceSide(
                    interf->index(), true);
                this->qDotRef().registerSideFieldsForInterfaceSide(
                    interf->index(), false);
            }
            else
            {
                this->qDotRef().registerSideFieldsForInterfaceSide(
                    interf->index(), interf->isMasterZone(domain->index()));
            }
        }
#endif /* HAS_INTERFACE */

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            this->qDotRef().registerSideField(domain->index(), iBoundary);
        }
    }
}

void fieldBroker::initializeVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().initialize(domain->index());
    }
}

void fieldBroker::initializePressure(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTemperature(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().initialize(domain->index());
    }
}

void fieldBroker::initializeSpecificEnthalpy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().initialize(domain->index());
    }
}

void fieldBroker::initializeSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().initialize(domain->index());
    }
}

void fieldBroker::initializeWallScale(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTransitionOnsetReynoldsNumber(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTurbulentIntermittency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().initialize(domain->index());
    }
}

void fieldBroker::initializeYoungModulus(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->smRealm_->E_)
    {
        ERef().initialize(domain->index());
    }
}

void fieldBroker::initializePoissonRatio(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->smRealm_->nu_)
    {
        nuRef().initialize(domain->index());
    }
}

void fieldBroker::initializeDensity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().initialize(domain->index());
    }
}

void fieldBroker::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().initialize(domain->index());
    }
}

void fieldBroker::initializeSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().initialize(domain->index());
    }
}

void fieldBroker::initializeThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().initialize(domain->index());
    }
}

void fieldBroker::initializeThermalExpansivity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->beta_)
    {
        betaRef().initialize(domain->index());
    }
}

void fieldBroker::initializeCompressibility(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->psi_)
    {
        psiRef().initialize(domain->index());
    }
}

void fieldBroker::initializeDisplacement(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().initialize(domain->index());
    }
}

void fieldBroker::initializeTotalDisplacement(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().initialize(domain->index());
    }
}

void fieldBroker::initializeMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Um_)
    {
        UmRef().initialize(domain->index());
    }
}

void fieldBroker::initializeDivMeshVelocity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->divUm_)
    {
        divUmRef().initialize(domain->index());
    }
}

void fieldBroker::initializeMassFlowRate(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mDot_)
    {
        // interior
        initializeMassFlowRateInterior_(domain);

#ifdef HAS_INTERFACE
        // Interfaces
        for (const interface* interf : domain->zonePtr()->interfacesRef())
        {
            if (interf->isConformalTreatment())
                continue;

            if (interf->isInternal())
            {
                initializeMassFlowRateInterfaceSideField_(
                    domain, interf->masterInfoPtr());
                initializeMassFlowRateInterfaceSideField_(
                    domain, interf->slaveInfoPtr());
            }
            else if (!interf->isFluidSolidType())
            {
                // get interface side that is sitting in this domain
                const auto* interfaceSideInfoPtr =
                    interf->interfaceSideInfoPtr(domain->index());

                initializeMassFlowRateInterfaceSideField_(domain,
                                                          interfaceSideInfoPtr);
            }
        }
#endif /* HAS_INTERFACE */

        // Boundary
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            initializeMassFlowRateBoundaryField_(
                domain, domain->zonePtr()->boundaryPtr(iBoundary));
        }
    }
}

void fieldBroker::resetVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetPressure(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTemperature(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetSpecificEnthalpy(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().initialize(domain->index(), true);
    }
}

void fieldBroker::resetWallScale(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTransitionOnsetReynoldsNumber(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetTurbulentIntermittency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetDensity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetDynamicViscosity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetThermalConductivity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetThermalExpansivity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->beta_)
    {
        betaRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetCompressibility(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->psi_)
    {
        psiRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetDisplacement(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().initializeField(domain->index());
    }
}

void fieldBroker::resetTotalDisplacement(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().initializeField(domain->index());
    }
}

void fieldBroker::resetMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Um_)
    {
        UmRef().initialize(domain->index(), true);
    }
}

void fieldBroker::resetDivMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->divUm_)
    {
        divUmRef().initialize(domain->index(), true);
    }
}

void fieldBroker::updateVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().update(domain->index());
    }
}

void fieldBroker::updatePressure(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().update(domain->index());
    }
}

void fieldBroker::updateTemperature(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().update(domain->index());
    }
}

void fieldBroker::updateSpecificEnthalpy(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().update(domain->index());
    }
}

void fieldBroker::updateSpecificTotalEnthalpy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().update(domain->index());
    }
}

void fieldBroker::updateWallScale(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().update(domain->index());
    }
}

void fieldBroker::updateTurbulentKineticEnergy(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().update(domain->index());
    }
}

void fieldBroker::updateTurbulentEddyFrequency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().update(domain->index());
    }
}

void fieldBroker::updateTurbulentDissipationRate(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().update(domain->index());
    }
}

void fieldBroker::updateTransitionOnsetReynoldsNumber(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().update(domain->index());
    }
}

void fieldBroker::updateTurbulentIntermittency(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().update(domain->index());
    }
}

void fieldBroker::updateYoungModulus(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->smRealm_->E_)
    {
        ERef().update(domain->index());
    }
}

void fieldBroker::updatePoissonRatio(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->smRealm_->nu_)
    {
        nuRef().update(domain->index());
    }
}

void fieldBroker::updateDensity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().update(domain->index());
    }
}

void fieldBroker::updateDynamicViscosity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().update(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().update(domain->index());
    }
}

void fieldBroker::updateThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().update(domain->index());
    }
}

void fieldBroker::updateThermalExpansivity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->beta_)
    {
        betaRef().update(domain->index());
    }
}

void fieldBroker::updateCompressibility(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->psi_)
    {
        psiRef().update(domain->index());
    }
}

void fieldBroker::updateDisplacement(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().update(domain->index());
    }
}

void fieldBroker::updateTotalDisplacement(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().update(domain->index());
    }
}

void fieldBroker::updateMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Um_)
    {
        UmRef().update(domain->index());
    }
}

void fieldBroker::updateDivMeshVelocity(const std::shared_ptr<domain> domain)
{
    if (realmPtr_->divUm_)
    {
        divUmRef().update(domain->index());
    }
}

void fieldBroker::updateMassFlowRate(const std::shared_ptr<domain> domain)
{
    // interior
    updateMassFlowRateInterior_(domain);

#ifdef HAS_INTERFACE
    // Interfaces
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isConformalTreatment())
            continue;

        if (interf->isInternal())
        {
            updateMassFlowRateInterfaceSideField_(domain,
                                                  interf->masterInfoPtr());
            updateMassFlowRateInterfaceSideField_(domain,
                                                  interf->slaveInfoPtr());
        }
        else if (!interf->isFluidSolidType())
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(domain->index());

            updateMassFlowRateInterfaceSideField_(domain, interfaceSideInfoPtr);
        }
    }
#endif /* HAS_INTERFACE */

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        updateMassFlowRateBoundaryField_(
            domain, domain->zonePtr()->boundaryPtr(iBoundary));
    }
}

void fieldBroker::updateVelocityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().updateGradientField(domain->index());
    }
}

void fieldBroker::updatePressureGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTemperatureGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateSpecificEnthalpyGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateSpecificTotalEnthalpyGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().updateGradientField(domain->index());
    }
}

void fieldBroker::updateWallScaleGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTurbulentKineticEnergyGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTurbulentEddyFrequencyGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTurbulentDissipationRateGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTransitionOnsetReynoldsNumberGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTurbulentIntermittencyGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateDensityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateDisplacementGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateTotalDisplacementGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateMeshVelocityGradientField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Um_)
    {
        UmRef().updateGradientField(domain->index());
    }
}

void fieldBroker::updateVelocityBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updatePressureBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTemperatureBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateSpecificEnthalpyBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateSpecificTotalEnthalpyBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTurbulentKineticEnergyBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTurbulentEddyFrequencyBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTurbulentDissipationRateBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTransitionOnsetReynoldsNumberBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateTurbulentIntermittencyBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateDensityBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityBlendingFactorField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateVelocityPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updatePressurePrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTemperaturePrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateSpecificEnthalpyPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateSpecificTotalEnthalpyPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateWallScalePrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTurbulentKineticEnergyPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTurbulentEddyFrequencyPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTurbulentDissipationRatePrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTransitionOnsetReynoldsNumberPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTurbulentIntermittencyPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateDensityPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateDisplacementPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateTotalDisplacementPrevIterField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateVelocityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->U_)
    {
        URef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updatePressurePrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->p_)
    {
        pRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTemperaturePrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->T_)
    {
        TRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateSpecificEnthalpyPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h_)
    {
        hRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateSpecificTotalEnthalpyPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->h0_)
    {
        h0Ref().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateWallScalePrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->yScale_)
    {
        yScaleRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTurbulentKineticEnergyPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->k_)
    {
        kRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTurbulentEddyFrequencyPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->omega_)
    {
        omegaRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTurbulentDissipationRatePrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->epsilon_)
    {
        epsilonRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTransitionOnsetReynoldsNumberPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->ReTheta_)
    {
        ReThetaRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTurbulentIntermittencyPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->tRealm_->gamma_)
    {
        gammaRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateDensityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->rho_)
    {
        rhoRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->mu_)
    {
        muRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->cp_)
    {
        cpRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->lambda_)
    {
        lambdaRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateCompressibilityPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->psi_)
    {
        psiRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateDisplacementPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->D_)
    {
        DRef().updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateTotalDisplacementPrevTimeField(
    const std::shared_ptr<domain> domain)
{
    if (realmPtr_->Dt_)
    {
        DtRef().updatePrevTimeField(domain->index());
    }
}

// field API

// Public access

velocity& fieldBroker::URef()
{
    RETURN_REF_ALLOC(
        realmPtr_->U_, velocity, realmPtr_, realm::U_ID, n_states, high_res);
}

const velocity& fieldBroker::URef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->U_, velocity, realmPtr_, realm::U_ID, n_states, high_res);
}

simpleVectorField& fieldBroker::UrRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Ur_,
                         simpleVectorField,
                         realmPtr_,
                         realm::Ur_ID,
                         stk::topology::NODE_RANK);
};

const simpleVectorField& fieldBroker::UrRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Ur_,
                         simpleVectorField,
                         realmPtr_,
                         realm::Ur_ID,
                         stk::topology::NODE_RANK);
};

density& fieldBroker::rhoRef()
{
    RETURN_REF_ALLOC(realmPtr_->rho_,
                     density,
                     realmPtr_,
                     realm::rho_ID,
                     n_states,
                     high_res,
                     true);
}

const density& fieldBroker::rhoRef() const
{
    RETURN_REF_ALLOC(realmPtr_->rho_,
                     density,
                     realmPtr_,
                     realm::rho_ID,
                     n_states,
                     high_res,
                     true);
}

massFlowRate& fieldBroker::mDotRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->mDot_, massFlowRate, realmPtr_, realm::mDot_ID, n_states);
}

const massFlowRate& fieldBroker::mDotRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->mDot_, massFlowRate, realmPtr_, realm::mDot_ID, n_states);
}

nodeVectorField& fieldBroker::UmRef()
{
    RETURN_REF_ALLOC(realmPtr_->Um_,
                     nodeVectorField,
                     realmPtr_,
                     realm::Um_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

const nodeVectorField& fieldBroker::UmRef() const
{
    RETURN_REF_ALLOC(realmPtr_->Um_,
                     nodeVectorField,
                     realmPtr_,
                     realm::Um_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

nodeScalarField& fieldBroker::divUmRef()
{
    RETURN_REF_ALLOC(realmPtr_->divUm_,
                     nodeScalarField,
                     realmPtr_,
                     realm::divUm_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

const nodeScalarField& fieldBroker::divUmRef() const
{
    RETURN_REF_ALLOC(realmPtr_->divUm_,
                     nodeScalarField,
                     realmPtr_,
                     realm::divUm_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

// Protected access

pressure& fieldBroker::pRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->p_, pressure, realmPtr_, realm::p_ID, n_states, high_res);
}

const pressure& fieldBroker::pRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->p_, pressure, realmPtr_, realm::p_ID, n_states, high_res);
}

simpleScalarField& fieldBroker::p0Ref()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->p0_,
                         simpleScalarField,
                         realmPtr_,
                         realm::p0_ID,
                         stk::topology::NODE_RANK);
};

const simpleScalarField& fieldBroker::p0Ref() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->p0_,
                         simpleScalarField,
                         realmPtr_,
                         realm::p0_ID,
                         stk::topology::NODE_RANK);
};

simpleScalarField& fieldBroker::MaRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Ma_,
                         simpleScalarField,
                         realmPtr_,
                         realm::Ma_ID,
                         stk::topology::NODE_RANK);
};

const simpleScalarField& fieldBroker::MaRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Ma_,
                         simpleScalarField,
                         realmPtr_,
                         realm::Ma_ID,
                         stk::topology::NODE_RANK);
};

simpleScalarField& fieldBroker::CoRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Co_,
                         simpleScalarField,
                         realmPtr_,
                         realm::Co_ID,
                         stk::topology::ELEM_RANK);
};

const simpleScalarField& fieldBroker::CoRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->Co_,
                         simpleScalarField,
                         realmPtr_,
                         realm::Co_ID,
                         stk::topology::ELEM_RANK);
};

temperature& fieldBroker::TRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->T_, temperature, realmPtr_, realm::T_ID, n_states, high_res);
}

const temperature& fieldBroker::TRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->T_, temperature, realmPtr_, realm::T_ID, n_states, high_res);
}

simpleScalarField& fieldBroker::T0Ref()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->T0_,
                         simpleScalarField,
                         realmPtr_,
                         realm::T0_ID,
                         stk::topology::NODE_RANK);
}

const simpleScalarField& fieldBroker::T0Ref() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->T0_,
                         simpleScalarField,
                         realmPtr_,
                         realm::T0_ID,
                         stk::topology::NODE_RANK);
}

specificEnthalpy& fieldBroker::hRef()
{
    RETURN_REF_ALLOC(realmPtr_->h_,
                     specificEnthalpy,
                     realmPtr_,
                     realm::h_ID,
                     n_states,
                     high_res);
}

const specificEnthalpy& fieldBroker::hRef() const
{
    RETURN_REF_ALLOC(realmPtr_->h_,
                     specificEnthalpy,
                     realmPtr_,
                     realm::h_ID,
                     n_states,
                     high_res);
}

specificTotalEnthalpy& fieldBroker::h0Ref()
{
    RETURN_REF_ALLOC(realmPtr_->h0_,
                     specificTotalEnthalpy,
                     realmPtr_,
                     realm::h0_ID,
                     n_states,
                     high_res);
}

const specificTotalEnthalpy& fieldBroker::h0Ref() const
{
    RETURN_REF_ALLOC(realmPtr_->h0_,
                     specificTotalEnthalpy,
                     realmPtr_,
                     realm::h0_ID,
                     n_states,
                     high_res);
}

youngModulus& fieldBroker::ERef()
{
    RETURN_REF_ALLOC(realmPtr_->smRealm_->E_,
                     youngModulus,
                     realmPtr_,
                     smRealm::E_ID,
                     n_states,
                     false,
                     true);
}

const youngModulus& fieldBroker::ERef() const
{
    RETURN_REF_ALLOC(realmPtr_->smRealm_->E_,
                     youngModulus,
                     realmPtr_,
                     smRealm::E_ID,
                     n_states,
                     false,
                     true);
}

poissonRatio& fieldBroker::nuRef()
{
    RETURN_REF_ALLOC(realmPtr_->smRealm_->nu_,
                     poissonRatio,
                     realmPtr_,
                     smRealm::nu_ID,
                     n_states,
                     false,
                     true);
}

const poissonRatio& fieldBroker::nuRef() const
{
    RETURN_REF_ALLOC(realmPtr_->smRealm_->nu_,
                     poissonRatio,
                     realmPtr_,
                     smRealm::nu_ID,
                     n_states,
                     false,
                     true);
}

dynamicViscosity& fieldBroker::muRef()
{
    RETURN_REF_ALLOC(realmPtr_->mu_,
                     dynamicViscosity,
                     realmPtr_,
                     realm::mu_ID,
                     n_states,
                     false,
                     false);
}

const dynamicViscosity& fieldBroker::muRef() const
{
    RETURN_REF_ALLOC(realmPtr_->mu_,
                     dynamicViscosity,
                     realmPtr_,
                     realm::mu_ID,
                     n_states,
                     false,
                     false);
}

specificHeatCapacity& fieldBroker::cpRef()
{
    RETURN_REF_ALLOC(realmPtr_->cp_,
                     specificHeatCapacity,
                     realmPtr_,
                     realm::cp_ID,
                     n_states,
                     false,
                     true);
}

const specificHeatCapacity& fieldBroker::cpRef() const
{
    RETURN_REF_ALLOC(realmPtr_->cp_,
                     specificHeatCapacity,
                     realmPtr_,
                     realm::cp_ID,
                     n_states,
                     false,
                     true);
}

thermalConductivity& fieldBroker::lambdaRef()
{
    RETURN_REF_ALLOC(realmPtr_->lambda_,
                     thermalConductivity,
                     realmPtr_,
                     realm::lambda_ID,
                     n_states,
                     false,
                     false);
}

const thermalConductivity& fieldBroker::lambdaRef() const
{
    RETURN_REF_ALLOC(realmPtr_->lambda_,
                     thermalConductivity,
                     realmPtr_,
                     realm::lambda_ID,
                     n_states,
                     false,
                     false);
}

thermalExpansivity& fieldBroker::betaRef()
{
    RETURN_REF_ALLOC(realmPtr_->beta_,
                     thermalExpansivity,
                     realmPtr_,
                     realm::beta_ID,
                     n_states,
                     false,
                     false);
}

const thermalExpansivity& fieldBroker::betaRef() const
{
    RETURN_REF_ALLOC(realmPtr_->beta_,
                     thermalExpansivity,
                     realmPtr_,
                     realm::beta_ID,
                     n_states,
                     false,
                     false);
}

compressibility& fieldBroker::psiRef()
{
    RETURN_REF_ALLOC(realmPtr_->psi_,
                     compressibility,
                     realmPtr_,
                     realm::psi_ID,
                     n_states,
                     false,
                     false);
}

const compressibility& fieldBroker::psiRef() const
{
    RETURN_REF_ALLOC(realmPtr_->psi_,
                     compressibility,
                     realmPtr_,
                     realm::psi_ID,
                     n_states,
                     false,
                     false);
}

nodeScalarField& fieldBroker::yScaleRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->yScale_,
        nodeScalarField,
        realmPtr_,
        realm::yScale_ID,
        n_states,
        true,
        false,
        true,
        this->controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.strongDirichletWallScale_
            ? true
            : false);
}

const nodeScalarField& fieldBroker::yScaleRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->yScale_,
        nodeScalarField,
        realmPtr_,
        realm::yScale_ID,
        n_states,
        true,
        false,
        true,
        this->controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.strongDirichletWallScale_
            ? true
            : false);
}

turbulentKineticEnergy& fieldBroker::kRef()
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->k_,
                          turbulentKineticEnergy,
                          realmPtr_,
                          turbRealm::k_ID,
                          n_states,
                          high_res);
}

const turbulentKineticEnergy& fieldBroker::kRef() const
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->k_,
                          turbulentKineticEnergy,
                          realmPtr_,
                          turbRealm::k_ID,
                          n_states,
                          high_res);
}

turbulentEddyFrequency& fieldBroker::omegaRef()
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->omega_,
                          turbulentEddyFrequency,
                          realmPtr_,
                          turbRealm::omega_ID,
                          n_states,
                          high_res);
}

const turbulentEddyFrequency& fieldBroker::omegaRef() const
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->omega_,
                          turbulentEddyFrequency,
                          realmPtr_,
                          turbRealm::omega_ID,
                          n_states,
                          high_res);
}

turbulentDissipationRate& fieldBroker::epsilonRef()
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->epsilon_,
                          turbulentDissipationRate,
                          realmPtr_,
                          turbRealm::epsilon_ID,
                          n_states,
                          high_res);
}

const turbulentDissipationRate& fieldBroker::epsilonRef() const
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->epsilon_,
                          turbulentDissipationRate,
                          realmPtr_,
                          turbRealm::epsilon_ID,
                          n_states,
                          high_res);
}

transitionOnsetReynoldsNumber& fieldBroker::ReThetaRef()
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->ReTheta_,
                          transitionOnsetReynoldsNumber,
                          realmPtr_,
                          turbRealm::ReTheta_ID,
                          n_states,
                          high_res);
}

const transitionOnsetReynoldsNumber& fieldBroker::ReThetaRef() const
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->ReTheta_,
                          transitionOnsetReynoldsNumber,
                          realmPtr_,
                          turbRealm::ReTheta_ID,
                          n_states,
                          high_res);
}

turbulentIntermittency& fieldBroker::gammaRef()
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->gamma_,
                          turbulentIntermittency,
                          realmPtr_,
                          turbRealm::gamma_ID,
                          n_states,
                          high_res);
}

const turbulentIntermittency& fieldBroker::gammaRef() const
{
    RETURN_REF_ALLOC_TURB(realmPtr_->tRealm_->gamma_,
                          turbulentIntermittency,
                          realmPtr_,
                          turbRealm::gamma_ID,

                          n_states,
                          high_res);
}

heatFlowRate& fieldBroker::qDotRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->qDot_, heatFlowRate, realmPtr_, realm::qDot_ID, n_states);
}

const heatFlowRate& fieldBroker::qDotRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->qDot_, heatFlowRate, realmPtr_, realm::qDot_ID, n_states);
}

turbulentViscosity& fieldBroker::mutRef()
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->mut_,
                     turbulentViscosity,
                     realmPtr_,
                     turbRealm::mut_ID,
                     n_states,
                     false,
                     false);
}

const turbulentViscosity& fieldBroker::mutRef() const
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->mut_,
                     turbulentViscosity,
                     realmPtr_,
                     turbRealm::mut_ID,
                     n_states,
                     false,
                     false);
}

turbulentThermalConductivity& fieldBroker::lambdatRef()
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->lambdat_,
                     turbulentThermalConductivity,
                     realmPtr_,
                     turbRealm::lambdat_ID,
                     n_states,
                     false,
                     false);
}

const turbulentThermalConductivity& fieldBroker::lambdatRef() const
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->lambdat_,
                     turbulentThermalConductivity,
                     realmPtr_,
                     turbRealm::lambdat_ID,
                     n_states,
                     false,
                     false);
}

nodeScalarField& fieldBroker::muEffRef()
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->muEff_,
                     nodeScalarField,
                     realmPtr_,
                     turbRealm::muEff_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

const nodeScalarField& fieldBroker::muEffRef() const
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->muEff_,
                     nodeScalarField,
                     realmPtr_,
                     turbRealm::muEff_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

nodeScalarField& fieldBroker::lambdaEffRef()
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->lambdaEff_,
                     nodeScalarField,
                     realmPtr_,
                     turbRealm::lambdaEff_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

const nodeScalarField& fieldBroker::lambdaEffRef() const
{
    RETURN_REF_ALLOC(realmPtr_->tRealm_->lambdaEff_,
                     nodeScalarField,
                     realmPtr_,
                     turbRealm::lambdaEff_ID,
                     n_states,
                     false,
                     false,
                     false,
                     false);
}

simpleScalarField& fieldBroker::yMinRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->minimumDistanceToWall_,
                         simpleScalarField,
                         realmPtr_,
                         realm::yMin_ID,
                         stk::topology::NODE_RANK);
};

const simpleScalarField& fieldBroker::yMinRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->minimumDistanceToWall_,
                         simpleScalarField,
                         realmPtr_,
                         realm::yMin_ID,
                         stk::topology::NODE_RANK);
};

simpleScalarField& fieldBroker::F1Ref()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->tRealm_->fOneBlending_,
                         simpleScalarField,
                         realmPtr_,
                         turbRealm::F1_ID,
                         stk::topology::NODE_RANK);
};

const simpleScalarField& fieldBroker::F1Ref() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->tRealm_->fOneBlending_,
                         simpleScalarField,
                         realmPtr_,
                         turbRealm::F1_ID,
                         stk::topology::NODE_RANK);
};

simpleScalarField& fieldBroker::PkRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->tRealm_->tkeProduction_,
                         simpleScalarField,
                         realmPtr_,
                         turbRealm::Pk_ID,
                         stk::topology::NODE_RANK);
};

const simpleScalarField& fieldBroker::PkRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->tRealm_->tkeProduction_,
                         simpleScalarField,
                         realmPtr_,
                         turbRealm::Pk_ID,
                         stk::topology::NODE_RANK);
};

simpleScalarField& fieldBroker::uTauRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uTau_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uTau_ID);
}

const simpleScalarField& fieldBroker::uTauRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uTau_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uTau_ID);
}

simpleScalarField& fieldBroker::yPlusRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->yPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::yPlus_ID);
};

const simpleScalarField& fieldBroker::yPlusRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->yPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::yPlus_ID);
};

simpleScalarField& fieldBroker::uPlusRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uPlus_ID);
};

const simpleScalarField& fieldBroker::uPlusRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uPlus_ID);
};

simpleScalarField& fieldBroker::TPlusRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->TPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::TPlus_ID);
};

const simpleScalarField& fieldBroker::TPlusRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->TPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::TPlus_ID);
};

simpleScalarField& fieldBroker::yStarRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->yStar_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::yStar_ID);
};

const simpleScalarField& fieldBroker::yStarRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->yStar_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::yStar_ID);
};

simpleScalarField& fieldBroker::uStarRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uStar_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uStar_ID);
};

const simpleScalarField& fieldBroker::uStarRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uStar_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uStar_ID);
};

simpleScalarField& fieldBroker::duPlusdyPlusRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->duPlusdyPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::duPlusdyPlus_ID);
};

const simpleScalarField& fieldBroker::duPlusdyPlusRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->duPlusdyPlus_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::duPlusdyPlus_ID);
};

simpleScalarField& fieldBroker::uWallCoeffsRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uWallCoeffs_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uWallCoeffs_ID);
};

const simpleScalarField& fieldBroker::uWallCoeffsRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->uWallCoeffs_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::uWallCoeffs_ID);
};

simpleScalarField& fieldBroker::TWallCoeffsRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->TWallCoeffs_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::TWallCoeffs_ID);
};

const simpleScalarField& fieldBroker::TWallCoeffsRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->TWallCoeffs_,
                              simpleScalarField,
                              realmPtr_,
                              turbRealm::TWallCoeffs_ID);
};

simpleVectorField& fieldBroker::wallShearStressRef()
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->wallShearStress_,
                              simpleVectorField,
                              realmPtr_,
                              turbRealm::wallShearStress_ID);
};

const simpleVectorField& fieldBroker::wallShearStressRef() const
{
    RETURN_REF_ALLOC_AUX_WALL(realmPtr_->tRealm_->wallShearStress_,
                              simpleVectorField,
                              realmPtr_,
                              turbRealm::wallShearStress_ID);
};

nodeVectorField& fieldBroker::DtRef()
{
    RETURN_REF_ALLOC(realmPtr_->Dt_,
                     nodeVectorField,
                     realmPtr_,
                     realm::Dt_ID,
                     n_states == 2 ? n_states + 1 : n_states,
                     false,
                     false,
                     false,
                     false);
}

const nodeVectorField& fieldBroker::DtRef() const
{
    RETURN_REF_ALLOC(realmPtr_->Dt_,
                     nodeVectorField,
                     realmPtr_,
                     realm::Dt_ID,
                     n_states == 2 ? n_states + 1 : n_states,
                     false,
                     false,
                     false,
                     false);
}

displacement& fieldBroker::DRef()
{
    RETURN_REF_ALLOC(
        realmPtr_->D_, displacement, realmPtr_, realm::D_ID, n_states);
}

const displacement& fieldBroker::DRef() const
{
    RETURN_REF_ALLOC(
        realmPtr_->D_, displacement, realmPtr_, realm::D_ID, n_states);
}

simpleTensorField& fieldBroker::stressRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->smRealm_->sigma_,
                         simpleTensorField,
                         realmPtr_,
                         smRealm::sigma_ID,
                         stk::topology::NODE_RANK);
};

const simpleTensorField& fieldBroker::stressRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->smRealm_->sigma_,
                         simpleTensorField,
                         realmPtr_,
                         smRealm::sigma_ID,
                         stk::topology::NODE_RANK);
};

simpleTensorField& fieldBroker::strainRef()
{
    RETURN_REF_ALLOC_AUX(realmPtr_->smRealm_->epsilon_,
                         simpleTensorField,
                         realmPtr_,
                         smRealm::epsilon_ID,
                         stk::topology::NODE_RANK);
};

const simpleTensorField& fieldBroker::strainRef() const
{
    RETURN_REF_ALLOC_AUX(realmPtr_->smRealm_->epsilon_,
                         simpleTensorField,
                         realmPtr_,
                         smRealm::epsilon_ID,
                         stk::topology::NODE_RANK);
};

void fieldBroker::setupBoundaryConditions_(
    const std::shared_ptr<domain> domain,
    std::function<void(const accel::domain* domain,
                       const label iBoundary,
                       const boundaryPhysicalType bc_type,
                       const YAML::Node bc_values)> setupBC)
{
    const YAML::Node bc_list = domain->getYAMLBoundaryConditions();
    for (label iBoundary = 0; iBoundary < bc_list.size(); iBoundary++)
    {
        const YAML::Node patch = bc_list[iBoundary];
        if (!patch["type"])
        {
            errorMsg("fieldBroker: boundary patch " +
                     std::to_string(iBoundary) + " for domain `" +
                     domain->name() + "` does not define a `type`");
        }
        const boundaryPhysicalType boundary_type =
            convertBoundaryPhysicalTypeFromString(
                patch["type"].template as<std::string>());
        YAML::Node bc_values;
        if (patch["boundary_details"])
        {
            bc_values = patch["boundary_details"];
        }

        setupBC(domain.get(), iBoundary, boundary_type, bc_values);
    }
}

} /* namespace accel */
