// File       : controlsIO.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "controls.h"
#include "messager.h"

// std
#include <algorithm>
#include <cassert>
#include <list>
#include <sstream>

namespace accel
{

// IO

void controls::read(YAML::Node inputNode)
{
    if (messager::master())
    {
        std::cout << "Reading controls .." << std::endl;
    }

    const auto& sim = inputNode;
    if (sim["physical_analysis"])
    {
        const auto& flowAna = sim["physical_analysis"];

        if (flowAna["analysis_type"])
        {
            const auto& anaType = flowAna["analysis_type"];

            if (anaType["option"].template as<std::string>() == "transient")
            {
                analysisType_.transient_ = true;

                analysisType_.totalTime_ =
                    anaType["total_time"].template as<scalar>();
                analysisType_.initialTimestep_ = 0.0;

                if (anaType["time_steps"])
                {
                    const auto& dt_info = anaType["time_steps"];

                    if (dt_info["option"])
                    {
                        analysisType_.timeSteps_.mode_ =
                            convertTimestepModeFromString(
                                dt_info["option"].template as<std::string>());
                    }

                    switch (analysisType_.timeSteps_.mode_)
                    {
                        case timestepMode::constant:
                            {
                                analysisType_.initialTimestep_ =
                                    dt_info["timestep"].template as<scalar>();
                            }
                            break;

                        case timestepMode::specifiedInterval:
                        case timestepMode::periodicInterval:
                            {
                                analysisType_.initialTimestep_ =
                                    dt_info["timestep"].template as<scalar>();

                                if (dt_info["start_time"])
                                {
                                    analysisType_.timeSteps_.startTime_ =
                                        dt_info["start_time"]
                                            .template as<std::list<scalar>>();
                                }

                                if (dt_info["interval_length"])
                                {
                                    analysisType_.timeSteps_.intervalLength_ =
                                        dt_info["interval_length"]
                                            .template as<std::list<scalar>>();
                                }

                                if (dt_info["interval_timestep"])
                                {
                                    analysisType_.timeSteps_.timestepInterval_ =
                                        dt_info["interval_timestep"]
                                            .template as<std::list<scalar>>();
                                }

                                if (dt_info["period"])
                                {
                                    analysisType_.timeSteps_.period_ =
                                        dt_info["period"].template as<scalar>();
                                }
                            }
                            break;

                        case timestepMode::adaptive:
                            {
                                if (dt_info["initial_timestep"])
                                {
                                    analysisType_.initialTimestep_ =
                                        dt_info["initial_timestep"]
                                            .template as<scalar>();
                                }

                                if (dt_info["timestep_update_frequency"])
                                {
                                    const auto freq =
                                        dt_info["timestep_update_frequency"]
                                            .template as<label>();
                                    analysisType_.timeSteps_
                                        .timestepUpdateFrequency_ =
                                        std::max<label>(1, freq);
                                }

                                if (dt_info["timestep_adaptation"])
                                {
                                    const auto& adapt =
                                        dt_info["timestep_adaptation"];

                                    if (adapt["option"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_.option_ =
                                            convertTimestepAdaptationTypeFromString(
                                                adapt["option"]
                                                    .template as<
                                                        std::string>());
                                    }

                                    if (adapt["courant_number"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_
                                            .courantNumber_ =
                                            adapt["courant_number"]
                                                .template as<scalar>();
                                    }

                                    if (adapt["min_timestep"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_.minTimestep_ =
                                            adapt["min_timestep"]
                                                .template as<scalar>();
                                    }

                                    if (adapt["max_timestep"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_.maxTimestep_ =
                                            adapt["max_timestep"]
                                                .template as<scalar>();
                                    }

                                    if (adapt["timestep_decrease_factor"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_
                                            .timestepDecreaseFactor_ =
                                            adapt["timestep_decrease_factor"]
                                                .template as<scalar>();
                                    }

                                    if (adapt["timestep_increase_factor"])
                                    {
                                        analysisType_.timeSteps_
                                            .timestepAdaptation_
                                            .timestepIncreaseFactor_ =
                                            adapt["timestep_increase_factor"]
                                                .template as<scalar>();
                                    }
                                }
                            }
                            break;
                    }
                }
                else
                {
                    errorMsg("time_steps block is required for transient "
                             "analysis");
                }

                if (analysisType_.timeSteps_.mode_ == timestepMode::adaptive &&
                    analysisType_.initialTimestep_ <= 0.0)
                {
                    errorMsg("initial_timestep is required for adaptive "
                             "time stepping");
                }

                for (label i = 0; i < analysisTypeDictionary::DT_ENTRIES; i++)
                {
                    analysisType_.timestep_[i] = analysisType_.initialTimestep_;
                }
            }
            else if (anaType["option"].template as<std::string>() ==
                     "steady_state")
            {
                analysisType_.transient_ = false;
            }
            else
            {
                errorMsg("unrecognized analysis_type option");
            }
        }
        else
        {
            errorMsg(
                "analysis_type block is not provided in the yaml input file");
        }
    }
    else
    {
        errorMsg(
            "physical_analysis block is not provided in the yaml input file");
    }

    if (sim["solver"])
    {
        const auto& solver = sim["solver"];

        if (solver["solver_control"]["basic_settings"])
        {
            const auto& basicSettings =
                solver["solver_control"]["basic_settings"];

            if (basicSettings["advection_scheme"])
            {
                solver_.solverControl_.basicSettings_.advectionScheme_ =
                    convertAdvectionSchemeTypeFromString(
                        basicSettings["advection_scheme"]
                            .template as<std::string>());
            }

            if (basicSettings["turbulence_numerics"])
            {
                solver_.solverControl_.basicSettings_.turbulenceNumerics_ =
                    convertAdvectionSchemeTypeFromString(
                        basicSettings["turbulence_numerics"]
                            .template as<std::string>());
            }

            if (analysisType_.transient_)
            {
                if (basicSettings["transient_scheme"])
                {
                    solver_.solverControl_.basicSettings_.transientScheme_ =
                        convertTransientSchemeTypeFromString(
                            basicSettings["transient_scheme"]
                                .template as<std::string>());
                }
                else
                {
                    errorMsg("transient_scheme key is not provided in the yaml "
                             "input file");
                }
            }

            if (basicSettings["convergence_controls"])
            {
                const auto& convCtrl = basicSettings["convergence_controls"];

                solver_.solverControl_.basicSettings_.convergenceControl_
                    .minIterations_ =
                    convCtrl["min_iterations"].template as<label>();
                solver_.solverControl_.basicSettings_.convergenceControl_
                    .maxIterations_ =
                    convCtrl["max_iterations"].template as<label>();

                if (!analysisType_.transient_)
                {
                    solver_.solverControl_.basicSettings_.convergenceControl_
                        .physicalTimescale_ =
                        convCtrl["physical_timescale"].template as<scalar>();
                }

                if (convCtrl["relaxation_parameters"])
                {
                    const auto& rlxPars = convCtrl["relaxation_parameters"];

                    if (rlxPars["relax_mass"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .relaxMass_ =
                            rlxPars["relax_mass"].template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .relaxMass_ <= 1.0);
                    }
                    if (rlxPars["wall_scale_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .wallScaleRelaxationFactor_ =
                            rlxPars["wall_scale_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .wallScaleRelaxationFactor_ <= 1.0);
                    }
                    if (rlxPars["energy_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .energyRelaxationFactor_ =
                            rlxPars["energy_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .wallScaleRelaxationFactor_ <= 1.0);
                    }
                    if (rlxPars["velocity_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .velocityRelaxationFactor_ =
                            rlxPars["velocity_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .velocityRelaxationFactor_ <= 1.0);
                    }
                    if (rlxPars["pressure_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .pressureRelaxationFactor_ =
                            rlxPars["pressure_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .pressureRelaxationFactor_ <= 1.0);
                    }
                    if (rlxPars["turbulence_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .turbulenceRelaxationFactor_ =
                            rlxPars["turbulence_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .turbulenceRelaxationFactor_ <= 1.0);
                    }
                    if (rlxPars["solid_displacement_relaxation_factor"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceControl_.relaxationParameters_
                            .solidDisplacementRelaxationFactor_ =
                            rlxPars["solid_displacement_relaxation_factor"]
                                .template as<scalar>();

                        assert(solver_.solverControl_.basicSettings_
                                   .convergenceControl_.relaxationParameters_
                                   .solidDisplacementRelaxationFactor_ <= 1.0);
                    }
                }
            }
            else
            {
                errorMsg("convergence_controls block is not provided in the "
                         "yaml input file");
            }

            if (basicSettings["convergence_criteria"])
            {
                const auto& convCriteria =
                    basicSettings["convergence_criteria"];

                solver_.solverControl_.basicSettings_.convergenceCriteria_
                    .residualType_ = convertResidualTypeFromString(
                    convCriteria["residual_type"].template as<std::string>());
                solver_.solverControl_.basicSettings_.convergenceCriteria_
                    .residualTarget_ =
                    convCriteria["residual_target"].template as<scalar>();

                if (convCriteria["physics_convergence"])
                {
                    const auto& physConv = convCriteria["physics_convergence"];

                    if (physConv["enabled"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceCriteria_.physicsConvergence_.enabled_ =
                            physConv["enabled"].template as<bool>();
                    }

                    if (physConv["write_residuals"])
                    {
                        solver_.solverControl_.basicSettings_
                            .convergenceCriteria_.physicsConvergence_
                            .writeResiduals_ =
                            physConv["write_residuals"].template as<bool>();
                    }

                    if (physConv["criteria"])
                    {
                        for (const auto& item : physConv["criteria"])
                        {
                            solver_.solverControl_.basicSettings_
                                .convergenceCriteria_.physicsConvergence_
                                .criteria_.push_back(
                                    convertPhysicsConvergenceTypeFromString(
                                        item.template as<std::string>()));
                        }
                    }

                    if (physConv["targets"])
                    {
                        const auto& targets = physConv["targets"];
                        if (targets["fsi_interface_residual"])
                        {
                            solver_.solverControl_.basicSettings_
                                .convergenceCriteria_.physicsConvergence_
                                .fsiInterfaceResidualTarget_ =
                                targets["fsi_interface_residual"]
                                    .template as<scalar>();
                        }
                        if (targets["fsi_force_residual"])
                        {
                            solver_.solverControl_.basicSettings_
                                .convergenceCriteria_.physicsConvergence_
                                .fsiForceResidualTarget_ =
                                targets["fsi_force_residual"]
                                    .template as<scalar>();
                        }
                    }
                }
            }
            else
            {
                errorMsg("convergence_controls block is not provided in the "
                         "yaml input file");
            }

            if (basicSettings["interpolation_scheme"])
            {
                const auto& interpScheme =
                    basicSettings["interpolation_scheme"];

                if (interpScheme["velocity_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_.velocityInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["velocity_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["pressure_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_.pressureInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["pressure_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["temperature_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .temperatureInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["temperature_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["turbulent_kinetic_energy_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentKineticEnergyInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["turbulent_kinetic_energy_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme
                        ["turbulent_dissipation_rate_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentDissipationRateInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["turbulent_dissipation_rate_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["turbulent_eddy_frequency_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentEddyFrequencyInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["turbulent_eddy_frequency_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme
                        ["transition_onset_reynolds_number_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .transitionOnsetReynoldsNumberInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["transition_onset_reynolds_number_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["turbulent_intermittency_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentIntermittencyInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["turbulent_intermittency_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme["wall_scale_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_.wallScaleInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["wall_scale_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["volume_fraction_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .volumeFractionInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["volume_fraction_interpolation_type"]
                                .template as<std::string>());
                }

                // Gradient interpolation scheme types
                if (interpScheme["velocity_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .velocityGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["velocity_gradient_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["pressure_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .pressureGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["pressure_gradient_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["temperature_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .temperatureGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["temperature_gradient_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme["turbulent_kinetic_energy_gradient_"
                                 "interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentKineticEnergyGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["turbulent_kinetic_energy_gradient_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["turbulent_dissipation_rate_gradient_"
                                 "interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentDissipationRateGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["turbulent_dissipation_rate_gradient_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["turbulent_eddy_frequency_gradient_"
                                 "interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentEddyFrequencyGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["turbulent_eddy_frequency_gradient_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["transition_onset_reynolds_number_gradient_"
                                 "interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .transitionOnsetReynoldsNumberGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["transition_onset_reynolds_number_"
                                         "gradient_interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme
                        ["turbulent_intermittency_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .turbulentIntermittencyGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme["turbulent_intermittency_gradient_"
                                         "interpolation_type"]
                                .template as<std::string>());
                }

                if (interpScheme["wall_scale_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .wallScaleGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["wall_scale_gradient_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme["volume_fraction_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .volumeFractionGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["volume_fraction_gradient_interpolation_type"]
                                    .template as<std::string>());
                }

                if (interpScheme["displacement_gradient_interpolation_type"])
                {
                    solver_.solverControl_.basicSettings_
                        .interpolationSchemeType_
                        .displacementGradientInterpolationType_ =
                        convertInterpolationSchemeTypeFromString(
                            interpScheme
                                ["displacement_gradient_interpolation_type"]
                                    .template as<std::string>());
                }
            }
        }
        else
        {
            errorMsg(
                "basic settings block is not provided in the yaml input file");
        }

        // optional block
        if (solver["solver_control"]["advanced_options"])
        {
            const auto& advancedOptions =
                solver["solver_control"]["advanced_options"];

            if (advancedOptions["pressure_level_information"])
            {
                const auto& pressureLevelInfo =
                    advancedOptions["pressure_level_information"];

                if (pressureLevelInfo["option"])
                {
                    solver_.solverControl_.advancedOptions_
                        .pressureLevelInformation_.option_ =
                        convertPressureLevelInformationSpecificationFromString(
                            pressureLevelInfo["option"]
                                .template as<std::string>());

                    if (solver_.solverControl_.advancedOptions_
                            .pressureLevelInformation_.option_ ==
                        pressureLevelInformationSpecification::
                            cartesianCoordinates)
                    {
                        solver_.solverControl_.advancedOptions_
                            .pressureLevelInformation_.cartesianCoordinates_ =
                            pressureLevelInfo["cartesian_coordinates"]
                                .template as<std::vector<scalar>>();
                    }
                }

                if (pressureLevelInfo["relative_pressure_level"])
                {
                    solver_.solverControl_.advancedOptions_
                        .pressureLevelInformation_.relativePressureLevel_ =
                        pressureLevelInfo["relative_pressure_level"]
                            .template as<scalar>();
                }
            }

#ifdef HAS_INTERFACE
            if (advancedOptions["interface_transfer"])
            {
                const auto& interfaceTransfer =
                    advancedOptions["interface_transfer"];

                if (interfaceTransfer["search_tolerance"])
                {
                    solver_.solverControl_.advancedOptions_.interfaceTransfer_
                        .searchTolerance_ =
                        interfaceTransfer["search_tolerance"]
                            .template as<scalar>();
                }

                if (interfaceTransfer["search_expansion_factor"])
                {
                    solver_.solverControl_.advancedOptions_.interfaceTransfer_
                        .searchExpansionFactor_ =
                        interfaceTransfer["search_expansion_factor"]
                            .template as<scalar>();
                }

                if (interfaceTransfer["force_research"])
                {
                    solver_.solverControl_.advancedOptions_.interfaceTransfer_
                        .forceResearch_ =
                        interfaceTransfer["force_research"].template as<bool>();
                }

                if (interfaceTransfer["verbose"])
                {
                    solver_.solverControl_.advancedOptions_.interfaceTransfer_
                        .verbose_ =
                        interfaceTransfer["verbose"].template as<label>();
                }

                if (interfaceTransfer["conservative_flux_transfer"])
                {
                    solver_.solverControl_.advancedOptions_.interfaceTransfer_
                        .conservativeFluxTransfer_ =
                        interfaceTransfer["conservative_flux_transfer"]
                            .template as<bool>();
                }
            }
#endif /* HAS_INTERFACE */

            if (advancedOptions["equation_controls"])
            {
                if (advancedOptions["equation_controls"]["sub_iterations"])
                {
                    if (advancedOptions["equation_controls"]["sub_iterations"]
                                       ["pressure_correction"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.subIterations_
                            .pressureCorrection_ =
                            advancedOptions["equation_controls"]
                                           ["sub_iterations"]
                                           ["pressure_correction"]
                                               .template as<label>();
                    }

                    if (advancedOptions["equation_controls"]["sub_iterations"]
                                       ["solid_displacement"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.subIterations_
                            .solidDisplacement_ =
                            advancedOptions["equation_controls"]
                                           ["sub_iterations"]
                                           ["solid_displacement"]
                                               .template as<label>();
                    }

                    if (advancedOptions["equation_controls"]["sub_iterations"]
                                       ["segregated_flow"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.subIterations_.segregatedFlow_ =
                            advancedOptions["equation_controls"]
                                           ["sub_iterations"]["segregated_flow"]
                                               .template as<label>();
                    }

                    if (advancedOptions["equation_controls"]["sub_iterations"]
                                       ["volume_fraction"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.subIterations_.volumeFraction_ =
                            advancedOptions["equation_controls"]
                                           ["sub_iterations"]["volume_fraction"]
                                               .template as<label>();
                    }
                }

                if (advancedOptions["equation_controls"]
                                   ["volume_fraction_smoothing"])
                {
                    if (advancedOptions["equation_controls"]
                                       ["volume_fraction_smoothing"]
                                       ["smooth_volume_fraction"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.volumeFractionSmoothing_
                            .smoothVolumeFraction_ =
                            advancedOptions["equation_controls"]
                                           ["volume_fraction_smoothing"]
                                           ["smooth_volume_fraction"]
                                               .template as<bool>();
                    }

                    if (advancedOptions["equation_controls"]
                                       ["volume_fraction_smoothing"]
                                       ["smoothing_iterations"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.volumeFractionSmoothing_
                            .smoothingIterations_ =
                            advancedOptions["equation_controls"]
                                           ["volume_fraction_smoothing"]
                                           ["smoothing_iterations"]
                                               .template as<label>();
                    }

                    if (advancedOptions["equation_controls"]
                                       ["volume_fraction_smoothing"]
                                       ["fourier_number"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.volumeFractionSmoothing_
                            .fourierNumber_ =
                            advancedOptions["equation_controls"]
                                           ["volume_fraction_smoothing"]
                                           ["fourier_number"]
                                               .template as<scalar>();
                    }
                }

                if (advancedOptions["equation_controls"]["mesh_motion"])
                {
                    if (advancedOptions["equation_controls"]["mesh_motion"]
                                       ["freeze_per_timestep"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.meshMotion_.freezePerTimestep_ =
                            advancedOptions["equation_controls"]["mesh_motion"]
                                           ["freeze_per_timestep"]
                                               .template as<bool>();
                    }

                    if (advancedOptions["equation_controls"]["mesh_motion"]
                                       ["max_smoothing_iters"])
                    {
                        solver_.solverControl_.advancedOptions_
                            .equationControls_.meshMotion_
                            .maxSmoothingIterations_ =
                            advancedOptions["equation_controls"]["mesh_motion"]
                                           ["max_smoothing_iters"]
                                               .template as<label>();
                    }
                }
            }
        }

        if (solver["solver_control"]["expert_parameters"])
        {
            const auto& expertParameters =
                solver["solver_control"]["expert_parameters"];

            if (expertParameters["disable_momentum_predictor"])
            {
                solver_.solverControl_.expertParameters_
                    .disableMomentumPredictor_ =
                    expertParameters["disable_momentum_predictor"]
                        .template as<bool>();
            }

            if (expertParameters["limit_gradients"])
            {
                solver_.solverControl_.expertParameters_.limitGradients_ =
                    expertParameters["limit_gradients"].template as<bool>();
            }

            if (expertParameters["correct_gradients"])
            {
                solver_.solverControl_.expertParameters_.correctGradients_ =
                    expertParameters["correct_gradients"].template as<bool>();
            }

            if (expertParameters["incremental_gradient_change"])
            {
                solver_.solverControl_.expertParameters_
                    .incrementalGradientChange_ =
                    expertParameters["incremental_gradient_change"]
                        .template as<bool>();
            }

            if (expertParameters["false_mass_accumulation"])
            {
                solver_.solverControl_.expertParameters_
                    .falseMassAccumulation_ =
                    expertParameters["false_mass_accumulation"]
                        .template as<bool>();
            }

            if (expertParameters["consistent"])
            {
                solver_.solverControl_.expertParameters_.consistent_ =
                    expertParameters["consistent"].template as<bool>();
            }

            if (expertParameters["fractional_step_method"])
            {
                solver_.solverControl_.expertParameters_.fractionalStepMethod_ =
                    expertParameters["fractional_step_method"]
                        .template as<bool>();
            }

            if (expertParameters["coriolis_production_turbulence"])
            {
                solver_.solverControl_.expertParameters_
                    .coriolisProductionTurbulence_ =
                    expertParameters["coriolis_production_turbulence"]
                        .template as<bool>();
            }

            // ensure that consistent and fractional step are not both enabled
            if (solver_.solverControl_.expertParameters_.consistent_ &&
                solver_.solverControl_.expertParameters_.fractionalStepMethod_)
            {
                errorMsg("consistent or fractional step method must be "
                         "enabled, not both");
            }

            if (expertParameters["body_force_redistribution"])
            {
                solver_.solverControl_.expertParameters_
                    .bodyForceRedistribution_ =
                    expertParameters["body_force_redistribution"]
                        .template as<bool>();
            }

            if (expertParameters["geometric_wall_distance_calculation"])
            {
                solver_.solverControl_.expertParameters_
                    .geometricWallDistanceCalculation_ =
                    expertParameters["geometric_wall_distance_calculation"]
                        .template as<bool>();
            }

            if (expertParameters["strong_dirichlet_wall_scale"])
            {
                solver_.solverControl_.expertParameters_
                    .strongDirichletWallScale_ =
                    expertParameters["strong_dirichlet_wall_scale"]
                        .template as<bool>();
            }

            if (expertParameters["volume_fraction_compressive_beta_max"])
            {
                solver_.solverControl_.expertParameters_
                    .volumeFractionBlendingFactorMax_ =
                    expertParameters["volume_fraction_compressive_beta_max"]
                        .template as<scalar>();
            }

            if (expertParameters["bandwidth_reduction"])
            {
                solver_.solverControl_.expertParameters_.bandwidthReduction_ =
                    expertParameters["bandwidth_reduction"].template as<bool>();
            }
        }

        if (solver["output_control"])
        {
            const auto& outputCtrl = solver["output_control"];

            solver_.outputControl_.filePath_ =
                outputCtrl["file_path"].template as<std::string>();
            if (outputCtrl["output_frequency"])
            {
                const auto& outputFreq = outputCtrl["output_frequency"];
                if (analysisType_.transient_)
                {
                    if (!outputFreq.IsMap())
                    {
                        errorMsg("output_frequency must be a map for "
                                 "transient simulations");
                    }

                    if (outputFreq["option"])
                    {
                        solver_.outputControl_.outputFrequency_
                            .option_ = convertOutputFrequencyTypeFromString(
                            outputFreq["option"].template as<std::string>());
                    }

                    if (outputFreq["time_interval"])
                    {
                        solver_.outputControl_.outputFrequency_.timeInterval_ =
                            outputFreq["time_interval"].template as<scalar>();
                    }

                    if (outputFreq["timestep_interval"])
                    {
                        solver_.outputControl_.outputFrequency_
                            .timestepInterval_ = outputFreq["timestep_interval"]
                                                     .template as<label>();
                    }
                }
                else
                {
                    if (outputFreq.IsMap())
                    {
                        errorMsg("output_frequency must be an integer for "
                                 "steady-state simulations");
                    }
                    solver_.outputControl_.outputFrequency_.option_ =
                        outputFrequencyType::timestepInterval;
                    solver_.outputControl_.outputFrequency_.timestepInterval_ =
                        outputFreq.template as<label>();
                }
            }
            solver_.outputControl_.outputFields_ =
                outputCtrl["output_fields"]
                    .template as<std::vector<std::string>>();
            if (outputCtrl["corrected_boundary_values"])
            {
                solver_.outputControl_.correctedBoundaryValues_ =
                    outputCtrl["corrected_boundary_values"].template as<bool>();
            }

            if (outputCtrl["restart_file_name"])
            {
                solver_.outputControl_.restartFileName_ =
                    outputCtrl["restart_file_name"].template as<std::string>();
            }
            if (outputCtrl["restart_frequency"])
            {
                solver_.outputControl_.restartFrequency_ =
                    outputCtrl["restart_frequency"].template as<label>();
            }
            if (outputCtrl["match_final_time"] &&
                outputCtrl["match_final_time"].template as<bool>())
            {
                solver_.outputControl_.matchFinalTime_ = true;
            }
            if (outputCtrl["write_timestep_info"])
            {
                solver_.outputControl_.writeTimestepInfo_ =
                    outputCtrl["write_timestep_info"].template as<bool>();
            }
        }
        else
        {
            errorMsg(
                "output_control block is not provided in the yaml input file");
        }

        if (solver["restart_control"])
        {
            const auto& restartCtrl = solver["restart_control"];
            solver_.restartControl_.isRestart_ = true;

            if (restartCtrl["file_path"])
            {
                solver_.restartControl_.inputFilePath_ =
                    restartCtrl["file_path"].template as<std::string>();
            }
            else
            {
                errorMsg("file_path for restart file in restart_control is "
                         "required");
            }

            if (restartCtrl["write_initial"] &&
                restartCtrl["write_initial"].template as<bool>())
            {
                solver_.restartControl_.writeInitial_ = true;
            }

            // if restart time is not specified, the last stored time in the
            // restart file will be used, otherwise STK will interpolate the
            // data to the specified restart time.
            if (restartCtrl["time"])
            {
                solver_.restartControl_.restartTime_ =
                    restartCtrl["time"].template as<scalar>();
            }

            if (restartCtrl["interpolation_type"])
            {
                const timeInterpolationSchemeType interp_type =
                    convertTimeInterpolationSchemeTypeFromString(
                        restartCtrl["interpolation_type"]
                            .template as<std::string>());

                if (interp_type == timeInterpolationSchemeType::piecewiseLinear)
                {
                    // FIXME: See related FIXME in
                    // nodeField<N, M>::initializeField()
                    errorMsg("interpolation_type piecewise_linear in "
                             "restart_control is currently not supported");
                    solver_.restartControl_.timeMatchOption_ =
                        stk::io::MeshField::LINEAR_INTERPOLATION;
                }
                else if (interp_type == timeInterpolationSchemeType::closest)
                {
                    solver_.restartControl_.timeMatchOption_ =
                        stk::io::MeshField::CLOSEST;
                }
                else
                {
                    errorMsg(
                        "interpolation_type in restart_control is unknown");
                }
            }

            if (restartCtrl["keep_snapshots"])
            {
                solver_.restartControl_.keepNRestartSnapshots_ =
                    restartCtrl["keep_snapshots"].template as<label>();
            }
        }

        // append results if results file exists and this is a restart
        std::ostringstream partition;
        if (messager::parallel())
        {
            partition << '.' << messager::nProcs() << '.'
                      << messager::myProcNo();
        }
        if (solver_.restartControl_.isRestart_ &&
            fs::exists(
                fs::path(solver_.outputControl_.filePath_ + partition.str())))
        {
            solver_.outputControl_.writeMode_ = stk::io::APPEND_RESULTS;
        }
    }
    else
    {
        errorMsg("solver block is not provided in the yaml input file");
    }

    if (messager::master())
    {
        std::cout << "Finished reading controls .." << std::endl << std::endl;
    }
}

} // namespace accel
