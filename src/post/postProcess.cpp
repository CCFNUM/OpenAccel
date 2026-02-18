// File : postProcess.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "postProcess.h"

namespace accel
{

postProcess::postProcess(realm* realm, fs::path directory)
    : fieldBroker(realm), directory_(directory)
{
}

void postProcess::read(const YAML::Node& inputNode)
{
    const auto& sim = inputNode;

    if (sim["solver"])
    {
        const auto& solver = sim["solver"];

        if (solver["output_control"])
        {
            const auto& outputCtrl = solver["output_control"];

            if (outputCtrl["post_process"])
            {
                const auto& postProcessorBlockArray =
                    outputCtrl["post_process"];
                for (label iPostProcess = 0;
                     iPostProcess <
                     static_cast<label>(postProcessorBlockArray.size());
                     iPostProcess++)
                {
                    const auto& postProcessBlock =
                        postProcessorBlockArray[iPostProcess];

                    std::string name =
                        postProcessBlock["name"].template as<std::string>();
                    std::vector<std::string> location =
                        postProcessBlock["location"]
                            .template as<std::vector<std::string>>();
                    label frequency =
                        postProcessBlock["frequency"].template as<label>();
                    bool writeToFile =
                        postProcessBlock["write_to_file"].template as<bool>();
                    postProcessType type = convertPostProcessTypeFromString(
                        postProcessBlock["type"].template as<std::string>());

                    switch (type)
                    {
                        case postProcessType::reduction:
                            {
                                if (postProcessBlock["options"])
                                {
                                    const auto& optionsBlock =
                                        postProcessBlock["options"];

                                    reductionType rType =
                                        convertStatisticsTypeFromString(
                                            optionsBlock["type"]
                                                .template as<std::string>());
                                    std::string field =
                                        optionsBlock["field"]
                                            .template as<std::string>();

                                    std::unique_ptr<reductionObject> obj =
                                        std::make_unique<reductionObject>(
                                            this,
                                            name,
                                            type,
                                            location,
                                            frequency,
                                            writeToFile,
                                            rType,
                                            field);

                                    objectPtrArray_.push_back(std::move(obj));
                                }
                                else
                                {
                                    errorMsg("options block is not "
                                             "provided in post_process");
                                }
                            }
                            break;

                        case postProcessType::force:
                            {
                                bool calculateMoment = false;
                                std::array<scalar, SPATIAL_DIM> momentCenter = {
                                    0.0};
                                bool totalPrint = true;

                                // check if additional options are queried
                                if (postProcessBlock["options"])
                                {
                                    const auto& optionsBlock =
                                        postProcessBlock["options"];

                                    if (optionsBlock["calculate_moment"])
                                    {
                                        calculateMoment =
                                            optionsBlock["calculate_moment"]
                                                .template as<bool>();

                                        if (calculateMoment)
                                        {
                                            momentCenter =
                                                optionsBlock["moment_center"]
                                                    .template as<std::array<
                                                        scalar,
                                                        SPATIAL_DIM>>();
                                        }
                                    }

                                    if (optionsBlock["total_print"])
                                    {
                                        totalPrint = optionsBlock["total_print"]
                                                         .template as<bool>();
                                    }
                                }

                                std::unique_ptr<forceObject> obj =
                                    std::make_unique<forceObject>(
                                        this,
                                        name,
                                        type,
                                        location,
                                        frequency,
                                        writeToFile,
                                        calculateMoment,
                                        momentCenter,
                                        totalPrint);

                                objectPtrArray_.push_back(std::move(obj));
                            }
                            break;

                        case postProcessType::probe:
                            {
                                if (postProcessBlock["options"])
                                {
                                    const auto& optionsBlock =
                                        postProcessBlock["options"];

                                    std::array<scalar, SPATIAL_DIM>
                                        probeLocation =
                                            optionsBlock["probe_location"]
                                                .template as<
                                                    std::array<scalar,
                                                               SPATIAL_DIM>>();
                                    std::string field =
                                        optionsBlock["field"]
                                            .template as<std::string>();

                                    std::unique_ptr<probeObject> obj =
                                        std::make_unique<probeObject>(
                                            this,
                                            name,
                                            type,
                                            location,
                                            frequency,
                                            writeToFile,
                                            probeLocation,
                                            field);

                                    objectPtrArray_.push_back(std::move(obj));
                                }
                                else
                                {
                                    errorMsg("options block is not "
                                             "provided in post_process");
                                }
                            }
                            break;
                    }
                }
            }
        }
        else
        {
            errorMsg(
                "output_control block is not provided in the yaml input file");
        }
    }
    else
    {
        errorMsg("solver block is not provided in the yaml input file");
    }
}

void postProcess::update()
{
    if (controlsRef().isTransient())
    {
        for (auto& obj : objectPtrArray_)
        {
            if (controlsRef().getTimeStepCount() % obj->frequency_ == 0)
            {
                iter_ = controlsRef().getTimeStepCount();
                time_ = controlsRef().time;
                obj->update();
            }
        }
    }
    else
    {
        for (auto& obj : objectPtrArray_)
        {
            if (controlsRef().iter % obj->frequency_ == 0)
            {
                iter_ = controlsRef().iter;
                time_ = controlsRef().time; // dummy
                obj->update();
            }
        }
    }
}

std::string postProcess::instance()
{
    if (controlsRef().isTransient())
    {
        return std::to_string(time_);
    }
    else
    {
        return std::to_string(iter_);
    }
}

std::string postProcess::instanceHeader()
{
    if (controlsRef().isTransient())
    {
        return "time";
    }
    else
    {
        return "iter";
    }
}

postProcessObject::postProcessObject(postProcess* postProcessManagerPtr,
                                     std::string name,
                                     postProcessType type,
                                     std::vector<std::string> location,
                                     label frequency,
                                     bool writeToFile)
    : postProcessPtr_(postProcessManagerPtr), name_(name), type_(type),
      location_(location), frequency_(frequency), writeToFile_(writeToFile)
{
}

} // namespace accel
