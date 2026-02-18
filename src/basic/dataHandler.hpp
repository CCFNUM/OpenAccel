// File : dataHandler.hpp
// Created : Tue Apr 30 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <size_t N>
void inputData<N>::query(const YAML::Node node,
                         std::string inputDataKey,
                         std::string inputTypeKey)
{
    // read and store input type
    std::string inputTypeString = "constant"; // default
    if (node[inputTypeKey])
    {
        inputTypeString = node[inputTypeKey].template as<std::string>();
    }
    ::accel::tolower(inputTypeString);
    type_ = convertInputTypeFromString(inputTypeString);

    // setup
    switch (type_)
    {
        case inputDataType::constant:
            queryConstantValue(node, inputDataKey);
            break;

        case inputDataType::expression:
            queryExpression(node, inputDataKey);
            break;

        case inputDataType::timeTable:
            queryTimeTable(node, "file_path");
            break;

        case inputDataType::profileData:
            queryProfileData(node, "file_path");
            break;

        default:
            errorMsg("input type `" + toString(type_) + "` not implemented");
            break;
    };
}

template <size_t N>
void inputData<N>::queryMulti(const YAML::Node node,
                              std::vector<std::string> inputDataKeys,
                              std::string inputTypeKey)
{
    // read and store input type
    std::string inputTypeString = "constant"; // default
    if (node[inputTypeKey])
    {
        inputTypeString = node[inputTypeKey].template as<std::string>();
    }
    ::accel::tolower(inputTypeString);
    type_ = convertInputTypeFromString(inputTypeString);

    // setup
    switch (type_)
    {
        case inputDataType::constant:
            {
                // create tmp node to fill with all components
                YAML::Node seq(YAML::NodeType::Sequence);
                for (label i = 0; i < N; i++)
                {
                    seq.push_back(node[inputDataKeys[i]]);
                }

                YAML::Node consolidatedNode;
                consolidatedNode["value"] = seq;

                queryConstantValue(consolidatedNode, "value");
            }
            break;

        case inputDataType::expression:
            {
                // create tmp node to fill with all components
                YAML::Node seq(YAML::NodeType::Sequence);
                for (label i = 0; i < N; i++)
                {
                    seq.push_back(node[inputDataKeys[i]]);
                }

                YAML::Node consolidatedNode;
                consolidatedNode["value"] = seq;

                queryExpression(consolidatedNode, "value");
            }
            break;

        case inputDataType::timeTable:
            queryTimeTable(node, "file_path");
            break;

        case inputDataType::profileData:
            queryProfileData(node, "file_path");
            break;

        default:
            errorMsg("input type `" + toString(type_) + "` not implemented");
            break;
    };
}

template <size_t N>
void inputData<N>::queryConstantValue(const YAML::Node node,
                                      std::string inputDataKey,
                                      const std::string internalMethod)
{
    if (!node[inputDataKey])
    {
        errorMsg("constant input-data: entry for key " + inputDataKey +
                 " is not provided");
    }

    try
    {
        if (node[inputDataKey].template as<std::string>() == internalMethod)
        {
            // Internal method specifies a default treatment for constant values
            // that are location dependent, where location may be a function of
            // time. An example would be the `flow_direction` vector computed
            // based on the surface normal on curved surfaces or moving meshes.
            useInternalMethod_ = true;
            rawStringValue_ = internalMethod;
            return;
        }
    }
    catch (...)
    {
        // nop
    }

    if (N == 1)
    {
        value_[0] = node[inputDataKey].template as<scalar>();
    }
    else
    {
        auto val = node[inputDataKey].template as<std::vector<scalar>>();

        for (label i = 0; i < N; i++)
        {
            value_[i] = val[i];
        }
    }
}

template <size_t N>
void inputData<N>::queryExpression(const YAML::Node node,
                                   std::string inputDataKey)
{
    if (!node[inputDataKey])
    {
        errorMsg("expression input-data: entry for key " + inputDataKey +
                 " is not provided");
    }

    if (N == 1)
    {
        expression_[0] = node[inputDataKey].template as<std::string>();
    }
    else
    {
        auto val = node[inputDataKey].template as<std::vector<std::string>>();

        for (label i = 0; i < N; i++)
        {
            expression_[i] = val[i];
        }
    }
}

template <size_t N>
void inputData<N>::queryTimeTable(const YAML::Node node,
                                  std::string filePathKey)
{
    if (!node[filePathKey])
    {
        errorMsg("`time_table` set but `" + filePathKey + "` field is missing");
    }

    std::string filePath = node[filePathKey].template as<std::string>();
    std::string hdf5_group = "dataset";
    if (node && node["hdf5_group"])
    {
        hdf5_group = node["hdf5_group"].template as<std::string>();
    }

    H5IO io;
    io.open_file(filePath);
    io.read_dataset(hdf5_group, timeTable_);
    io.close_file();

    this->setupInterpolator_(&node);
}

template <size_t N>
void inputData<N>::queryProfileData(const YAML::Node node,
                                    std::string filePathKey)
{
    if (!node[filePathKey])
    {
        errorMsg("`profile_data` set but `" + filePathKey +
                 "` field is missing");
    }

    std::string filePath = node[filePathKey].template as<std::string>();
    std::string coordsDataset =
        node["coords_dataset"].template as<std::string>();
    std::string fieldDataset = node["field_dataset"].template as<std::string>();

    H5IO io;
    io.open_file(filePath);
    io.read_dataset(coordsDataset, scatterPoints_);
    io.read_dataset(fieldDataset, scatterValues_);
    io.close_file();

    // other parameters
    if (node["interpolation_options"])
    {
        const auto& interpolationOptionsNode = node["interpolation_options"];

        donorCount_ =
            interpolationOptionsNode["donor_points_count"].template as<label>();
        distancePower_ = interpolationOptionsNode["distance_power_parameter"]
                             .template as<scalar>();
    }
}

} // namespace accel
