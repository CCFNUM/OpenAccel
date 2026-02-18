// File : dataHandler.h
// Created : Tue Apr 30 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Templated input data storage with interpolation and time table
// support
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DATAHANDLER_H
#define DATAHANDLER_H

// external
#include "tabular_props/BSpline.h"
#include "tabular_props/H5IO.h"

// code
#include "types.h"

// return flat index of a table of data structured as follows (considering
// column size is m): t0, t1, ..., tm-1, x0, x1, ..., xm-1, y0, y1, ..., ym-1
// given:
// RIDX: row index
// CIDX: col index
// NROWS: column size
#define FIDX(RIDX, CIDX, NROWS) (CIDX) * (NROWS) + (RIDX)

namespace accel
{

// Bring external tabular_props types into accel namespace
using sierra::nalu::BSpline1D;
using sierra::nalu::H5IO;

class baseInputData
{
public:
    baseInputData() = default;
    virtual ~baseInputData() = default;
};

template <size_t N>
class inputData : public baseInputData
{
public:
    static constexpr size_t NComponents = N;

    inputData()
        : type_(inputDataType::null),
          interpType_(timeInterpolationSchemeType::BSpline), interpOrder_(2),
          value_{0}, useInternalMethod_(false)
    {
    }

    // Access
    bool useInternalMethod() const
    {
        return useInternalMethod_;
    }

    timeInterpolationSchemeType interpolationType() const
    {
        return interpType_;
    }

    label interpolationOrder() const
    {
        return interpOrder_;
    }

    const inputDataType& type() const
    {
        return type_;
    }

    scalar* value()
    {
        assert(!this->useInternalMethod_);
        return value_.data();
    };

    const scalar* value() const
    {
        assert(!this->useInternalMethod_);
        return value_.data();
    };

    std::string* expression()
    {
        assert(!expression_.empty());
        return expression_.data();
    };

    const std::string* expression() const
    {
        assert(!expression_.empty());
        return expression_.data();
    };

    std::vector<scalar>& timeTable()
    {
        return timeTable_;
    };

    const std::vector<scalar>& timeTable() const
    {
        return timeTable_;
    };

    std::vector<scalar>& scatterPoints()
    {
        return scatterPoints_;
    };

    const std::vector<scalar>& scatterPoints() const
    {
        return scatterPoints_;
    };

    std::vector<scalar>& scatterValues()
    {
        return scatterValues_;
    };

    const std::vector<scalar>& scatterValues() const
    {
        return scatterValues_;
    };

    scalar distancePower() const
    {
        return distancePower_;
    }

    label donorCount() const
    {
        return donorCount_;
    }

    label rawLabelValue() const
    {
        return rawLabelValue_;
    }

    scalar rawScalarValue() const
    {
        assert(!this->useInternalMethod_);
        return rawScalarValue_;
    }

    std::string rawStringValue() const
    {
        return rawStringValue_;
    }

    bool rawBoolValue() const
    {
        return rawBoolValue_;
    }

    // Operations

    void query(const YAML::Node node,
               std::string inputDataKey = "value",
               std::string inputTypeKey = "input_type");

    void queryMulti(const YAML::Node node,
                    std::vector<std::string> inputDataKeys,
                    std::string inputTypeKey = "input_type");

    void queryConstantValue(const YAML::Node node,
                            std::string inputDataKey,
                            const std::string internalMethod = "");

    void queryExpression(const YAML::Node node, std::string inputDataKey);

    void queryTimeTable(const YAML::Node node, std::string filePathKey);

    void queryProfileData(const YAML::Node node, std::string filePathKey);

    std::array<scalar, N> interpolate(const scalar t0)
    {
        assert(!timeTable_.empty());
        assert(timeTable_.size() % (N + 1) == 0);

        const label columnSize = timeTable_.size() / (N + 1);

        if (interpType_ == timeInterpolationSchemeType::BSpline)
        {
            assert(!splines_.empty());
            label i = 0;
            for (const auto& spline : splines_)
            {
                value_[i++] = spline->value(t0);
            }
        }
        else if (interpType_ == timeInterpolationSchemeType::piecewiseLinear)
        {
            if (t0 < timeTable_[FIDX(0, 0, columnSize)]) // Lower bound
            {
                for (label i = 0; i < N; i++)
                {
                    value_[i] = timeTable_[FIDX(0, i + 1, columnSize)];
                }
            }
            else if (t0 >= timeTable_[FIDX(
                               columnSize - 1, 0, columnSize)]) // Upper bound
            {
                for (label i = 0; i < N; i++)
                {
                    value_[i] =
                        timeTable_[FIDX(columnSize - 1, i + 1, columnSize)];
                }
            }
            else
            {
                for (label i = 0; i < columnSize - 1; i++)
                {
                    if (t0 >= timeTable_[FIDX(i, 0, columnSize)] &&
                        t0 < timeTable_[FIDX(i + 1, 0, columnSize)])
                    {
                        const scalar t1 = timeTable_[FIDX(i, 0, columnSize)];
                        const scalar t2 =
                            timeTable_[FIDX(i + 1, 0, columnSize)];

                        for (label j = 0; j < N; j++)
                        {
                            const scalar y1 =
                                timeTable_[FIDX(i, j + 1, columnSize)];
                            const scalar y2 =
                                timeTable_[FIDX(i + 1, j + 1, columnSize)];

                            const scalar slope = (y2 - y1) / (t2 - t1);
                            const scalar y_intercept = y1 - (slope * t1);

                            value_[j] = (slope * t0) + y_intercept;
                        }

                        break;
                    }
                }
            }
        }
        // TODO: Implement closest time interpolation scheme
        // timeInterpolationSchemeType::closest
        else
        {
            errorMsg(
                "time interpolation method for time table is not implemented");
        }

        return value_;
    }

    // Setters

    void setType(inputDataType t)
    {
        type_ = t;
    };

    void setValue(scalar* valuePtr)
    {
        assert(!this->useInternalMethod_);
        for (label i = 0; i < N; i++)
        {
            value_[i] = valuePtr[i];
        }
    }

    void setValue(std::vector<scalar> value)
    {
        assert(!this->useInternalMethod_);
        assert(N == value.size());

        label i = 0;
        for (auto cv : value)
        {
            value_[i++] = cv;
        }
    }

    void setValue(std::initializer_list<scalar> value)
    {
        assert(!this->useInternalMethod_);
        assert(N == value.size());

        label i = 0;
        for (auto cv : value)
        {
            value_[i++] = cv;
        }
    }

    void setExpression(std::string* expressionPtr)
    {
        assert(N == expression_.size());
        for (label i = 0; i < N; i++)
        {
            expression_[i] = expressionPtr[i];
        }
    }

    void setExpression(std::vector<std::string> expression)
    {
        assert(N == expression.size());
        assert(N == expression_.size());

        label i = 0;
        for (auto cv : expression)
        {
            expression_[i++] = cv;
        }
    }

    void setExpression(std::initializer_list<std::string> expression)
    {
        assert(N == expression.size());
        assert(N == expression_.size());

        label i = 0;
        for (auto cv : expression)
        {
            expression_[i++] = cv;
        }
    }

    void setTimeTable(std::vector<scalar>& timeTable)
    {
        timeTable_ = timeTable;
        this->setupInterpolator_();
    }

    void setTimeTable(std::initializer_list<scalar> timeTable)
    {
        timeTable_.resize(timeTable.size());

        label i = 0;
        for (auto cv : timeTable)
        {
            timeTable_[i++] = cv;
        }
        this->setupInterpolator_();
    }

    void setScatterPoints(std::vector<scalar>& scatterPoints)
    {
        scatterPoints_ = scatterPoints;
    }

    void setScatterPoints(std::initializer_list<scalar> scatterPoints)
    {
        scatterPoints_.resize(scatterPoints.size());

        label i = 0;
        for (auto cv : scatterPoints)
        {
            scatterPoints_[i++] = cv;
        }
    }

    void setScatterValues(std::vector<scalar>& scatterValues)
    {
        scatterValues_ = scatterValues;
    }

    void setScatterValues(std::initializer_list<scalar> scatterValues)
    {
        scatterValues_.resize(scatterValues.size());

        label i = 0;
        for (auto cv : scatterValues)
        {
            scatterValues_[i++] = cv;
        }
    }

    void setProfileData(std::vector<scalar>& scatterPoints,
                        std::vector<scalar>& scatterValues)
    {
        scatterPoints_ = scatterPoints;
        scatterValues_ = scatterValues;
    }

    void setProfileData(std::initializer_list<scalar> scatterPoints,
                        std::initializer_list<scalar> scatterValues)
    {
        scatterPoints_.resize(scatterPoints.size());
        scatterValues_.resize(scatterValues.size());

        label i = 0;
        for (auto sp : scatterPoints)
        {
            scatterPoints_[i++] = sp;
        }

        i = 0;
        for (auto sv : scatterValues)
        {
            scatterValues_[i++] = sv;
        }
    }

    void setRawValue(label value)
    {
        rawLabelValue_ = value;
    }

    void setRawValue(scalar value)
    {
        rawScalarValue_ = value;
    }

    void setRawValue(std::string value)
    {
        rawStringValue_ = value;
    }

    void setRawValue(bool value)
    {
        rawBoolValue_ = value;
    }

private:
    inputDataType type_;

    // time interpolation
    timeInterpolationSchemeType interpType_;
    label interpOrder_;

    // mandatory and optional data members

    // constant value
    std::array<scalar, N> value_;

    // expressions
    std::array<std::string, N> expression_;

    // time table
    bool useInternalMethod_;
    std::vector<scalar> timeTable_;
    std::vector<std::unique_ptr<BSpline1D>> splines_;

    // profile data and idw parameters
    std::vector<scalar> scatterPoints_;
    std::vector<scalar> scatterValues_;
    scalar distancePower_ = 2.0;
    label donorCount_ = 4;

    // raw data containers
    label rawLabelValue_;
    scalar rawScalarValue_;
    std::string rawStringValue_;
    bool rawBoolValue_;

    void setupInterpolator_(const YAML::Node* node = nullptr)
    {
        if (node)
        {
            if ((*node)["interpolation_type"])
            {
                interpType_ = convertTimeInterpolationSchemeTypeFromString(
                    (*node)["interpolation_type"].template as<std::string>());
            }
            if ((*node)["interpolation_order"])
            {
                interpOrder_ =
                    (*node)["interpolation_order"].template as<label>();
            }
        }

        if (interpType_ == timeInterpolationSchemeType::BSpline)
        {
            assert(!timeTable_.empty());
            assert(timeTable_.size() % (N + 1) == 0);

            splines_.clear();
            const label n_rows = timeTable_.size() / (N + 1);
            for (label i = 0; i < N; i++)
            {
                auto indepSpan =
                    std::span<const scalar>(timeTable_).subspan(0, n_rows);
                auto depSpan = std::span<const scalar>(timeTable_)
                                   .subspan(n_rows * (i + 1), n_rows);
                splines_.push_back(std::make_unique<BSpline1D>(
                    interpOrder_,
                    std::vector<double>(indepSpan.begin(), indepSpan.end()),
                    std::vector<double>(depSpan.begin(), depSpan.end())));
            }
        }
        else if (interpType_ == timeInterpolationSchemeType::closest)
        {
            interpOrder_ = 0;
        }
        else if (interpType_ == timeInterpolationSchemeType::piecewiseLinear)
        {
            interpOrder_ = 1;
        }
    }
};

class dataHandler
{
private:
    std::unordered_map<std::string, baseInputData*> inputDataMap_;

    template <size_t N>
    void addInputData(std::string inputDataKey)
    {
        // Check if key already available
        if (inputDataMap_.find(inputDataKey) != inputDataMap_.end())
        {
            errorMsg("input key: " + inputDataKey + " already available");
        }

        // Instantiate and put in the map
        inputDataMap_[inputDataKey] = new inputData<N>();
    }

public:
    // Constructors
    dataHandler() = default;

    virtual ~dataHandler()
    {
        for (auto it : inputDataMap_)
        {
            if (it.second)
            {
                delete it.second;
            }
        }
    }

    // Access data base by key
    template <size_t N>
    const inputData<N>& data(std::string inputDataKey) const
    {
        // Check if key already available
        if (inputDataMap_.find(inputDataKey) == inputDataMap_.end())
        {
            errorMsg("input key: " + inputDataKey + " not available");
        }
        return dynamic_cast<const inputData<N>&>(
            *inputDataMap_.find(inputDataKey)->second);
    }

    // Access data base by key
    template <size_t N>
    inputData<N>& data(std::string inputDataKey)
    {
        // Check if key already available
        if (inputDataMap_.find(inputDataKey) == inputDataMap_.end())
        {
            errorMsg("input key: " + inputDataKey + " not available");
        }
        return dynamic_cast<inputData<N>&>(
            *inputDataMap_.find(inputDataKey)->second);
    }

    // Register data

    template <size_t N>
    void query(const YAML::Node node,
               std::string registeredKeyName = "value",
               std::string inputDataKey = "value",
               std::string inputTypeKey = "input_type")
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // query for values
        data<N>(registeredKeyName).query(node, inputDataKey, inputTypeKey);
    }

    template <size_t N>
    void queryMulti(const YAML::Node node,
                    std::string registeredKeyName,
                    std::vector<std::string> inputDataKeys,
                    std::string inputTypeKey = "input_type")
    {
        assert(inputDataKeys.size() == N);

        // add input data
        addInputData<N>(registeredKeyName);

        // query for values
        data<N>(registeredKeyName)
            .queryMulti(node, inputDataKeys, inputTypeKey);
    }

    template <size_t N>
    void queryConstantValue(const YAML::Node node,
                            std::string registeredKeyName,
                            std::string inputDataKey,
                            const std::string internalMethod = "")
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::constant);

        // query for values
        data<N>(registeredKeyName)
            .queryConstantValue(node, inputDataKey, internalMethod);
    }

    template <size_t N>
    void queryExpression(const YAML::Node node,
                         std::string registeredKeyName,
                         std::string inputDataKey)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::expression);

        // query for values
        data<N>(registeredKeyName).queryExpression(node, inputDataKey);
    }

    template <size_t N>
    void queryTimeTable(const YAML::Node node,
                        std::string registeredKeyName,
                        std::string filePathKey)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::timeTable);

        // query for values
        data<N>(registeredKeyName).queryTimeTable(node, filePathKey);
    }

    template <size_t N>
    void setNull(std::string registeredKeyName)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::null);
    }

    template <size_t N>
    void setConstantValue(std::string registeredKeyName,
                          std::vector<scalar> value)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::constant);

        // query for values
        data<N>(registeredKeyName).setValue(value);
    }

    template <size_t N>
    void setConstantValue(std::string registeredKeyName,
                          std::initializer_list<scalar> value)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::constant);

        // query for values
        data<N>(registeredKeyName).setValue(value);
    }

    template <size_t N>
    void addExpression(std::string registeredKeyName,
                       std::vector<std::string> expression)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::expression);

        // query for values
        data<N>(registeredKeyName).setExpression(expression);
    }

    template <size_t N>
    void addExpression(std::string registeredKeyName,
                       std::initializer_list<std::string> expression)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::expression);

        // query for values
        data<N>(registeredKeyName).setExpression(expression);
    }

    template <size_t N>
    void addTimeTable(std::string registeredKeyName,
                      std::vector<scalar>& timeTable)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::timeTable);

        // query for values
        data<N>(registeredKeyName).setTimeTable(timeTable);
    }

    template <size_t N>
    void addTimeTable(std::string registeredKeyName,
                      std::initializer_list<scalar> timeTable)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::timeTable);

        // query for values
        data<N>(registeredKeyName).setTimeTable(timeTable);
    }

    template <size_t N>
    void addProfileData(std::string registeredKeyName,
                        std::vector<scalar>& scatterPoints,
                        std::vector<scalar>& scatterValues)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::profileData);

        // query for values
        data<N>(registeredKeyName).setProfileData(scatterPoints, scatterValues);
    }

    template <size_t N>
    void addProfileData(std::string registeredKeyName,
                        std::initializer_list<scalar> scatterPoints,
                        std::initializer_list<scalar> scatterValues)
    {
        // add input data
        addInputData<N>(registeredKeyName);

        // set type
        data<N>(registeredKeyName).setType(inputDataType::profileData);

        // query for values
        data<N>(registeredKeyName).setProfileData(scatterPoints, scatterValues);
    }

    void addRawData(std::string registeredKeyName, label value)
    {
        // add input data
        addInputData<1>(registeredKeyName);

        // set type
        data<1>(registeredKeyName).setType(inputDataType::null);

        // store value
        data<1>(registeredKeyName).setRawValue(value);
    }

    void addRawData(std::string registeredKeyName, scalar value)
    {
        // add input data
        addInputData<1>(registeredKeyName);

        // set type
        data<1>(registeredKeyName).setType(inputDataType::null);

        // store value
        data<1>(registeredKeyName).setRawValue(value);
    }

    void addRawData(std::string registeredKeyName, std::string value)
    {
        // add input data
        addInputData<1>(registeredKeyName);

        // set type
        data<1>(registeredKeyName).setType(inputDataType::null);

        // store value
        data<1>(registeredKeyName).setRawValue(value);
    }

    void addRawData(std::string registeredKeyName, const char* value)
    {
        // add input data
        addInputData<1>(registeredKeyName);

        // set type
        data<1>(registeredKeyName).setType(inputDataType::null);

        // store value
        data<1>(registeredKeyName).setRawValue(std::string(value));
    }

    void addRawData(std::string registeredKeyName, bool value)
    {
        // add input data
        addInputData<1>(registeredKeyName);

        // set type
        data<1>(registeredKeyName).setType(inputDataType::null);

        // store value
        data<1>(registeredKeyName).setRawValue(value);
    }

    label rawLabelValue(std::string registeredKeyName) const
    {
        return data<1>(registeredKeyName).rawLabelValue();
    }

    scalar rawScalarValue(std::string registeredKeyName) const
    {
        return data<1>(registeredKeyName).rawScalarValue();
    }

    std::string rawStringValue(std::string registeredKeyName) const
    {
        return data<1>(registeredKeyName).rawStringValue();
    }

    bool rawBoolValue(std::string registeredKeyName) const
    {
        return data<1>(registeredKeyName).rawBoolValue();
    }

    bool empty() const
    {
        return inputDataMap_.empty();
    }

    bool isInputDataAdded(std::string inputDataKey) const
    {
        return (inputDataMap_.find(inputDataKey) != inputDataMap_.end());
    }
};

} // namespace accel

#undef FIDX

#include "dataHandler.hpp"

#endif // DATAHANDLER_H
