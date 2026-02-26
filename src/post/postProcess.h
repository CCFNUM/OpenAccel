// File       : postProcess.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Post-processing manager for reductions, forces, and field probes
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "fieldBroker.h"

namespace accel
{

class postProcessObject;

class postProcess : public fieldBroker
{
private:
    fs::path directory_;

    std::vector<std::unique_ptr<postProcessObject>> objectPtrArray_;

    label iter_ = 0;

    scalar time_ = 0.0;

public:
    // Constructors
    postProcess(realm* realm, fs::path directory);

    // IO

    void read(const YAML::Node& inputNode);

    void update();

    // Access

    fs::path directory() const
    {
        return directory_;
    }

    label iter() const
    {
        return iter_;
    }

    scalar time() const
    {
        return time_;
    }

    // Methods

    std::string instance();

    std::string instanceHeader();

    // Enable public access for required fields
    using fieldBroker::pRef;
    using fieldBroker::wallShearStressRef;
};

class postProcessObject
{
public:
    postProcessObject(postProcess* postProcessManagerPtr,
                      std::string name,
                      postProcessType type,
                      std::vector<std::string> location,
                      label frequency,
                      bool writeToFile);

    virtual void update() = 0;

protected:
    friend class postProcess;

    postProcess* postProcessPtr_;

    std::string name_;

    postProcessType type_;

    std::vector<std::string> location_;

    label frequency_;

    bool writeToFile_;
};

class reductionObject : public postProcessObject
{
private:
    reductionType reductionType_;

    std::string field_;

public:
    reductionObject(postProcess* postProcessManagerPtr,
                    std::string name,
                    postProcessType type,
                    std::vector<std::string> location,
                    label frequency,
                    bool writeToFile,
                    reductionType rType,
                    std::string field);

    void update() override;
};

class forceObject : public postProcessObject
{
private:
    bool calculateMoment_ = false;

    std::array<scalar, SPATIAL_DIM> momentCenter_ = {0.0};

    // total print vs. viscous and pressure force components
    bool totalPrint_ = true;

public:
    forceObject(postProcess* postProcessManagerPtr,
                std::string name,
                postProcessType type,
                std::vector<std::string> location,
                label frequency,
                bool writeToFile,
                bool calculateMoment,
                std::array<scalar, SPATIAL_DIM> momentCenter,
                bool totalPrint);

    void update() override;
};

class probeObject : public postProcessObject
{
private:
    std::array<scalar, SPATIAL_DIM> probeLocation_;

    std::string field_;

    stk::mesh::EntityId encapsulatingElementIdent_ = stk::mesh::InvalidEntityId;

public:
    probeObject(postProcess* postProcessManagerPtr,
                std::string name,
                postProcessType type,
                std::vector<std::string> location,
                label frequency,
                bool writeToFile,
                std::array<scalar, SPATIAL_DIM> probeLocation,
                std::string field);

    void update() override;
};

} // namespace accel

#endif // POSTPROCESS_H
