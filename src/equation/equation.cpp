// File       : equation.cpp
// Created    : Fri Jan 26 2024 09:07:47 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "equation.h"
#include "boundary.h"
#include "convergenceAcceleration.h"
#include "mesh.h"
#include "simulation.h"
#include "zone.h"

#include <cctype>

namespace accel
{

namespace
{
std::string equationKeyFromName_(const std::string& name)
{
    std::string key;
    key.reserve(name.size());
    for (unsigned char c : name)
    {
        if (std::isalnum(c))
        {
            key.push_back(static_cast<char>(std::tolower(c)));
        }
        else if (c == ' ' || c == '-' || c == '_')
        {
            if (key.empty() || key.back() != '_')
            {
                key.push_back('_');
            }
        }
    }
    if (!key.empty() && key.back() == '_')
    {
        key.pop_back();
    }
    return key;
}
} // namespace

stk::mesh::PartVector equation::collectInactiveInteriorParts()
{
    // get mesh
    auto& mesh = domainVector_[0].get()->meshRef();

    // collect first all parts associated to active domains
    stk::mesh::PartVector incPartVec, excPartVec;
    for (const auto& domain : domainVector_)
    {
        for (auto part : domain->zonePtr()->interiorParts())
        {
            incPartVec.push_back(part);
        }
    }

    // collect all other remaining parts
    for (auto part : mesh.interiorActiveParts())
    {
        auto it = std::find(incPartVec.begin(), incPartVec.end(), part);
        if (it == incPartVec.end())
        {
            excPartVec.push_back(const_cast<stk::mesh::Part*>(part));
        }
    }

    return excPartVec;
}

stk::mesh::PartVector equation::collectInteriorParts()
{
    // get mesh
    auto& mesh = domainVector_[0].get()->meshRef();

    // collect first all parts associated to active domains
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
        for (auto part : domain->zonePtr()->interiorParts())
        {
            incPartVec.push_back(part);
        }
    }

    return incPartVec;
}

stk::mesh::PartVector equation::collectBoundaryParts()
{
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            for (auto part : boundaryRef.parts())
            {
                incPartVec.push_back(part);
            }
        }
    }

    return incPartVec;
}

stk::mesh::PartVector equation::collectStationaryParts()
{
    // get mesh
    auto& mesh = domainVector_[0].get()->meshRef();

    // collect first all parts associated to active domains
    stk::mesh::PartVector statPartVec;
    for (const auto& domain : domainVector_)
    {
        for (auto part : domain->zonePtr()->stationaryParts())
        {
            statPartVec.push_back(part);
        }
    }

    return statPartVec;
}

void equation::initializeAcceleration_()
{
    if (accelerationInitialized_)
    {
        return;
    }
    accelerationInitialized_ = true;

    if (domainVector_.empty())
    {
        return;
    }

    const auto& sim = domainVector_[0]->simulationRef();
    const auto solverNode = sim.getYAMLSolverControlNode();
    if (!solverNode || !solverNode["advanced_options"] ||
        !solverNode["advanced_options"]["equation_controls"] ||
        !solverNode["advanced_options"]["equation_controls"]["acceleration"])
    {
        return;
    }

    const auto accelNode =
        solverNode["advanced_options"]["equation_controls"]["acceleration"];
    const std::string eqKey = equationKeyFromName_(name_);
    if (!accelNode[eqKey])
    {
        return;
    }

    const auto eqAccel = accelNode[eqKey];
    if (!eqAccel["option"])
    {
        return;
    }

    convergenceAcceleration::Config cfg;
    cfg.type = convertAccelerationTypeFromString(
        eqAccel["option"].template as<std::string>());
    if (cfg.type == accelerationType::none)
    {
        return;
    }

    if (eqAccel["initial_omega"])
    {
        cfg.aitkenInitialOmega = eqAccel["initial_omega"].template as<scalar>();
    }
    if (eqAccel["omega_min"])
    {
        cfg.aitkenOmegaMin = eqAccel["omega_min"].template as<scalar>();
    }
    if (eqAccel["omega_max"])
    {
        cfg.aitkenOmegaMax = eqAccel["omega_max"].template as<scalar>();
    }
    if (eqAccel["iqn_ils_window"])
    {
        cfg.iqnIlsWindow = eqAccel["iqn_ils_window"].template as<label>();
    }
    if (eqAccel["iqn_ils_regularization"])
    {
        cfg.iqnIlsRegularization =
            eqAccel["iqn_ils_regularization"].template as<scalar>();
    }

    accelerationPtr_ = std::make_unique<convergenceAcceleration>(cfg);
}

const Vector& equation::applyAcceleration_(const Vector& correction,
                                           const scalar relaxValue,
                                           scalar& outRelaxValue)
{
    outRelaxValue = relaxValue;

    initializeAcceleration_();
    if (!accelerationPtr_ || !accelerationPtr_->enabled())
    {
        return correction;
    }

    if (domainVector_.empty())
    {
        return correction;
    }

    auto& sim = domainVector_[0]->simulationRef();
    const label timeStep = sim.controlsRef().getTimeStepCount();
    if (timeStep != lastAccelTimeStep_)
    {
        accelerationPtr_->resetForTimeStep();
        lastAccelTimeStep_ = timeStep;
        lastAccelIter_ = -1;
        lastAccelCorrectionPtr_ = nullptr;
        lastAccelUsesScratch_ = false;
    }

    const label iter = sim.controlsRef().iter;
    if (lastAccelCorrectionPtr_ == &correction && lastAccelIter_ == iter)
    {
        outRelaxValue = lastAccelRelaxValue_;
        return lastAccelUsesScratch_ ? accelerationScratch_[0] : correction;
    }

    const Vector& result = accelerationPtr_->apply(
        correction, relaxValue, accelerationScratch_, outRelaxValue);

    lastAccelCorrectionPtr_ = &correction;
    lastAccelIter_ = iter;
    lastAccelRelaxValue_ = outRelaxValue;
    lastAccelUsesScratch_ = (&result != &correction);

    return result;
}

} /* namespace accel */
