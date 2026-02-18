// File : multiphaseModel.h
// Created : Sun Jan 26 2025 22:02:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base multiphase model managing phase definitions and volume
// fractions
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MULTIPHASEMODEL_H
#define MULTIPHASEMODEL_H

// code
#include "flowModel.h"

namespace accel
{

struct phase
{
    phase() = default;

    phase(label index, std::string name) : index_(index), name_(name)
    {
    }

    phase(label index, std::string name, bool primaryPhase)
        : index_(index), name_(name), primaryPhase_(primaryPhase)
    {
    }

    // global index in global set of materials in the simulation
    label index_ = -1;

    std::string name_;

    bool primaryPhase_ = true;
};

class multiphaseModel : public flowModel
{
protected:
    std::vector<phase> phases_;

public:
    // Constructors

    multiphaseModel(realm* realm);

    // Enable public use

    using fieldBroker::alphaRef;

    // Methods

    label phaseIndex(label iPhase) const
    {
        assert(iPhase < nPhases());
        return phases_[iPhase].index_;
    }

    label nPhases() const
    {
        return phases_.size();
    }

    phase& phaseRef(label iPhase)
    {
        return phases_[iPhase];
    }

    const phase& phaseRef(label iPhase) const
    {
        return phases_[iPhase];
    }
};

} /* namespace accel */

#endif // MULTIPHASEMODEL_H
