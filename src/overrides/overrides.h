// File       : overrides.h
// Created    : Mon Feb 09 2026 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Runtime configuration overrides including fluid-structure
//              interaction
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef OVERRIDES_H
#define OVERRIDES_H

// code
#include "masking.h"
#include "types.h"

#include <memory>

namespace accel
{

class realm;

class overrides
{
private:
    realm* realmPtr_;

    bool fsi_ = false;

    // masking a region of the mesh is an override of what the solver would
    // otherwise do there, so the bodies live with the other overrides
    std::unique_ptr<masking> maskingPtr_ = nullptr;

public:
    overrides(realm* realm);

    // operations

    void read(const YAML::Node& inputNode);

    void initialize();

    // access

    bool fsi() const
    {
        return fsi_;
    }

    YAML::Node getFSINode() const;

    masking* maskingPtr()
    {
        return maskingPtr_.get();
    }

    const masking* maskingPtr() const
    {
        return maskingPtr_.get();
    }
};

} // namespace accel

#endif // OVERRIDES_H
