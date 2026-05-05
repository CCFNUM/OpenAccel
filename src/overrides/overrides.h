// File       : overrides.h
// Created    : Mon Feb 09 2026 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Runtime configuration overrides including fluid-structure
//              interaction
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef OVERRIDES_H
#define OVERRIDES_H

// code
#include "types.h"

namespace accel
{

class realm;

class overrides
{
private:
    realm* realmPtr_;

public:
    overrides(realm* realm);

    // operations

    void read(const YAML::Node& inputNode);
};

} // namespace accel

#endif // OVERRIDES_H
