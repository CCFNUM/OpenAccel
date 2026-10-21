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

#ifdef HAS_UTOPIA
    bool fsi_ = false;
#endif /* HAS_UTOPIA */

public:
    overrides(realm* realm);

    // operations

    void read(const YAML::Node& inputNode);

#ifdef HAS_UTOPIA
    // access

    bool fsi() const
    {
        return fsi_;
    }

    YAML::Node getFSINode() const;
#endif /* HAS_UTOPIA */
};

} // namespace accel

#endif // OVERRIDES_H
