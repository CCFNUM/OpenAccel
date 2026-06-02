// File       : compressibility.h
// Created    : Fri Apr 11 2025 18:10:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for fluid compressibility
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef COMPRESSIBILITY_H
#define COMPRESSIBILITY_H

#include "property.h"

namespace accel
{

class compressibility : public property
{
public:
    // Constructors

    compressibility(realm* realmPtr,
                    const std::string name,
                    unsigned numberOfStates,
                    bool highResolution,
                    bool computeGradient);
};

} // namespace accel

#endif // COMPRESSIBILITY_H
