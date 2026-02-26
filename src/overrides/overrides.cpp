// File       : overrides.cpp
// Created    : Mon Feb 09 2026 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "overrides.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

overrides::overrides(realm* realm) : realmPtr_(realm)
{
}

void overrides::read(const YAML::Node& inputNode)
{
}

} // namespace accel
