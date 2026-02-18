// File : equation.cpp
// Created : Fri Jan 26 2024 09:07:47 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "equation.h"
#include "boundary.h"
#include "mesh.h"
#include "zone.h"

namespace accel
{

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

} /* namespace accel */
