// File : displacement.h
// Created : Sat Dec 06 2025 10:22:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Mesh displacement vector field with deforming mesh coordinate
// handling
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include "nodeVectorField.h"

namespace accel
{

class displacement : public nodeVectorField
{
public:
    // Constructors

    displacement(realm* realmPtr,
                 const std::string name,
                 unsigned numberOfStates);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;

protected:
    std::string getCoordinatesID_(label iZone) const override
    {
        auto* zonePtr = this->meshPtr()->zonePtr(iZone);

        bool relDisp = true;
        if (zonePtr->meshDeforming() &&
            !zonePtr->deformationRef().displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::coordinates_ID : mesh::original_coordinates_ID;
    }

    std::string getDualNodalVolumeID_(label iZone) const override
    {
        auto* zonePtr = this->meshPtr()->zonePtr(iZone);

        bool relDisp = true;
        if (zonePtr->meshDeforming() &&
            !zonePtr->deformationRef().displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::dual_nodal_volume_ID
                       : mesh::original_dual_nodal_volume_ID;
    }

    std::string getExposedAreaVectorID_(label iZone) const override
    {
        auto* zonePtr = this->meshPtr()->zonePtr(iZone);

        bool relDisp = true;
        if (zonePtr->meshDeforming() &&
            !zonePtr->deformationRef().displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::exposed_area_vector_ID
                       : mesh::original_exposed_area_vector_ID;
    }
};

} // namespace accel

#endif // DISPLACEMENT_H
