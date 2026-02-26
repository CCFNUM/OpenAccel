// File       : model.h
// Created    : Wed Jan 07 2026
// Author     : Mhamad Mahdi Alloush
// Description: Base model class providing mesh coordinate and geometry access
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MODEL_H
#define MODEL_H

// code
#include "fieldBroker.h"

namespace accel
{

class model : public fieldBroker
{
public:
    model(realm* realm);

protected:
    virtual std::string
    getCoordinatesID_(const std::shared_ptr<domain> /*domain*/) const
    {
        return mesh::coordinates_ID;
    }

    virtual std::string
    getDualNodalVolumeID_(const std::shared_ptr<domain> /*domain*/) const
    {
        return mesh::dual_nodal_volume_ID;
    }

    virtual std::string
    getExposedAreaVectorID_(const std::shared_ptr<domain> /*domain*/) const
    {
        return mesh::exposed_area_vector_ID;
    }
};

} /* namespace accel */

#endif // MODEL_H
