// File : sideField.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Field defined at integration points on boundary sides
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIDEFIELD_H
#define SIDEFIELD_H

// code
#include "field.h"
#include "master_element/MasterElement.h"

namespace accel
{

template <class T, size_t N>
class nodeSideField;

template <class T, size_t N>
class sideField : public field<T, N>
{
protected:
    std::unique_ptr<sideField<T, N>> prevTimeFieldPtr_ = nullptr;

    // set the values given an array (T* arr): internal use only

    void setToValue_(const T* val);

    void setToValue_(const T* val, const stk::mesh::Part& part);

    void setToValue_(const T* val, const stk::mesh::PartVector parts);

    virtual std::string getCoordinatesID_(label iZone) const
    {
        return mesh::coordinates_ID;
    }

    virtual std::string getDualNodalVolumeID_(label iZone) const
    {
        return mesh::dual_nodal_volume_ID;
    }

    virtual std::string getExposedAreaVectorID_(label iZone) const
    {
        return mesh::exposed_area_vector_ID;
    }

public:
    // Constructors

    sideField(mesh* meshPtr, std::string name, unsigned numberOfStates = 1);

    sideField(mesh* meshPtr, stk::mesh::Field<T>* stkField_ptr);

    // Operations

    void putFieldOnPart(const stk::mesh::Part& part);

    bool definedOn(const stk::mesh::Part& part) const;

    bool definedOn(const stk::mesh::PartVector parts) const;

    bool definedOnBoundary(label iZone, label iBoundary) const;

    // set the values given various inputs

    void setToValue(std::initializer_list<T> val);

    void setToValue(std::initializer_list<T> val, const stk::mesh::Part& part);

    void setToValue(std::initializer_list<T> val,
                    const stk::mesh::PartVector parts);

    void setToValue(std::vector<T> val);

    void setToValue(std::vector<T> val, const stk::mesh::Part& part);

    void setToValue(std::vector<T> val, const stk::mesh::PartVector parts);

    void setToValue(std::array<T, N> val);

    void setToValue(std::array<T, N> val, const stk::mesh::Part& part);

    void setToValue(std::array<T, N> val, const stk::mesh::PartVector parts);

    // interpolate a node-side field to the current side field
    void interpolate(const nodeSideField<T, N>& nsf,
                     label iZone,
                     label iBoundary,
                     bool shifted);

    void interpolate(const nodeSideField<T, N>& nsf,
                     label iInterface,
                     bool master,
                     bool shifted);

    // Access

    sideField<T, N>& prevTimeRef();

    const sideField<T, N>& prevTimeRef() const;
};

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const sideField<T, N>& field);

// Specializations

template <>
void sideField<scalar, 1>::interpolate(const nodeSideField<scalar, 1>& nsf,
                                       label iZone,
                                       label iBoundary,
                                       bool shifted);

template <>
void sideField<scalar, SPATIAL_DIM>::interpolate(
    const nodeSideField<scalar, SPATIAL_DIM>& nsf,
    label iZone,
    label iBoundary,
    bool shifted);

} // namespace accel

#include "sideField.hpp"

#endif // SIDEFIELD_H
