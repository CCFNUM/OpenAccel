// File       : nodeSideField.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Field defined at nodes on boundary sides with interpolation
// support
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NODESIDEFIELD_H
#define NODESIDEFIELD_H

// code
#include "field.h"
#include "master_element/MasterElement.h"

namespace accel
{

template <class T, size_t N>
class sideField;

template <class T, size_t N>
class nodeSideField : public field<T, N>
{
protected:
    std::unique_ptr<nodeSideField<T, N>> prevTimeFieldPtr_ = nullptr;

    std::unique_ptr<simpleScalarField> localAssembledAreaFieldPtr_ = nullptr;

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

    nodeSideField(mesh* meshPtr, std::string name, unsigned numberOfStates = 1);

    nodeSideField(mesh* meshPtr, stk::mesh::Field<T>* stkField_ptr);

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

    void interpolate(const sideField<T, N>& sf, label iZone, label iBoundary);

#ifdef HAS_INTERFACE
    void interpolate(const sideField<T, N>& sf, label iInterface, bool master);
#endif /* HAS_INTERFACE */

    // Access

    nodeSideField<T, N>& prevTimeRef();

    const nodeSideField<T, N>& prevTimeRef() const;
};

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const nodeSideField<T, N>& field);

// Specializations

template <>
void nodeSideField<scalar, 1>::interpolate(const sideField<scalar, 1>& sf,
                                           label iZone,
                                           label iBoundary);

#ifdef HAS_INTERFACE
template <>
void nodeSideField<scalar, 1>::interpolate(const sideField<scalar, 1>& sf,
                                           label iInterface,
                                           bool master);
#endif /* HAS_INTERFACE */

template <>
void nodeSideField<scalar, SPATIAL_DIM>::interpolate(
    const sideField<scalar, SPATIAL_DIM>& sf,
    label iZone,
    label iBoundary);

#ifdef HAS_INTERFACE
template <>
void nodeSideField<scalar, SPATIAL_DIM>::interpolate(
    const sideField<scalar, SPATIAL_DIM>& sf,
    label iInterface,
    bool master);
#endif /* HAS_INTERFACE */

} // namespace accel

#include "nodeSideField.hpp"

#endif // nodalSideField_h
