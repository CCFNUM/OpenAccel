// File : elementField.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Templated base class for element fields with side fields and
// time history
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied
// Sciences and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef ELEMENTFIELD_H
#define ELEMENTFIELD_H

// code
#include "boundary.h"
#include "field.h"
#include "master_element/MasterElement.h"
#include "mesh.h"
#include "messager.h"
#include "sideField.h"
#include "zone.h"

namespace accel
{

class mesh;

template <class T, size_t N>
class elementField : public field<T, N>
{
protected:
    std::unique_ptr<elementField<T, N>> prevTimeFieldPtr_ = nullptr;

    std::unique_ptr<sideField<T, N>> sideFieldPtr_ = nullptr;

    // flag to mark initialzed field at a zone
    std::vector<bool> isInitialized_;

    virtual void putFieldOnRegisteredParts_();

    // set the values given an array (T* arr): internal use only

    void setToValue_(const T* val);

    void setToValue_(const T* val, const stk::mesh::Part& part);

    void setToValue_(const T* val, const stk::mesh::PartVector parts);

public:
    // Constructors

    elementField(mesh* meshPtr, std::string name, unsigned numberOfStates = 1);

    elementField(mesh* meshPtr, stk::mesh::Field<T>* stkField_ptr);

    // Operations

    void putFieldOnPart(const stk::mesh::Part& part);

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

    virtual void registerSideField(label iZone, label iBoundary);

    elementField& operator=(const elementField& fld);

    // Initialize

    virtual void initialize(label iZone, bool force);

    virtual void initializeField(label iZone);

    virtual void initializeSideField(label iZone);

    virtual void initializeBoundarySideField(label iZone, label iBoundary);

    // Update

    virtual void update(label iZone);

    virtual void updateField(label iZone);

    virtual void updateSideFields(label iZone);

    virtual void updateBoundarySideField(label iZone, label iBoundary);

    // Getters and setters

    bool isInitialized(label iZone) const
    {
        assert(isInitialized_.size() > iZone);
        assert(this->isZoneSet(iZone));
        return isInitialized_[iZone];
    }

    void setIsInitialized(label iZone)
    {
        assert(isInitialized_.size() > iZone);
        assert(this->isZoneSet(iZone));
        isInitialized_[iZone] = true;
    }

    // Access

    elementField<T, N>& prevTimeRef();

    const elementField<T, N>& prevTimeRef() const;

    sideField<T, N>& sideFieldRef();

    const sideField<T, N>& sideFieldRef() const;

    sideField<T, N>* sideFieldPtr();

    const sideField<T, N>* sideFieldPtr() const;
};

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const elementField<T, N>& field);

} // namespace accel

#include "elementField.hpp"

#endif // ELEMENTFIELD_H
