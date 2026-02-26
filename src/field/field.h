// File       : field.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base templated field class wrapping STK mesh field data
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIELD_H
#define FIELD_H

// code
#include "mesh.h"
#include "messager.h"
#include "zone.h"

namespace accel
{

// non-templated base field class for the purpose
// of puting fields in a single map
class baseField
{
protected:
    // assembly/solution-related values: only useful for physical scalar-based
    // fields

    scalar urf_;

    scalar clipFactor_;

public:
    baseField() : urf_(1.0), clipFactor_(0.75)
    {
    }

    virtual ~baseField() = default;

    // Getters and setters

    scalar urf() const
    {
        return urf_;
    };

    void setURF(scalar urf)
    {
        urf_ = urf;
    };

    scalar clipFactor() const
    {
        return clipFactor_;
    };

    void setClipFactor(scalar fac)
    {
        clipFactor_ = fac;
    };
};

class zone;

template <class T, size_t N>
class field : public baseField
{
protected:
    mesh* meshPtr_;

    std::vector<label> activeZones_;

    std::string name_;

    // The original stk field, stored in stk data structures
    stk::mesh::Field<T>* stkFieldPtr_ = nullptr;

    // stats-related values

    T scale_;

    T maxScale_;

    T minScale_;

    T minAcceptedValue_;

    T maxAcceptedValue_;

public:
    static constexpr size_t NComponents = N;
    using DataType = T;

    // Constructors

    // main constructors
    field(mesh* meshPtr,
          stk::topology::rank_t rank,
          std::string name,
          unsigned numberOfStates = 1);

    // special constructor to create a field from an already existing stk field:
    // this is useful if the functionalities of the class field are to be
    // employed to an stk field, i.e. overloaded operators
    field(mesh* meshPtr, stk::mesh::Field<T>* stkField_ptr);

    // Destructors
    ~field();

    // Operators

    inline T& operator[](const stk::mesh::Entity& entity)
    {
        return *stk::mesh::field_data(*stkFieldPtr_, entity);
    };

    inline const T& operator[](const stk::mesh::Entity& entity) const
    {
        return *stk::mesh::field_data(*stkFieldPtr_, entity);
    };

    // Methods

    void synchronizeGhostedEntities(label iZone);

    void synchronizeGhostedEntities(label iZone, label iBoundary);

    // Getters and setters
    virtual void updateScale(void)
    {
        errorMsg("Should not reach here");
    };

    // Access

    std::string name() const
    {
        return name_;
    };

    inline stk::mesh::Field<T>* stkFieldPtr()
    {
        return stkFieldPtr_;
    }

    inline const stk::mesh::Field<T>* stkFieldPtr() const
    {
        return stkFieldPtr_;
    }

    inline stk::mesh::Field<T>& stkFieldRef()
    {
        return *stkFieldPtr_;
    }

    inline const stk::mesh::Field<T>& stkFieldRef() const
    {
        return *stkFieldPtr_;
    }

    stk::mesh::BulkData& bulkDataRef();

    const stk::mesh::BulkData& bulkDataRef() const;

    stk::mesh::MetaData* metaDataPtr();

    const stk::mesh::MetaData* metaDataPtr() const;

    stk::mesh::MetaData& metaDataRef();

    const stk::mesh::MetaData& metaDataRef() const;

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;

    // Getters and setters

    T scale() const
    {
        return scale_;
    };

    void setScale(T scale)
    {
        scale_ = scale;
    };

    T maxScale() const
    {
        return maxScale_;
    };

    void setMaxScale(T max)
    {
        maxScale_ = max;
    };

    T minScale() const
    {
        return minScale_;
    };

    void setMinScale(T min)
    {
        minScale_ = min;
    };

    T minAcceptedValue() const
    {
        return minAcceptedValue_;
    };

    void setMinAcceptedValue(T min)
    {
        minAcceptedValue_ = min;
    };

    T maxAcceptedValue() const
    {
        return maxAcceptedValue_;
    };

    void setMaxAcceptedValue(T max)
    {
        maxAcceptedValue_ = max;
    };

    virtual void setZone(label iZone)
    {
        assert(activeZones_.size() > iZone);
        activeZones_[iZone] = 1;
    }

    bool isZoneSet(label iZone) const
    {
        assert(activeZones_.size() > iZone);
        return activeZones_[iZone] == 1;
    }

    bool isZoneUnset(label iZone) const
    {
        assert(activeZones_.size() > iZone);
        return activeZones_[iZone] == 0;
    }
};

using simpleScalarField = field<scalar, 1>;
using simpleVectorField = field<scalar, SPATIAL_DIM>;
using simpleTensorField = field<scalar, SPATIAL_DIM * SPATIAL_DIM>;

} // namespace accel

#include "field.hpp"

#endif // FIELD_H
