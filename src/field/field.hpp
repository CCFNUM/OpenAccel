// File : field.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

// Constructors

template <class T, size_t N>
field<T, N>::field(mesh* meshPtr,
                   stk::topology::rank_t rank,
                   std::string name,
                   unsigned numberOfStates)
    : meshPtr_(meshPtr), name_(name), scale_(zero<T>()), maxScale_(zero<T>()),
      minScale_(zero<T>()), minAcceptedValue_(std::numeric_limits<T>::lowest()),
      maxAcceptedValue_(std::numeric_limits<T>::max()),
      activeZones_(std::vector<label>(meshPtr->nZones(), 0))
{
    static_assert(std::is_arithmetic<T>::value, "T must be a numeric type");

    // Create the stk field
    stkFieldPtr_ = &(meshPtr_->metaDataPtr()->template declare_field<T>(
        rank, name, numberOfStates));

    // Set field output type
    stk::io::set_field_output_type(this->stkFieldRef(), fieldType[N]);
}

template <class T, size_t N>
field<T, N>::field(mesh* meshPtr, stk::mesh::Field<T>* stkField_ptr)
    : meshPtr_(meshPtr), stkFieldPtr_(stkField_ptr),
      name_(stkField_ptr->name()), scale_(zero<T>()), maxScale_(zero<T>()),
      minScale_(zero<T>()), minAcceptedValue_(std::numeric_limits<T>::lowest()),
      maxAcceptedValue_(std::numeric_limits<T>::max()),
      activeZones_(std::vector<label>(meshPtr->nZones(), 0))
{
    // Set field output type
    stk::io::set_field_output_type(this->stkFieldRef(), fieldType[N]);
}

// Destructors

template <class T, size_t N>
field<T, N>::~field()
{
}

// Operations

template <class T, size_t N>
void field<T, N>::synchronizeGhostedEntities(label iZone)
{
    // ensure zone is active
    assert(activeZones_[iZone] == 1);

    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(this->bulkDataRef(),
                                          {this->stkFieldPtr_});
    }
}

template <class T, size_t N>
void field<T, N>::synchronizeGhostedEntities(label iZone, label iBoundary)
{
    // ensure zone is active
    assert(activeZones_[iZone] == 1);

    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(this->bulkDataRef(),
                                          {this->stkFieldPtr_});
    }
}

// Access

template <class T, size_t N>
stk::mesh::BulkData& field<T, N>::bulkDataRef()
{
    return meshPtr_->metaDataPtr()->mesh_bulk_data();
}

template <class T, size_t N>
const stk::mesh::BulkData& field<T, N>::bulkDataRef() const
{
    return meshPtr_->metaDataPtr()->mesh_bulk_data();
}

template <class T, size_t N>
stk::mesh::MetaData* field<T, N>::metaDataPtr()
{
    return meshPtr_->metaDataPtr();
}

template <class T, size_t N>
const stk::mesh::MetaData* field<T, N>::metaDataPtr() const
{
    return meshPtr_->metaDataPtr();
}

template <class T, size_t N>
stk::mesh::MetaData& field<T, N>::metaDataRef()
{
    return *meshPtr_->metaDataPtr();
}

template <class T, size_t N>
const stk::mesh::MetaData& field<T, N>::metaDataRef() const
{
    return *meshPtr_->metaDataPtr();
}

template <class T, size_t N>
mesh* field<T, N>::meshPtr()
{
    return meshPtr_;
}

template <class T, size_t N>
const mesh* field<T, N>::meshPtr() const
{
    return meshPtr_;
}

template <class T, size_t N>
mesh& field<T, N>::meshRef()
{
    return *meshPtr_;
}

template <class T, size_t N>
const mesh& field<T, N>::meshRef() const
{
    return *meshPtr_;
}

} // namespace accel
