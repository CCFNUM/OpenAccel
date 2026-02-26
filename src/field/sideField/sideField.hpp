// File       : sideField.hpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Template and inline implementations for side field data access
// and operations.
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <class T, size_t N>
sideField<T, N>::sideField(mesh* meshPtr,
                           std::string name,
                           unsigned numberOfStates)
    : field<T, N>(meshPtr,
                  meshPtr->metaDataRef().side_rank(),
                  name,
                  numberOfStates)
{
    // Populate previous time step fields if required: special constructors are
    // employed here
    if (numberOfStates > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<sideField<T, N>>(
            meshPtr, &this->stkFieldRef().field_of_state(stk::mesh::StateN));
    }
}

template <class T, size_t N>
sideField<T, N>::sideField(mesh* meshPtr, stk::mesh::Field<T>* stkFieldPtr)
    : field<T, N>(meshPtr, stkFieldPtr)
{
    if (stkFieldPtr->number_of_states() - int(stkFieldPtr->state()) > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<sideField<T, N>>(
            meshPtr,
            &this->stkFieldRef().field_of_state(
                (enum stk::mesh::FieldState)((int)stkFieldPtr->state() + 1)));
    }
}

template <class T, size_t N>
void sideField<T, N>::putFieldOnPart(const stk::mesh::Part& part)
{
    // A part could have different topologies, by getting the subparts
    // we make sure we take into account all topologies.
    for (const stk::mesh::Part* subPart : part.subsets())
    {
        // Determine number of integration points in the face
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(subPart->topology());
        const label numScsBip = meFC->numIntPoints_;

        // initial value
        T zeroVal[N * numScsBip];
        for (label i = 0; i < N * numScsBip; i++)
        {
            zeroVal[i] = static_cast<T>(0);
        }

        // Put the field on mesh
        stk::mesh::put_field_on_mesh(
            this->stkFieldRef(), *subPart, N * numScsBip, &zeroVal[0]);
    }
}

template <class T, size_t N>
bool sideField<T, N>::definedOn(const stk::mesh::Part& part) const
{
    for (const stk::mesh::Part* subPart : part.subsets())
    {
        if (this->stkFieldRef().defined_on(*subPart))
        {
            return true;
        }
    }
    return false;
}

template <class T, size_t N>
bool sideField<T, N>::definedOn(const stk::mesh::PartVector parts) const
{
    for (const stk::mesh::Part* part : parts)
    {
        if (!this->definedOn(*part))
        {
            return false;
        }
    }
    return true;
}

template <class T, size_t N>
bool sideField<T, N>::definedOnBoundary(label iZone, label iBoundary) const
{
    if (!this->definedOn(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts()))
    {
        return false;
    }
    return true;
}

template <class T, size_t N>
void sideField<T, N>::setToValue_(const T* val)
{
    // select all sides relevant to the field
    stk::mesh::Selector selAllSides =
        this->metaDataRef().universal_part() &
        stk::mesh::selectField(this->stkFieldRef());
    stk::mesh::BucketVector const& sideBuckets =
        this->bulkDataRef().get_buckets(this->metaDataRef().side_rank(),
                                        selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;
        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        // extract master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());

        // extract master element specifics
        const label numScsBip = meFC->numIntPoints_;

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            T* value =
                stk::mesh::field_data(this->stkFieldRef(), sideBucket, iSide);
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                for (label i = 0; i < N; i++)
                {
                    value[ip * N + i] = val[i];
                }
            }
        }
    }
}

template <class T, size_t N>
void sideField<T, N>::setToValue_(const T* val, const stk::mesh::Part& part)
{
    // select all sides
    stk::mesh::Selector selAllSides =
        this->metaDataRef().universal_part() & part;
    stk::mesh::BucketVector const& sideBuckets =
        this->bulkDataRef().get_buckets(this->metaDataRef().side_rank(),
                                        selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;
        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        // extract master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());

        // extract master element specifics
        const label numScsBip = meFC->numIntPoints_;

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            T* value =
                stk::mesh::field_data(this->stkFieldRef(), sideBucket, iSide);
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                for (label i = 0; i < N; i++)
                {
                    value[ip * N + i] = val[i];
                }
            }
        }
    }
}

template <class T, size_t N>
void sideField<T, N>::setToValue_(const T* val,
                                  const stk::mesh::PartVector parts)
{
    // select all sides
    stk::mesh::Selector selAllSides =
        this->metaDataRef().universal_part() & stk::mesh::selectUnion(parts);
    stk::mesh::BucketVector const& sideBuckets =
        this->bulkDataRef().get_buckets(this->metaDataRef().side_rank(),
                                        selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;
        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        // extract master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());

        // extract master element specifics
        const label numScsBip = meFC->numIntPoints_;

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            T* value =
                stk::mesh::field_data(this->stkFieldRef(), sideBucket, iSide);
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                for (label i = 0; i < N; i++)
                {
                    value[ip * N + i] = val[i];
                }
            }
        }
    }
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::initializer_list<T> val)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::initializer_list<T> val,
                                 const stk::mesh::Part& part)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::initializer_list<T> val,
                                 const stk::mesh::PartVector parts)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::vector<T> val)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::vector<T> val,
                                 const stk::mesh::Part& part)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::vector<T> val,
                                 const stk::mesh::PartVector parts)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::array<T, N> val)
{
    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::array<T, N> val,
                                 const stk::mesh::Part& part)
{
    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <class T, size_t N>
void sideField<T, N>::setToValue(std::array<T, N> val,
                                 const stk::mesh::PartVector parts)
{
    T value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <class T, size_t N>
void sideField<T, N>::interpolate(const nodeSideField<T, N>& nsf,
                                  label iZone,
                                  label iBoundary,
                                  bool shifted)
{
    errorMsg("Must not reach here");
}

#ifdef HAS_INTERFACE
template <class T, size_t N>
void sideField<T, N>::interpolate(const nodeSideField<T, N>& nsf,
                                  label iInterface,
                                  bool master,
                                  bool shifted)
{
    errorMsg("Must not reach here");
}

template <class T, size_t N>
void sideField<T, N>::transfer(label iInterface, bool reverse, bool shifted)
{
    errorMsg("Must not reach here");
}
#endif /* HAS_INTERFACE */

template <class T, size_t N>
sideField<T, N>& sideField<T, N>::prevTimeRef()
{
    return *prevTimeFieldPtr_.get();
}

template <class T, size_t N>
const sideField<T, N>& sideField<T, N>::prevTimeRef() const
{
    return *prevTimeFieldPtr_.get();
}

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const sideField<T, N>& field)
{
    if (messager::master())
    {
        os << std::endl << field.stkFieldRef().name() << std::endl;
    }

    for (label iProc = 0; iProc < messager::nProcs(); iProc++)
    {
        if (messager::myProcNo() == iProc)
        {
            if (messager::parallel())
            {
                std::cout << "Proc: " << iProc << std::endl;
            }
            std::cout << "{" << std::endl;

            // select all sides relevant to the field
            stk::mesh::Selector selAllSides =
                field.metaDataRef().universal_part() &
                stk::mesh::selectField(field.stkFieldRef());
            stk::mesh::BucketVector const& sideBuckets =
                field.bulkDataRef().get_buckets(field.metaDataRef().side_rank(),
                                                selAllSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;
                const stk::mesh::Bucket::size_type nSidesPerBucket =
                    sideBucket.size();

                // extract master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                // extract master element specifics
                const label numScsBip = meFC->numIntPoints_;

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    const T* value = stk::mesh::field_data(
                        field.stkFieldRef(), sideBucket, iSide);

                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        for (label i = 0; i < N; i++)
                        {
                            os << std::scientific << std::setprecision(14)
                               << value[N * ip + i] << " ";
                        }

                        os << std::endl;
                    }
                }
            }

            os << "}" << std::endl;
        }
        messager::barrier();
    }

    return os;
}

} // namespace accel
