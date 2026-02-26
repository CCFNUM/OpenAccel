// File       : nodeSideField.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Template and inline implementations for node-side field
// containers.
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <class T, size_t N>
nodeSideField<T, N>::nodeSideField(mesh* meshPtr,
                                   std::string name,
                                   unsigned numberOfStates)
    : field<T, N>(meshPtr, stk::topology::NODE_RANK, name, numberOfStates)
{
    // Populate previous time step fields if required: special constructors are
    // employed here
    if (numberOfStates > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<nodeSideField<T, N>>(
            meshPtr, &this->stkFieldRef().field_of_state(stk::mesh::StateN));
    }
}

template <class T, size_t N>
nodeSideField<T, N>::nodeSideField(mesh* meshPtr,
                                   stk::mesh::Field<T>* stkFieldPtr)
    : field<T, N>(meshPtr, stkFieldPtr)
{
    if (stkFieldPtr->number_of_states() - int(stkFieldPtr->state()) > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<nodeSideField<T, N>>(
            meshPtr,
            &this->stkFieldRef().field_of_state(
                (enum stk::mesh::FieldState)((int)stkFieldPtr->state() + 1)));
    }
}

template <class T, size_t N>
void nodeSideField<T, N>::putFieldOnPart(const stk::mesh::Part& part)
{
    // initial value
    T zeroVal[N];
    for (label i = 0; i < N; i++)
    {
        zeroVal[i] = static_cast<T>(0);
    }

    stk::mesh::put_field_on_mesh(this->stkFieldRef(), part, N, &zeroVal[0]);
}

template <class T, size_t N>
bool nodeSideField<T, N>::definedOn(const stk::mesh::Part& part) const
{
    return (this->stkFieldRef().defined_on(part));
}

template <class T, size_t N>
bool nodeSideField<T, N>::definedOn(const stk::mesh::PartVector parts) const
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
bool nodeSideField<T, N>::definedOnBoundary(label iZone, label iBoundary) const
{
    if (!this->definedOn(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts()))
    {
        return false;
    }
    return true;
}

template <class T, size_t N>
void nodeSideField<T, N>::setToValue_(const T* val)
{
    // get field
    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes relevant to the field
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() &
        stk::mesh::selectField(this->stkFieldRef());
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();
        T* valueb = stk::mesh::field_data(stkFieldRef, sideNodeBucket);
        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * iSideNode + i] = val[i];
            }
        }
    }
}

template <class T, size_t N>
void nodeSideField<T, N>::setToValue_(const T* val, const stk::mesh::Part& part)
{
    if (!this->definedOn(part))
    {
        errorMsg("Can't set value: field " + this->name() +
                 " is not defined on " + part.name());
    }

    // get field
    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() & part;
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();
        T* valueb = stk::mesh::field_data(stkFieldRef, sideNodeBucket);
        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * iSideNode + i] = val[i];
            }
        }
    }
}

template <class T, size_t N>
void nodeSideField<T, N>::setToValue_(const T* val,
                                      const stk::mesh::PartVector parts)
{
    for (auto* part : parts)
    {
        if (!this->definedOn(*part))
        {
            errorMsg("Can't set value: field " + this->name() +
                     " is not defined on " + part->name());
        }
    }

    // get field
    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() & stk::mesh::selectUnion(parts);
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();
        T* valueb = stk::mesh::field_data(stkFieldRef, sideNodeBucket);
        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * iSideNode + i] = val[i];
            }
        }
    }
}

template <class T, size_t N>
void nodeSideField<T, N>::setToValue(std::initializer_list<T> val)
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
void nodeSideField<T, N>::setToValue(std::initializer_list<T> val,
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
void nodeSideField<T, N>::setToValue(std::initializer_list<T> val,
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
void nodeSideField<T, N>::setToValue(std::vector<T> val)
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
void nodeSideField<T, N>::setToValue(std::vector<T> val,
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
void nodeSideField<T, N>::setToValue(std::vector<T> val,
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
void nodeSideField<T, N>::setToValue(std::array<T, N> val)
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
void nodeSideField<T, N>::setToValue(std::array<T, N> val,
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
void nodeSideField<T, N>::setToValue(std::array<T, N> val,
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
void nodeSideField<T, N>::interpolate(const sideField<T, N>& sf,
                                      label iZone,
                                      label iBoundary)
{
    errorMsg("Must not reach here");
}

#ifdef HAS_INTERFACE
template <class T, size_t N>
void nodeSideField<T, N>::interpolate(const sideField<T, N>& sf,
                                      label iInterface,
                                      bool master)
{
    errorMsg("Must not reach here");
}
#endif /* HAS_INTERFACE */

template <class T, size_t N>
nodeSideField<T, N>& nodeSideField<T, N>::prevTimeRef()
{
    return *prevTimeFieldPtr_.get();
}

template <class T, size_t N>
const nodeSideField<T, N>& nodeSideField<T, N>::prevTimeRef() const
{
    return *prevTimeFieldPtr_.get();
}

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const nodeSideField<T, N>& field)
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

            // select all locally owned nodes relevant to the field
            stk::mesh::Selector selOwnedNodes =
                field.metaDataRef().locally_owned_part() &
                stk::mesh::selectField(field.stkFieldRef());
            stk::mesh::BucketVector const& sideNodeBuckets =
                field.bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                selOwnedNodes);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideNodeBuckets.begin();
                 ib != sideNodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideNodeBucket = **ib;
                const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                    sideNodeBucket.size();
                T* valueb =
                    stk::mesh::field_data(field.stkFieldRef(), sideNodeBucket);
                for (stk::mesh::Bucket::size_type iSideNode = 0;
                     iSideNode < nSideNodesPerBucket;
                     ++iSideNode)
                {
                    for (label i = 0; i < N; i++)
                    {
                        os << std::scientific << std::setprecision(14)
                           << valueb[N * iSideNode + i] << " ";
                    }

                    os << std::endl;
                }
            }
            os << "}" << std::endl;
        }
        messager::barrier();
    }

    return os;
}

} // namespace accel
