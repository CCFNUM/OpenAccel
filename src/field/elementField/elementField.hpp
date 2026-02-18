// File : elementField.hpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <class T, size_t N>
elementField<T, N>::elementField(mesh* meshPtr,
                                 std::string name,
                                 unsigned numberOfStates)
    : field<T, N>(meshPtr, stk::topology::ELEMENT_RANK, name, numberOfStates),
      isInitialized_(this->meshPtr()->nZones(), false)
{
    // Put the stk field on interior mesh parts
    putFieldOnRegisteredParts_();

    // Populate previous time step fields if required: special constructors are
    // employed here
    if (numberOfStates > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<elementField<T, N>>(
            meshPtr, &this->stkFieldRef().field_of_state(stk::mesh::StateN));
    }
}

template <class T, size_t N>
elementField<T, N>::elementField(mesh* meshPtr,
                                 stk::mesh::Field<T>* stkFieldPtr)
    : field<T, N>(meshPtr, stkFieldPtr),
      isInitialized_(this->meshPtr()->nZones(), false)
{
    if (stkFieldPtr->number_of_states() - int(stkFieldPtr->state()) > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<elementField<T, N>>(
            meshPtr,
            &this->stkFieldRef().field_of_state(
                (enum stk::mesh::FieldState)((int)stkFieldPtr->state() + 1)));
    }
}

template <class T, size_t N>
void elementField<T, N>::putFieldOnRegisteredParts_()
{
    // apply zone loop
    for (zone* zonePtr : this->meshPtr()->zoneVector())
    {
        // Put the stk field on registered parts
        for (const stk::mesh::Part* part : zonePtr->interiorParts())
        {
            putFieldOnPart(*part);
        }
    }
}

template <class T, size_t N>
void elementField<T, N>::putFieldOnPart(const stk::mesh::Part& part)
{
    // Determine number of integration points in the element
    const stk::topology theTopo = part.topology();
    MasterElement* meSCS =
        MasterElementRepo::get_surface_master_element(theTopo);
    const label numScsIp = meSCS->numIntPoints_;

    // initial value
    T zeroVal[N * numScsIp];
    for (label i = 0; i < N * numScsIp; i++)
    {
        zeroVal[i] = static_cast<T>(0);
    }

    // Put the field on mesh
    stk::mesh::put_field_on_mesh(
        this->stkFieldRef(), part, N * numScsIp, &zeroVal[0]);
}

template <class T, size_t N>
void elementField<T, N>::setToValue_(const T* val)
{
    // select all elements relevant to the field
    stk::mesh::Selector selAllElements =
        this->metaDataRef().universal_part() &
        stk::mesh::selectField(this->stkFieldRef());
    stk::mesh::BucketVector const& elementBuckets =
        this->bulkDataRef().get_buckets(stk::topology::ELEMENT_RANK,
                                        selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label numScsIp = meSCS->numIntPoints_;

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            T* value = stk::mesh::field_data(
                this->stkFieldRef(), elementBucket, iElement);

            for (label ip = 0; ip < numScsIp; ++ip)
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
void elementField<T, N>::setToValue_(const T* val, const stk::mesh::Part& part)
{
    // select all elements
    stk::mesh::Selector selAllElements =
        this->metaDataRef().universal_part() & part;
    stk::mesh::BucketVector const& elementBuckets =
        this->bulkDataRef().get_buckets(stk::topology::ELEMENT_RANK,
                                        selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label numScsIp = meSCS->numIntPoints_;

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            T* value = stk::mesh::field_data(
                this->stkFieldRef(), elementBucket, iElement);

            for (label ip = 0; ip < numScsIp; ++ip)
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
void elementField<T, N>::setToValue_(const T* val,
                                     const stk::mesh::PartVector parts)
{
    // select all elements
    stk::mesh::Selector selAllElements =
        this->metaDataRef().universal_part() & stk::mesh::selectUnion(parts);
    stk::mesh::BucketVector const& elementBuckets =
        this->bulkDataRef().get_buckets(stk::topology::ELEMENT_RANK,
                                        selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type length = elementBucket.size();

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label numScsIp = meSCS->numIntPoints_;

        for (stk::mesh::Bucket::size_type iElement = 0; iElement < length;
             ++iElement)
        {
            T* value = stk::mesh::field_data(
                this->stkFieldRef(), elementBucket, iElement);

            for (label ip = 0; ip < numScsIp; ++ip)
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
void elementField<T, N>::setToValue(std::initializer_list<T> val)
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
void elementField<T, N>::setToValue(std::initializer_list<T> val,
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
void elementField<T, N>::setToValue(std::initializer_list<T> val,
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
void elementField<T, N>::setToValue(std::vector<T> val)
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
void elementField<T, N>::setToValue(std::vector<T> val,
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
void elementField<T, N>::setToValue(std::vector<T> val,
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
void elementField<T, N>::setToValue(std::array<T, N> val)
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
void elementField<T, N>::setToValue(std::array<T, N> val,
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
void elementField<T, N>::setToValue(std::array<T, N> val,
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
void elementField<T, N>::registerSideField(label iZone, label iBoundary)
{
    // Instantiate if not yet
    if (sideFieldPtr_ == nullptr)
    {
        sideFieldPtr_ = std::make_unique<sideField<T, N>>(
            this->meshPtr(), this->name() + "_side", 1);

        // copy some properties to side field
        sideFieldPtr_->setURF(this->urf());
    }

    // Put the side field on the corresponding boundary part
    for (auto* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        this->sideFieldRef().putFieldOnPart(*part);
    }
}

#ifdef HAS_INTERFACE
template <class T, size_t N>
void elementField<T, N>::registerSideFieldsForInterfaceSide(
    label iInterface,
    bool master,
    bool onlyIfNonoverlap)
{
    if (master)
    {
        const auto& interf = this->meshRef().interfaceRef(iInterface);

        if (onlyIfNonoverlap && !interf.masterInfoRef().hasNonoverlap_)
            return;

        // Instantiate if not yet
        if (sideFieldPtr_ == nullptr)
        {
            sideFieldPtr_ = std::make_unique<sideField<T, N>>(
                this->meshPtr(), this->name() + "_side", 1);

            // copy some properties to side field
            sideFieldPtr_->setURF(this->urf());
        }

        // Put the side field on the corresponding parts
        for (auto* part : interf.masterInfoRef().currentPartVec_)
        {
            this->sideFieldRef().putFieldOnPart(*part);
        }
    }
    else
    {
        const auto& interf = this->meshRef().interfaceRef(iInterface);

        if (onlyIfNonoverlap && !interf.slaveInfoRef().hasNonoverlap_)
            return;

        // Instantiate if not yet
        if (sideFieldPtr_ == nullptr)
        {
            sideFieldPtr_ = std::make_unique<sideField<T, N>>(
                this->meshPtr(), this->name() + "_side", 1);

            // copy some properties to side field
            sideFieldPtr_->setURF(this->urf());
        }

        // Put the side field on the corresponding parts
        for (auto* part : interf.slaveInfoRef().currentPartVec_)
        {
            this->sideFieldRef().putFieldOnPart(*part);
        }
    }
}
#endif /* HAS_INTERFACE */

template <class T, size_t N>
elementField<T, N>& elementField<T, N>::operator=(const elementField<T, N>& fld)
{
    // select all locally owned nodes relevant to the field
    stk::mesh::Selector selAllElements =
        this->metaDataRef().universal_part() &
        stk::mesh::selectField(fld.stkFieldRef());
    stk::mesh::BucketVector const& elementBuckets =
        this->bulkDataRef().get_buckets(stk::topology::ELEMENT_RANK,
                                        selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label numScsIp = meSCS->numIntPoints_;

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            T* value1 = stk::mesh::field_data(
                this->stkFieldRef(), elementBucket, iElement);
            T* value2 = stk::mesh::field_data(
                fld.stkFieldRef(), elementBucket, iElement);

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                for (label i = 0; i < N; i++)
                {
                    value1[ip * N + i] = value2[ip * N + i];
                }
            }
        }
    }
}

// Initialize

template <class T, size_t N>
void elementField<T, N>::initialize(label iZone, bool force)
{
    if (this->isInitialized(iZone) && !force)
        return;

    initializeField(iZone);

    initializeSideField(iZone);

    this->synchronizeGhostedEntities(iZone);

    if (sideFieldPtr_ != nullptr)
    {
        this->sideFieldRef().synchronizeGhostedEntities(iZone);
    }

    this->setIsInitialized(iZone);
}

template <class T, size_t N>
void elementField<T, N>::initializeField(label iZone)
{
}

template <class T, size_t N>
void elementField<T, N>::initializeSideField(label iZone)
{
    if (sideFieldPtr_ == nullptr)
    {
        return;
    }

    // ensure current zone is active
    assert((this->isZoneSet(iZone)));

    zone* zonePtr = this->meshPtr()->zonePtr(iZone);

#ifdef HAS_INTERFACE
    // Interfaces
    for (const interface* interf : zonePtr->interfacesRef())
    {
        if (interf->isInternal())
        {
            initializeInterfaceSideField(interf->masterInfoPtr());
            initializeInterfaceSideField(interf->slaveInfoPtr());
        }
        else
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(iZone);

            initializeInterfaceSideField(interfaceSideInfoPtr);
        }
    }
#endif /* HAS_INTERFACE */

    // Boundaries
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        // consider only if mDot is defined on the boundary
        if (!this->sideFieldRef().definedOn(
                zonePtr->boundaryPtr(iBoundary)->parts()))
            continue;

        initializeBoundarySideField(iZone, iBoundary);
    }
}

#ifdef HAS_INTERFACE
template <class T, size_t N>
void elementField<T, N>::initializeInterfaceSideField(
    const interfaceSideInfo* interfaceSideInfoPtr)
{
}
#endif /* HAS_INTERFACE */

template <class T, size_t N>
void elementField<T, N>::initializeBoundarySideField(label iZone,
                                                     label iBoundary)
{
}

// Update

template <class T, size_t N>
void elementField<T, N>::update(label iZone)
{
    updateField(iZone);

    updateSideFields(iZone);

    this->synchronizeGhostedEntities(iZone);

    if (sideFieldPtr_ != nullptr)
    {
        this->sideFieldRef().synchronizeGhostedEntities(iZone);
    }
}

template <class T, size_t N>
void elementField<T, N>::updateField(label iZone)
{
}

template <class T, size_t N>
void elementField<T, N>::updateSideFields(label iZone)
{
    if (sideFieldPtr_ == nullptr)
    {
        return;
    }

    // ensure current zone is active
    assert((this->isZoneSet(iZone)));

    zone* zonePtr = this->meshPtr()->zonePtr(iZone);

#ifdef HAS_INTERFACE
    // Interfaces
    for (const interface* interf : zonePtr->interfacesRef())
    {
        if (interf->isInternal())
        {
            updateInterfaceSideField(interf->index(), true);
            updateInterfaceSideField(interf->index(), false);
        }
        else
        {
            updateInterfaceSideField(interf->index(),
                                     interf->isMasterZone(iZone));
        }
    }
#endif /* HAS_INTERFACE */

    // Boundaries
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        updateBoundarySideField(iZone, iBoundary);
    }
}

#ifdef HAS_INTERFACE
template <class T, size_t N>
void elementField<T, N>::updateInterfaceSideField(label iInterface, bool master)
{
}
#endif /* HAS_INTERFACE */

template <class T, size_t N>
void elementField<T, N>::updateBoundarySideField(label iZone, label iBoundary)
{
}

// Access

template <class T, size_t N>
elementField<T, N>& elementField<T, N>::prevTimeRef()
{
    return *prevTimeFieldPtr_.get();
}

template <class T, size_t N>
const elementField<T, N>& elementField<T, N>::prevTimeRef() const
{
    return *prevTimeFieldPtr_.get();
}

template <class T, size_t N>
sideField<T, N>& elementField<T, N>::sideFieldRef()
{
    return *sideFieldPtr_.get();
}

template <class T, size_t N>
const sideField<T, N>& elementField<T, N>::sideFieldRef() const
{
    return *sideFieldPtr_.get();
}

template <class T, size_t N>
sideField<T, N>* elementField<T, N>::sideFieldPtr()
{
    return sideFieldPtr_.get();
}

template <class T, size_t N>
const sideField<T, N>* elementField<T, N>::sideFieldPtr() const
{
    return sideFieldPtr_.get();
}

// Out-of-line definitions

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const elementField<T, N>& field)
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
            stk::mesh::Selector selOwnedElements =
                field.metaDataRef().universal_part() &
                stk::mesh::selectField(field.stkFieldRef());
            stk::mesh::BucketVector const& elementBuckets =
                field.bulkDataRef().get_buckets(stk::topology::ELEMENT_RANK,
                                                selOwnedElements);
            for (stk::mesh::BucketVector::const_iterator ib =
                     elementBuckets.begin();
                 ib != elementBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& elementBucket = **ib;

                std::cout << "Bucket topology: " << elementBucket.topology()
                          << std::endl;

                const stk::mesh::Bucket::size_type nElementsPerBucket =
                    elementBucket.size();

                // extract master element
                MasterElement* meSCS =
                    MasterElementRepo::get_surface_master_element(
                        elementBucket.topology());

                // extract master element specifics
                const label numScsIp = meSCS->numIntPoints_;

                for (stk::mesh::Bucket::size_type iElement = 0;
                     iElement < nElementsPerBucket;
                     ++iElement)
                {
                    T* value = stk::mesh::field_data(
                        field.stkFieldRef(), elementBucket, iElement);

                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        for (label i = 0; i < N; i++)
                        {
                            std::cout << " " << value[ip * N + i];
                        }
                    }

                    std::cout << std::endl;
                }
            }
            os << "}" << std::endl;
        }
        messager::barrier();
    }
    return os;
}

} // namespace accel
