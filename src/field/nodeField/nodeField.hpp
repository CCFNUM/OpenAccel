// File : nodeField.hpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <size_t N, size_t M>
nodeField<N, M>::nodeField(mesh* meshPtr,
                           std::string name,
                           unsigned numberOfStates,
                           bool prevIter,
                           bool highResolution,
                           bool computeGradient,
                           bool correctedBoundaryNodeValues)
    : field<scalar, N>(meshPtr, stk::topology::NODE_RANK, name, numberOfStates),
      isInitialized_(this->meshPtr()->nZones(), false),
      correctedBoundaryNodeValues_(correctedBoundaryNodeValues)
{
    // Put the stk field on interior mesh parts
    putFieldOnRegisteredParts_();

    // Populate previous iteration field if required
    if (prevIter)
    {
        prevIterFieldPtr_ = std::make_unique<nodeField<N, M>>(
            meshPtr, name + "_prev_iter", 1, false);
    }

    // Populate previous time step fields if required: special constructors are
    // employed here
    if (numberOfStates > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<nodeField<N, M>>(
            meshPtr, &this->stkFieldRef().field_of_state(stk::mesh::StateN));
    }

    // Boundary conditions array

    // set size
    boundaryConditionsDictionaryArray_.resize(this->meshPtr()->nZones());

    for (label iZone = 0; iZone < this->meshPtr()->nZones(); iZone++)
    {
        // set size
        boundaryConditionsDictionaryArray_[iZone].resize(
            this->meshPtr()->zonePtr(iZone)->nBoundaries());

        // set index and name only for now. Value will be later filled when
        // required
        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            boundaryConditionsDictionaryArray_[iZone][iBoundary].setIndex(
                iBoundary);
            boundaryConditionsDictionaryArray_[iZone][iBoundary].setName(
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).name());
        }
    }

    // Initial conditions array
    initialConditionsDictionaryArray_.resize(this->meshPtr()->nZones());

    // Gradient required?
    if (computeGradient)
    {
        // limit the gradient? enable boolean if yes
        limitGradient_ = this->meshRef()
                             .controlsRef()
                             .solverRef()
                             .solverControl_.expertParameters_.limitGradients_;

        // correct the gradient at symmetry planes? enable boolean if yes
        correctGradient_ =
            this->meshRef()
                .controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.correctGradients_;

        // use incremental gradient changes for better accuracy? enable boolean
        // if yes
        incrementalGradientChange_ =
            this->meshRef()
                .controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.incrementalGradientChange_;

        // setup the gradient
        this->setupGradientField();
    }

    // High-resolution: check if gradient is enabled since required
    if (highResolution)
    {
        // Turn HR flag ON
        advectionScheme_ = advectionSchemeType::highResolution;

        // Update high-resolution blending coefficient field
        this->setupBlendingFactorField();
    }

    // If gradient limiter or a high-resolution is required, min/max fields must
    // be instantiated
    if (limitGradient_ ||
        advectionScheme_ == advectionSchemeType::highResolution)
    {
        this->setupMinMaxFields();
    }
}

template <size_t N, size_t M>
nodeField<N, M>::nodeField(mesh* meshPtr, STKScalarField* stkFieldPtr)
    : field<scalar, N>(meshPtr, stkFieldPtr),
      isInitialized_(this->meshPtr()->nZones(), false)
{
    if (stkFieldPtr->number_of_states() - int(stkFieldPtr->state()) > 1)
    {
        prevTimeFieldPtr_ = std::make_unique<nodeField<N, M>>(
            meshPtr,
            &this->stkFieldRef().field_of_state(
                (enum stk::mesh::FieldState)((int)stkFieldPtr->state() + 1)));
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setZone(label iZone)
{
    // enable zone for the current field
    field<scalar, N>::setZone(iZone);

    // enable zone for prev iter field
    if (prevIterFieldPtr_)
    {
        this->prevIterRef().setZone(iZone);
    }

    // enable zone for prev time field
    if (prevTimeFieldPtr_)
    {
        this->prevTimeRef().setZone(iZone);
    }

    // enable zone for gradient field
    if (gradFieldPtr_)
    {
        this->gradRef().setZone(iZone);
    }

    // enable zone for min field
    if (minValueFieldPtr_)
    {
        this->minValueRef().setZone(iZone);
    }

    // enable zone for max field
    if (maxValueFieldPtr_)
    {
        this->maxValueRef().setZone(iZone);
    }

    // enable zone for beta field
    if (blendingFactorFieldPtr_)
    {
        this->blendingFactorRef().setZone(iZone);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::putFieldOnRegisteredParts_()
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

template <size_t N, size_t M>
void nodeField<N, M>::putFieldOnPart(const stk::mesh::Part& part)
{
    // initial value
    scalar zeroVal[N];
    for (label i = 0; i < N; i++)
    {
        zeroVal[i] = 0;
    }

    stk::mesh::put_field_on_mesh(this->stkFieldRef(), part, N, &zeroVal[0]);
}

template <size_t N, size_t M>
bool nodeField<N, M>::definedOn(const stk::mesh::Part& part) const
{
    return (this->stkFieldRef().defined_on(part));
}

template <size_t N, size_t M>
bool nodeField<N, M>::definedOn(const stk::mesh::PartVector parts) const
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

// Operations

template <size_t N, size_t M>
void nodeField<N, M>::setToValue_(const scalar* val)
{
    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes relevant to the field
    stk::mesh::Selector selAllNodes = this->metaDataRef().universal_part() &
                                      stk::mesh::selectField(stkFieldRef);
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* valueb = stk::mesh::field_data(stkFieldRef, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * k + i] = val[i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue_(const scalar* val,
                                  const stk::mesh::Part& part)
{
    if (!this->definedOn(part))
    {
        errorMsg("Can't set value: field " + this->name() +
                 " is not defined on " + part.name());
    }

    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() & part;
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* valueb = stk::mesh::field_data(stkFieldRef, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * k + i] = val[i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue_(const scalar* val,
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

    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() & stk::mesh::selectUnion(parts);
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* valueb = stk::mesh::field_data(stkFieldRef, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * k + i] = val[i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::initializer_list<scalar> val)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::initializer_list<scalar> val,
                                 const stk::mesh::Part& part)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::initializer_list<scalar> val,
                                 const stk::mesh::PartVector parts)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::vector<scalar> val)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::vector<scalar> val,
                                 const stk::mesh::Part& part)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::vector<scalar> val,
                                 const stk::mesh::PartVector parts)
{
    if (N != val.size())
    {
        errorMsg("Setting a value of size " + std::to_string(val.size()) +
                 " to field of size " + std::to_string(N));
    }

    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::array<scalar, N> val)
{
    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0]);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::array<scalar, N> val,
                                 const stk::mesh::Part& part)
{
    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], part);
}

template <size_t N, size_t M>
void nodeField<N, M>::setToValue(std::array<scalar, N> val,
                                 const stk::mesh::PartVector parts)
{
    scalar value[N];
    label i = 0;
    for (auto cv : val)
    {
        value[i++] = cv;
    }

    this->setToValue_(&value[0], parts);
}

template <size_t N, size_t M>
void nodeField<N, M>::registerSideFields(label iZone, label iBoundary)
{
    registerNodeSideField(iZone, iBoundary);
    registerSideField(iZone, iBoundary);
}

template <size_t N, size_t M>
void nodeField<N, M>::registerNodeSideField(label iZone, label iBoundary)
{
    // Instantiate if not yet
    if (nodeSideFieldPtr_ == nullptr)
    {
        nodeSideFieldPtr_ = std::make_unique<nodeSideField<scalar, N>>(
            this->meshPtr(), this->name() + "_node_side", 1);
    }

    if (this->nodeSideFieldRef().isZoneUnset(iZone))
    {
        this->nodeSideFieldRef().setZone(iZone);
    }

    // Put the node side field on the corresponding boundary part
    for (auto* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        this->nodeSideFieldRef().putFieldOnPart(*part);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::registerSideField(label iZone, label iBoundary)
{
    // Instantiate if not yet
    if (sideFieldPtr_ == nullptr)
    {
        sideFieldPtr_ = std::make_unique<sideField<scalar, N>>(
            this->meshPtr(), this->name() + "_side", 1);
    }

    if (this->sideFieldRef().isZoneUnset(iZone))
    {
        this->sideFieldRef().setZone(iZone);
    }

    // Put the side field on the corresponding boundary part
    for (auto* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        this->sideFieldRef().putFieldOnPart(*part);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::registerSideFluxField(label iZone, label iBoundary)
{
    // Instantiate if not yet
    if (sideFluxFieldPtr_ == nullptr)
    {
        sideFluxFieldPtr_ = std::make_unique<sideField<scalar, N>>(
            this->meshPtr(), this->name() + "_side_flux", 1);
    }

    // Put the side field on the corresponding boundary part
    for (auto* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        this->sideFluxFieldRef().putFieldOnPart(*part);
    }
}

template <size_t N, size_t M>
nodeField<N, M>& nodeField<N, M>::operator=(const nodeField& fld)
{
    auto& stkFieldRef = this->stkFieldRef();

    // select all nodes relevant to the field
    stk::mesh::Selector selAllNodes = this->metaDataRef().universal_part() &
                                      stk::mesh::selectField(fld.stkFieldRef());
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* valueb = stk::mesh::field_data(stkFieldRef, b);
        const scalar* rhsValueb = stk::mesh::field_data(fld.stkFieldRef(), b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * k + i] = rhsValueb[N * k + i];
            }
        }
    }
    return (*this);
}

template <size_t N, size_t M>
void nodeField<N, M>::setupGradientField()
{
    assert(N * SPATIAL_DIM == M);
    gradFieldPtr_ = std::make_unique<nodeField<M>>(
        this->meshPtr(), this->name() + "_gradient", 1, false);

    if (limitGradient_)
    {
        this->setupGradientLimiterField();
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setupGradientLimiterField()
{
    // populate gradient limiter field if required
    gradientLimiterFieldPtr_ = std::make_unique<nodeField<N>>(
        this->meshPtr(), this->name() + "_gradient_limiter", 1, true);
    this->gradientLimiterRef().setToValue(std::vector<scalar>(N, 0));

    // Gradient field must be enabled already
    assert(gradFieldPtr_ != nullptr);
}

template <size_t N, size_t M>
void nodeField<N, M>::setupBlendingFactorField(bool dummy)
{
    // populate blending function field if required: if not required we still
    // populate with 0 values
    blendingFactorFieldPtr_ = std::make_unique<nodeField<N>>(
        this->meshPtr(), this->name() + "_blending_factor", 1, true);
    this->blendingFactorRef().setToValue(std::vector<scalar>(N, 0));

    if (!dummy)
    {
        assert(gradFieldPtr_ != nullptr);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::setupMinMaxFields()
{
    // auxiliary fields for high resolution / gradient limiter
    minValueFieldPtr_ = std::make_unique<nodeField<N, M>>(
        this->meshPtr(), this->name() + "_min_value", 1, false);
    maxValueFieldPtr_ = std::make_unique<nodeField<N, M>>(
        this->meshPtr(), this->name() + "_max_value", 1, false);
}

// Initialize

template <size_t N, size_t M>
void nodeField<N, M>::initialize(label iZone, bool force)
{
    if (this->isInitialized(iZone) && !force)
        return;

    this->initializeField(iZone);

    this->initializeSideFields(iZone);

    this->setIsInitialized(iZone);

    // supplementary tasks for transient simulations, primarily for second order
    // transient schemes
    this->updatePrevTimeField(iZone);
}

template <size_t N, size_t M>
void nodeField<N, M>::initializeField(label iZone)
{
    // ensure field is enabled on the zone
    assert(this->isZoneSet(iZone));

    bool initialize_field = true;
    const auto& restart_ctrl =
        this->meshRef().controlsRef().solverRef().restartControl_;

    // restore field data for explicit restarts
    if (restart_ctrl.isRestart_)
    {
        stk::mesh::FieldBase* theField = stk::mesh::get_field_by_name(
            this->name(), this->meshRef().metaDataRef());
        stk::io::MeshField mf(
            theField, theField->name(), restart_ctrl.timeMatchOption_);

        scalar restart_time = restart_ctrl.restartTime_;
        if (restart_time == 0.0)
        {
            restart_time = this->meshRef().ioBrokerRef().get_max_time();
        }
        // FIXME: [2025-01-14] The call to
        // mf.set_read_time() will always set
        // stk::io::MeshField::TimeMatchOption to SPECIFIED. This effectively
        // disables use of LINEAR_INTERPOLATION. STK API does not support this
        // as a per field operation. Instead, add_input_field() methods should
        // be called for a field followed by bulk read_defined_input_fields*
        // calls rather than individual field processing. Old code did this, the
        // new code does not.
        mf.set_read_time(restart_time);
        mf.set_single_state(false);

        // limit initialized parts to current zone's
        for (const stk::mesh::Part* part :
             this->meshPtr()->zonePtr(iZone)->interiorParts())
        {
            mf.add_subset(*part);
        }

        // restore
        this->meshRef().ioBrokerRef().read_input_field(mf);

        if (mf.field_restored())
        {
            initialize_field = false;
            if (messager::master())
            {
                std::cout << "Field " + this->name() +
                                 " has been restored at time "
                          << std::scientific << mf.time_restored() << std::endl;
            }

            // FIXME: [2025-01-14] Since we operate on fields individually
            // rather than calling bulk field restore method from STK API, it is
            // possible that individual fields may have inconsistent restored
            // time (typically this should not happen based on tests I did, but
            // it is a possibility and it will happen here).
            //
            // WARNING: The assignments to the global
            // controlsRef().time/globalIter variables below will be executed
            // for EVERY restored field, it is assumed these calls will be
            // directed to the correct internal STK database and hence they all
            // produce the same value for mf.time_restored() and global
            // parameter fetch.
            if (this->meshRef().controlsRef().isTransient())
            {
                this->meshRef().controlsRef().time = mf.time_restored();
            }
            else
            {
                assert(mf.time_restored() -
                           static_cast<label>(mf.time_restored()) ==
                       0.0);
                assert(this->meshRef().controlsRef().globalIter ==
                       static_cast<label>(mf.time_restored()));
                this->meshRef().controlsRef().globalIter =
                    static_cast<label>(mf.time_restored());
            }
        }
    }

    // apply zone loop
    if (initialize_field)
    {
        const auto& initCond = this->initialConditionRef(iZone);
        if (initCond.type() == initialConditionOption::value)
        {
            const auto& data = initCond.template data<N>(this->name());

            switch (data.type())
            {
                case inputDataType::constant:
                    {
                        std::array<scalar, N> value{};
                        std::copy_n(data.value(), N, value.begin());

                        this->setToValue(
                            value,
                            this->meshPtr()->zonePtr(iZone)->interiorParts());

                        if (messager::master())
                        {
                            std::cout << "Field " + this->name() +
                                             " is initialized with a uniform "
                                             "value over zone "
                                      << this->meshPtr()->zonePtr(iZone)->name()
                                      << std::endl;
                        }
                    }
                    break;

                case inputDataType::expression:
                    {
                        typedef exprtk::symbol_table<scalar> symbol_table_t;
                        typedef exprtk::expression<scalar> expression_t;
                        typedef exprtk::parser<scalar> parser_t;
                        std::vector<expression_t> expression_list;

                        symbol_table_t symbol_table;
                        symbol_table.add_constants();

                        // declare function vars
                        scalar t, x, y, z;
                        symbol_table.add_variable("t", t);
                        symbol_table.add_variable("x", x);
                        symbol_table.add_variable("y", y);
                        symbol_table.add_variable("z", z);

                        expression_t componentExpression;
                        componentExpression.register_symbol_table(symbol_table);

                        parser_t parser;

                        for (label i = 0; i < N; i++)
                        {
                            if (parser.compile(data.expression()[i],
                                               componentExpression))
                            {
                                expression_list.push_back(componentExpression);
                            }
                            else
                            {
                                errorMsg("Error in the expression provided "
                                         "for field " +
                                         this->name() + ": " +
                                         data.expression()[i]);
                            }
                        }

                        // set time value
                        t = this->meshRef().controlsRef().time;

                        // get field
                        auto& stkFieldRef = this->stkFieldRef();

                        // extract geometric fields
                        const STKScalarField* coordsSTKFieldPtr =
                            this->metaDataRef().template get_field<scalar>(
                                stk::topology::NODE_RANK,
                                this->getCoordinatesID_(iZone));

                        // select all nodes relevant to the node side field
                        stk::mesh::Selector selAllNodes =
                            this->metaDataRef().universal_part() &
                            stk::mesh::selectUnion(this->meshPtr()
                                                       ->zonePtr(iZone)
                                                       ->interiorParts());
                        stk::mesh::BucketVector const& nodeBuckets =
                            this->bulkDataRef().get_buckets(
                                stk::topology::NODE_RANK, selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 nodeBuckets.begin();
                             ib != nodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& nodeBucket = **ib;
                            const stk::mesh::Bucket::size_type nNodesPerBucket =
                                nodeBucket.size();
                            scalar* valueb =
                                stk::mesh::field_data(stkFieldRef, nodeBucket);
                            for (stk::mesh::Bucket::size_type iNode = 0;
                                 iNode < nNodesPerBucket;
                                 ++iNode)
                            {
                                const scalar* coords = stk::mesh::field_data(
                                    *coordsSTKFieldPtr, nodeBucket, iNode);

#if SPATIAL_DIM == 3
                                x = coords[0];
                                y = coords[1];
                                z = coords[2];
#elif SPATIAL_DIM == 2
                                x = coords[0];
                                y = coords[1];
#endif
                                for (label i = 0; i < N; i++)
                                {
                                    valueb[N * iNode + i] =
                                        expression_list[i].value();
                                }
                            }
                        }

                        if (messager::master())
                        {
                            std::cout
                                << "Field " + this->name() +
                                       " is initialized from formula over "
                                       "zone "
                                << this->meshPtr()->zonePtr(iZone)->name()
                                << std::endl;
                        }
                    }
                    break;

                case inputDataType::profileData:
                    {
                        errorMsg("profile data not provided yet");
                    }
                    break;
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::initializeSideFields(label iZone)
{
    this->updateSideFields(iZone);
}

// Update

template <size_t N, size_t M>
void nodeField<N, M>::update(label iZone)
{
    this->updateField(iZone);

    this->updateSideFields(iZone);
}

template <size_t N, size_t M>
void nodeField<N, M>::updateField(label iZone)
{
    // ensure field is enabled on the zone
    assert(this->isZoneSet(iZone));
}

template <size_t N, size_t M>
void nodeField<N, M>::updateSideFields(label iZone)
{
    // ensure field is enabled on the zone
    assert(this->isZoneSet(iZone));

    zone* zonePtr = this->meshPtr()->zonePtr(iZone);

    // Boundaries
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        updateBoundarySideField(iZone, iBoundary);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updateBoundarySideField(label iZone, label iBoundary)
{
    boundaryPhysicalType physicalType =
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();
    const auto& bcType = this->boundaryConditionRef(iZone, iBoundary).type();

    switch (physicalType)
    {
        case boundaryPhysicalType::symmetry:
            break;

        case boundaryPhysicalType::inlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedValue:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        break;

                    case boundaryConditionType::zeroGradient:
                        break;

                    default:
                        errorMsg("Boundary condition " + toString(bcType) +
                                 " not implemented for field " + this->name());
                        break;
                }
            }
            break;

        case boundaryPhysicalType::outlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedValue:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        break;

                    case boundaryConditionType::zeroGradient:
                        break;

                    default:
                        errorMsg("Boundary condition " + toString(bcType) +
                                 " not implemented for field " + this->name());
                        break;
                }
            }
            break;

        case boundaryPhysicalType::opening:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedValue:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        break;

                    case boundaryConditionType::zeroGradient:
                        break;

                    default:
                        errorMsg("Boundary condition " + toString(bcType) +
                                 " not implemented for field " + this->name());
                        break;
                }
            }
            break;

        case boundaryPhysicalType::wall:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedValue:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        break;

                    case boundaryConditionType::specifiedFlux:
                        updateSideFluxField(iZone, iBoundary);
                        break;

                    case boundaryConditionType::zeroGradient:
                        break;

                    case boundaryConditionType::mixed:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        updateSideFluxField(iZone, iBoundary);
                        break;

                    default:
                        errorMsg("Boundary condition " + toString(bcType) +
                                 " not implemented for field " + this->name());
                        break;
                }
            }
            break;

        default:
            errorMsg("Boundary type " + toString(physicalType) +
                     " not implemented");
            break;
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updateBoundarySideFieldSpecifiedValue(label iZone,
                                                            label iBoundary)
{
    assert(nodeSideFieldPtr_);
    assert(sideFieldPtr_);

    // Consider only if the side field is defined on the boundary.
    if ((!this->sideFieldRef().definedOn(
            this->meshPtr()->zonePtr(iZone)->boundaryPtr(iBoundary)->parts())))
        errorMsg("Side field not defined on boundary");

    // get data
    auto& bc = this->boundaryConditionRef(iZone, iBoundary);
    auto& data = bc.template data<N>("value");

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (data.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                std::array<scalar, N> inputValue;
                std::copy(data.value(), data.value() + N, inputValue.begin());

                // skip if already initialized: a constant value may only be
                // updated once at initialization. Correct field if required.
                if (!isInitialized(iZone))
                {
                    this->nodeSideFieldRef().setToValue(
                        inputValue,
                        this->meshPtr()
                            ->zonePtr(iZone)
                            ->boundaryRef(iBoundary)
                            .parts());
                    this->sideFieldRef().setToValue(inputValue,
                                                    this->meshPtr()
                                                        ->zonePtr(iZone)
                                                        ->boundaryRef(iBoundary)
                                                        .parts());
                }

                if (correctedBoundaryNodeValues_)
                {
                    this->correctBoundaryNodes(iZone, iBoundary);
                }
            }
            break;

        case inputDataType::timeTable:
            {
                std::array<scalar, N> inputValue =
                    data.interpolate(this->meshRef().controlsRef().time);

                this->nodeSideFieldRef().setToValue(inputValue,
                                                    this->meshPtr()
                                                        ->zonePtr(iZone)
                                                        ->boundaryRef(iBoundary)
                                                        .parts());
                this->sideFieldRef().setToValue(inputValue,
                                                this->meshPtr()
                                                    ->zonePtr(iZone)
                                                    ->boundaryRef(iBoundary)
                                                    .parts());

                if (correctedBoundaryNodeValues_)
                {
                    this->correctBoundaryNodes(iZone, iBoundary);
                }
            }
            break;

        case inputDataType::expression:
            {
                typedef exprtk::symbol_table<scalar> symbol_table_t;
                typedef exprtk::expression<scalar> expression_t;
                typedef exprtk::parser<scalar> parser_t;
                std::vector<expression_t> expression_list;

                symbol_table_t symbol_table;
                symbol_table.add_constants();

                // declare function vars
                scalar t, x, y, z;
                symbol_table.add_variable("t", t);
                symbol_table.add_variable("x", x);
                symbol_table.add_variable("y", y);
                symbol_table.add_variable("z", z);

                expression_t componentExpression;
                componentExpression.register_symbol_table(symbol_table);

                parser_t parser;

                for (label i = 0; i < N; i++)
                {
                    if (parser.compile(data.expression()[i],
                                       componentExpression))
                    {
                        expression_list.push_back(componentExpression);
                    }
                    else
                    {
                        errorMsg("Error in the expression provided for field " +
                                 this->name() + ": " + data.expression()[i]);
                    }
                }

                // set time value
                t = this->meshRef().controlsRef().time;

                auto& nodeSideSTKFieldRef =
                    this->nodeSideFieldRef().stkFieldRef();
                auto& stkFieldRef = this->stkFieldRef();

                const auto& coordsSTKFieldRef =
                    *metaData.template get_field<scalar>(
                        stk::topology::NODE_RANK,
                        this->getCoordinatesID_(iZone));

                // select all nodes relevant to the node side field
                stk::mesh::Selector selAllNodes =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(this->meshPtr()
                                               ->zonePtr(iZone)
                                               ->boundaryRef(iBoundary)
                                               .parts());
                stk::mesh::BucketVector const& sideNodeBuckets =
                    bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         sideNodeBuckets.begin();
                     ib != sideNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& sideNodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                        sideNodeBucket.size();
                    scalar* snvalueb = stk::mesh::field_data(
                        nodeSideSTKFieldRef, sideNodeBucket);
                    scalar* valueb =
                        stk::mesh::field_data(stkFieldRef, sideNodeBucket);
                    for (stk::mesh::Bucket::size_type iSideNode = 0;
                         iSideNode < nSideNodesPerBucket;
                         ++iSideNode)
                    {
                        const scalar* coords = stk::mesh::field_data(
                            coordsSTKFieldRef, sideNodeBucket, iSideNode);

#if SPATIAL_DIM == 3
                        x = coords[0];
                        y = coords[1];
                        z = coords[2];
#elif SPATIAL_DIM == 2
                        x = coords[0];
                        y = coords[1];
#endif
                        for (label i = 0; i < N; i++)
                        {
                            snvalueb[N * iSideNode + i] =
                                expression_list[i].value();
                        }

                        // override internal node values
                        if (correctedBoundaryNodeValues_)
                        {
                            for (label i = 0; i < N; i++)
                            {
                                valueb[N * iSideNode + i] =
                                    snvalueb[N * iSideNode + i];
                            }
                        }
                    }
                }

                // Update side field on the current boundary
                this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                 iZone,
                                                 iBoundary,
                                                 this->isShifted());
            }
            break;

        case inputDataType::profileData:
            {
                // skip if already initialized: a profile data may only be
                // updated once at initialization. Correct field if required.
                if (!isInitialized(iZone))
                {
                    // invoke scatter-to-surface functionality and fill the
                    // node-side field
                    utils::scatterToSurface sts(data.scatterPoints(),
                                                this->meshPtr()
                                                    ->zonePtr(iZone)
                                                    ->boundaryRef(iBoundary)
                                                    .parts(),
                                                data.distancePower(),
                                                data.donorCount());
                    sts.interpolateToField(
                        data.scatterValues(),
                        this->nodeSideFieldRef().stkFieldPtr());

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }

                if (correctedBoundaryNodeValues_)
                {
                    this->correctBoundaryNodes(iZone, iBoundary);
                }
            }
            break;
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updateSideFluxField(label iZone, label iBoundary)
{
    assert(sideFluxFieldPtr_);

    // Consider only if a side field is defined on the boundary.
    if ((!this->sideFluxFieldRef().definedOn(
            this->meshPtr()->zonePtr(iZone)->boundaryPtr(iBoundary)->parts())))
        errorMsg("Side flux field not defined on boundary");

    auto& bc = this->boundaryConditionRef(iZone, iBoundary);
    auto& data = bc.template data<N>("flux");

    // common
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    switch (data.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
        case inputDataType::timeTable:
            {
                std::array<scalar, N> flux;
                if (data.type() == inputDataType::constant)
                {
                    std::copy(data.value(), data.value() + N, flux.begin());
                }
                else
                {
                    flux = data.interpolate(this->meshRef().controlsRef().time);
                }

                this->sideFluxFieldRef().setToValue(flux,
                                                    this->meshPtr()
                                                        ->zonePtr(iZone)
                                                        ->boundaryRef(iBoundary)
                                                        .parts());
            }
            break;

        case inputDataType::expression:
            {
                typedef exprtk::symbol_table<scalar> symbol_table_t;
                typedef exprtk::expression<scalar> expression_t;
                typedef exprtk::parser<scalar> parser_t;
                std::vector<expression_t> expression_list;

                symbol_table_t symbol_table;
                symbol_table.add_constants();

                // declare function vars
                scalar t, x, y, z;
                symbol_table.add_variable("t", t);
                symbol_table.add_variable("x", x);
                symbol_table.add_variable("y", y);
                symbol_table.add_variable("z", z);

                expression_t componentExpression;
                componentExpression.register_symbol_table(symbol_table);

                parser_t parser;

                for (label i = 0; i < N; i++)
                {
                    if (parser.compile(data.expression()[i],
                                       componentExpression))
                    {
                        expression_list.push_back(componentExpression);
                    }
                    else
                    {
                        errorMsg("Error in the expression provided for field " +
                                 this->name() + ": " + data.expression()[i]);
                    }
                }

                // set time value
                t = this->meshRef().controlsRef().time;

                // nodal fields to gather
                std::vector<scalar> ws_coordinates;

                // master element
                std::vector<scalar> ws_face_shape_function;
                std::vector<scalar> ws_coordinate_face_shape_function;

                // ip values
                std::vector<scalar> coordBip(SPATIAL_DIM);

                // pointers for fast access
                scalar* p_coordBip = &coordBip[0];

                auto& sideFluxSTKFieldRef =
                    this->sideFluxFieldRef().stkFieldRef();

                const auto& coordsSTKFieldRef =
                    *metaData.template get_field<scalar>(
                        stk::topology::NODE_RANK,
                        this->getCoordinatesID_(iZone));

                // shifted ip?
                bool isShifted = this->isShifted();

                // select all nodes relevant to the node side field
                stk::mesh::Selector selAllSides =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(this->meshPtr()
                                               ->zonePtr(iZone)
                                               ->boundaryRef(iBoundary)
                                               .parts());
                stk::mesh::BucketVector const& sideBuckets =
                    bulkData.get_buckets(metaData.side_rank(), selAllSides);
                for (stk::mesh::BucketVector::const_iterator ib =
                         sideBuckets.begin();
                     ib != sideBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& sideBucket = **ib;
                    const stk::mesh::Bucket::size_type nSidesPerBucket =
                        sideBucket.size();

                    // face master element
                    MasterElement* meFC =
                        accel::MasterElementRepo::get_surface_master_element(
                            sideBucket.topology());
                    const label nodesPerSide =
                        sideBucket.topology().num_nodes();
                    const label numScsBip = meFC->numIntPoints_;

                    // algorithm related; element
                    ws_coordinates.resize(nodesPerSide * SPATIAL_DIM);
                    ws_face_shape_function.resize(numScsBip * nodesPerSide);
                    ws_coordinate_face_shape_function.resize(numScsBip *
                                                             nodesPerSide);

                    // pointers
                    scalar* p_coordinates = &ws_coordinates[0];
                    scalar* p_face_shape_function = &ws_face_shape_function[0];
                    scalar* p_coordinate_face_shape_function =
                        &ws_coordinate_face_shape_function[0];

                    // shape functions; boundary
                    if (isShifted)
                    {
                        meFC->shifted_shape_fcn(&p_face_shape_function[0]);
                    }
                    else
                    {
                        meFC->shape_fcn(&p_face_shape_function[0]);
                    }

                    // Always use trilinear (standard) shape functions for
                    // coordinates
                    meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

                    for (stk::mesh::Bucket::size_type iSide = 0;
                         iSide < nSidesPerBucket;
                         ++iSide)
                    {
                        stk::mesh::Entity side = sideBucket[iSide];

                        scalar* flux =
                            stk::mesh::field_data(sideFluxSTKFieldRef, side);

                        //======================================
                        // gather nodal data off of element
                        //======================================
                        stk::mesh::Entity const* sideNodeRels =
                            bulkData.begin_nodes(side);
                        label numNodes = bulkData.num_nodes(side);

                        // sanity check on num nodes
                        STK_ThrowAssert(numNodes == nodesPerSide);
                        for (label ni = 0; ni < numNodes; ++ni)
                        {
                            stk::mesh::Entity node = sideNodeRels[ni];

                            // gather vectors
                            const scalar* coords =
                                stk::mesh::field_data(coordsSTKFieldRef, node);
                            const label offset = ni * SPATIAL_DIM;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_coordinates[offset + j] = coords[j];
                            }
                        }

                        for (label ip = 0; ip < numScsBip; ++ip)
                        {
                            // zero out vector quantities; form aMag
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_coordBip[j] = 0.0;
                            }

                            // interpolate to bip
                            const label offSetSF_face = ip * nodesPerSide;
                            for (label ic = 0; ic < nodesPerSide; ++ic)
                            {
                                const scalar r =
                                    p_face_shape_function[offSetSF_face + ic];
                                const scalar r_coord =
                                    p_coordinate_face_shape_function
                                        [offSetSF_face + ic];

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_coordBip[j] +=
                                        r_coord *
                                        p_coordinates[ic * SPATIAL_DIM + j];
                                }
                            }

#if SPATIAL_DIM == 3
                            x = p_coordBip[0];
                            y = p_coordBip[1];
                            z = p_coordBip[2];
#elif SPATIAL_DIM == 2
                            x = p_coordBip[0];
                            y = p_coordBip[1];
#endif

                            for (label i = 0; i < N; i++)
                            {
                                flux[N * ip + i] = expression_list[i].value();
                            }
                        }
                    }
                }
            }
            break;

        case inputDataType::profileData:
            {
                errorMsg("profile data not provided yet");
            }
            break;
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updatePrevIterField(label iZone)
{
    assert(this->isZoneSet(iZone));

    // get fields
    const auto& currentSTKFieldRef = this->stkFieldRef();
    auto& prevIterSTKFieldRef = this->prevIterRef().stkFieldRef();

    // fill
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        const scalar* valueb = stk::mesh::field_data(currentSTKFieldRef, b);
        scalar* prevValueb = stk::mesh::field_data(prevIterSTKFieldRef, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                prevValueb[N * k + i] = valueb[N * k + i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updatePrevTimeField(label iZone)
{
    assert(this->isZoneSet(iZone));

    if (prevTimeFieldPtr_ != nullptr)
    {
        // update prev of prev
        this->prevTimeRef().updatePrevTimeField(iZone);

        // update prev

        // get fields
        const auto& currentSTKFieldRef = this->stkFieldRef();
        auto& prevTimeSTKFieldRef = this->prevTimeRef().stkFieldRef();

        // fill
        stk::mesh::Selector selAllNodes =
            this->metaDataRef().universal_part() &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());
        stk::mesh::BucketVector const& nodeBuckets =
            this->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                            selAllNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            const scalar* valueb = stk::mesh::field_data(currentSTKFieldRef, b);
            scalar* prevValueb = stk::mesh::field_data(prevTimeSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    prevValueb[N * k + i] = valueb[N * k + i];
                }
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::restorePrevIterField(label iZone)
{
    assert(this->isZoneSet(iZone));

    // get fields
    const auto& prevIterSTKFieldRef = this->prevIterRef().stkFieldRef();
    auto& currentSTKFieldRef = this->stkFieldRef();

    // fill
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());
    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();
        scalar* valueb = stk::mesh::field_data(currentSTKFieldRef, b);
        const scalar* prevValueb =
            stk::mesh::field_data(prevIterSTKFieldRef, b);
        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            for (label i = 0; i < N; i++)
            {
                valueb[N * k + i] = prevValueb[N * k + i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::restorePrevTimeField(label iZone)
{
    assert(this->isZoneSet(iZone));

    if (prevTimeFieldPtr_ != nullptr)
    {
        // get fields
        const auto& prevTimeSTKFieldRef = this->prevTimeRef().stkFieldRef();
        auto& currentSTKFieldRef = this->stkFieldRef();

        // fill
        stk::mesh::Selector selAllNodes =
            this->metaDataRef().universal_part() &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());
        stk::mesh::BucketVector const& nodeBuckets =
            this->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                            selAllNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            scalar* value = stk::mesh::field_data(currentSTKFieldRef, b);
            const scalar* prevValue =
                stk::mesh::field_data(prevTimeSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    value[N * k + i] = prevValue[N * k + i];
                }
            }
        }

        // re-store
        this->prevTimeRef().restorePrevTimeField(iZone);
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updateMinMaxFields(label iZone,
                                         bool applyExtremaExpansion,
                                         scalar extremaCoeff,
                                         scalar smoothLimiter)
{
    // ensure zone is active
    assert(this->isZoneSet(iZone));

    stk::mesh::MetaData& metaData = this->meshPtr_->metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshPtr_->bulkDataRef();

    // Get fields
    const auto& phiSTKFieldRef = this->stkFieldRef();
    auto& minValueSTKFieldRef = this->minValueRef().stkFieldRef();
    auto& maxValueSTKFieldRef = this->maxValueRef().stkFieldRef();

    // Set values of min/max fields equal to those of the field, or in case of
    // extrema expansion, start with absolute min and max of the field
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    if (applyExtremaExpansion)
    {
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            scalar* minValue = stk::mesh::field_data(minValueSTKFieldRef, b);
            scalar* maxValue = stk::mesh::field_data(maxValueSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    minValue[N * k + i] = this->maxAcceptedValue_;
                    maxValue[N * k + i] = this->minAcceptedValue_;
                }
            }
        }
    }
    else
    {
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            const scalar* value = stk::mesh::field_data(phiSTKFieldRef, b);
            scalar* minValue = stk::mesh::field_data(minValueSTKFieldRef, b);
            scalar* maxValue = stk::mesh::field_data(maxValueSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    minValue[N * k + i] = value[N * k + i];
                    maxValue[N * k + i] = value[N * k + i];
                }
            }
        }
    }

    // Interior ip only
    {
        // nodal fields to gather
        std::vector<scalar> ws_phi;

        stk::mesh::Selector selAllElements =
            metaData.universal_part() &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            // extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            // extract master element specifics
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;

            // algorithm related
            ws_phi.resize(nodesPerElement * N);

            scalar* p_phi = &ws_phi[0];

            const label* lrscv = meSCS->adjacentNodes();

            for (size_t iElement = 0; iElement < nElementsPerBucket; ++iElement)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);
                const label nNodesPerElement =
                    elementBucket.num_nodes(iElement);

                for (label ni = 0; ni < nNodesPerElement; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    const scalar* phi =
                        stk::mesh::field_data(phiSTKFieldRef, node);

                    // gather values
                    for (label i = 0; i < N; i++)
                    {
                        p_phi[N * ni + i] = phi[i];
                    }
                }

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    scalar* fldMinValL =
                        stk::mesh::field_data(minValueSTKFieldRef, nodeL);
                    scalar* fldMaxValL =
                        stk::mesh::field_data(maxValueSTKFieldRef, nodeL);

                    scalar* fldMinValR =
                        stk::mesh::field_data(minValueSTKFieldRef, nodeR);
                    scalar* fldMaxValR =
                        stk::mesh::field_data(maxValueSTKFieldRef, nodeR);

                    for (label i = 0; i < N; i++)
                    {
                        fldMinValL[i] =
                            std::min(fldMinValL[i], p_phi[N * ir + i]);
                        fldMaxValL[i] =
                            std::max(fldMaxValL[i], p_phi[N * ir + i]);

                        fldMinValR[i] =
                            std::min(fldMinValR[i], p_phi[N * il + i]);
                        fldMaxValR[i] =
                            std::max(fldMaxValR[i], p_phi[N * il + i]);
                    }
                }
            }
        }
    }

    if (applyExtremaExpansion && extremaCoeff > 0.0)
    {
        const scalar extrema =
            extremaCoeff * (this->maxAcceptedValue_ - this->minAcceptedValue_);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            scalar* minValue = stk::mesh::field_data(minValueSTKFieldRef, b);
            scalar* maxValue = stk::mesh::field_data(maxValueSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    const size_t offset = N * k + i;
                    maxValue[offset] = std::min(maxValue[offset] + extrema,
                                                this->maxAcceptedValue_);
                    minValue[offset] = std::max(minValue[offset] - extrema,
                                                this->minAcceptedValue_);
                }
            }
        }
    }

    if (applyExtremaExpansion && smoothLimiter > 0.0)
    {
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            const scalar* value = stk::mesh::field_data(phiSTKFieldRef, b);
            scalar* minValue = stk::mesh::field_data(minValueSTKFieldRef, b);
            scalar* maxValue = stk::mesh::field_data(maxValueSTKFieldRef, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                for (label i = 0; i < N; i++)
                {
                    const size_t offset = N * k + i;
                    const scalar blendedMax =
                        smoothLimiter * value[offset] +
                        (1.0 - smoothLimiter) * maxValue[offset];
                    const scalar blendedMin =
                        smoothLimiter * value[offset] +
                        (1.0 - smoothLimiter) * minValue[offset];

                    maxValue[offset] =
                        std::min(blendedMax, this->maxAcceptedValue_);
                    minValue[offset] =
                        std::max(blendedMin, this->minAcceptedValue_);
                }
            }
        }
    }

    // synchronize in case of parallel
    this->minValueRef().synchronizeGhostedEntities(iZone);
    this->maxValueRef().synchronizeGhostedEntities(iZone);
}

template <size_t N, size_t M>
void nodeField<N, M>::updateBlendingFactorField(label iZone)
{
    // ensure zone is active
    assert(this->isZoneSet(iZone));

    // this call assumes the following:
    // 1) all field nodes on the current zone are up-to-date, including all
    // ghosted nodes
    // 2) in the case of an existing interface on the current zone, and given
    // the fact that this interface connects another zone, then the other zone's
    // field nodes has to be update-to-date including the ghosted nodes
    // 3) all min/max fields' nodes on the current zone are up-to-date,
    // including all ghosted nodes

    if (advectionScheme_ != advectionSchemeType::highResolution)
        return;

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // update min/max value fields if not yet updated. They are already updated
    // in case limited gradient is required and no need to do it here
    if (!limitGradient_)
    {
        this->updateMinMaxFields(iZone);
    }

    // store current beta for the purpose of relaxing it
    this->blendingFactorRef().updatePrevIterField(iZone);

    // Set the upper limit for beta: this is decisive of how compressive the
    // field is to be
    this->blendingFactorRef().setToValue(
        std::vector<scalar>(N, blendingFactorMax_),
        this->meshPtr()->zonePtr(iZone)->interiorParts());

    //
    // Evaluate the blending factor
    //

    // get fields
    const STKScalarField* phiSTKFieldPtr = this->stkFieldPtr();
    const STKScalarField* gradPhiSTKFieldPtr = this->gradRef().stkFieldPtr();

    const STKScalarField* minValueSTKFieldPtr =
        this->minValueRef().stkFieldPtr();
    const STKScalarField* maxValueSTKFieldPtr =
        this->maxValueRef().stkFieldPtr();

    STKScalarField* blendingFactorSTKFieldPtr =
        this->blendingFactorRef().stkFieldPtr();

    // Interior ip
    {
        // temporary arrays
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scs_areav;

        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;

        // define some common selectors
        stk::mesh::Selector selAllElements =
            metaData.universal_part() &
            stk::mesh::selectField(this->stkFieldRef()) &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());

        // extract geometric fields
        STKScalarField* coordsSTKFieldPtr = metaData.template get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));

        // fixed-size arrays
        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        // get pointers
        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            // extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            ws_shape_function.resize(numScsIp * nodesPerElement);

            scalar* p_shape_function = &ws_shape_function[0];

            if (interpolationScheme_ == interpolationSchemeType::linearLinear)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            // set sizes
            ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);
            ws_phi.resize(nodesPerElement * N);
            ws_gradPhi.resize(nodesPerElement * N * SPATIAL_DIM);
            ws_min.resize(nodesPerElement * N);
            ws_max.resize(nodesPerElement * N);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_scs_areav.resize(numScsIp * SPATIAL_DIM);

            // pointers
            scalar* p_coordinate_shape_function =
                &ws_coordinate_shape_function[0];
            scalar* p_phi = &ws_phi[0];
            scalar* p_gradPhi = &ws_gradPhi[0];
            scalar* p_min = &ws_min[0];
            scalar* p_max = &ws_max[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_scs_areav = &ws_scs_areav[0];

            // Coordinate shape function - always use trilinear (standard) shape
            // functions
            meSCS->shape_fcn(&p_coordinate_shape_function[0]);

            const label* lrscv = meSCS->adjacentNodes();

            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElementsPerBucket;
                 ++iElement)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);
                const label nNodesPerElement =
                    elementBucket.num_nodes(iElement);

                for (label iNode = 0; iNode < nNodesPerElement; ++iNode)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    const scalar* phi =
                        stk::mesh::field_data(*phiSTKFieldPtr, node);
                    const scalar* min =
                        stk::mesh::field_data(*minValueSTKFieldPtr, node);
                    const scalar* max =
                        stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                    const scalar* gradPhi =
                        stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr, node);

                    // gather vectors

                    const label offset = iNode * N;
                    for (label i = 0; i < N; ++i)
                    {
                        p_phi[offset + i] = phi[i];
                        p_min[offset + i] = min[i];
                        p_max[offset + i] = max[i];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_gradPhi[SPATIAL_DIM * offset + SPATIAL_DIM * i +
                                      j] = gradPhi[SPATIAL_DIM * i + j];
                        }
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordinates[iNode * SPATIAL_DIM + i] = coords[i];
                    }
                }

                // compute geometry
                scalar scs_error = 0.0;
                meSCS->determinant(
                    1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    // interpolate to scs ip
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        p_coordIp[j] = 0.0;
                    }
                    const label offset = ip * nodesPerElement;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r_coord =
                            p_coordinate_shape_function[offset + ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_coordIp[j] +=
                                r_coord * p_coordinates[ic * SPATIAL_DIM + j];
                        }
                    }

                    // left
                    {
                        // position vector from left node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            dxj[j] = p_coordIp[j] -
                                     p_coordinates[il * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, nodeL);

                        for (label i = 0; i < N; i++)
                        {
                            scalar rgrad = 0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                rgrad += p_gradPhi[il * N * SPATIAL_DIM +
                                                   i * SPATIAL_DIM + j] *
                                         dxj[j];
                            }

                            // stabilisation
                            if (rgrad > 0)
                            {
                                rgrad = rgrad + SMALL;
                            }
                            else
                            {
                                rgrad = rgrad - SMALL;
                            }

                            // local vals
                            scalar betaLocal = std::max(
                                (p_max[N * il + i] - p_phi[N * il + i]) / rgrad,
                                (p_min[N * il + i] - p_phi[N * il + i]) /
                                    rgrad);
                            scalar beta12Local = std::max(
                                -(p_max[N * il + i] - p_phi[N * il + i]) /
                                    rgrad,
                                -(p_min[N * il + i] - p_phi[N * il + i]) /
                                    rgrad);

                            scalar y = std::min(betaLocal, beta12Local);

                            if (blendingFactorMax_ < 2.0 - SMALL)
                            {
                                y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                            }

                            if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                            {
                                beta[i] = std::min(beta[i], y);
                            }
                        }
                    }

                    // right
                    {
                        // position vector from left node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[ir * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, nodeR);

                        for (label i = 0; i < N; i++)
                        {
                            scalar rgrad = 0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                rgrad += p_gradPhi[ir * N * SPATIAL_DIM +
                                                   i * SPATIAL_DIM + j] *
                                         dxj[j];
                            }

                            // stabilisation
                            if (rgrad > 0)
                            {
                                rgrad = rgrad + SMALL;
                            }
                            else
                            {
                                rgrad = rgrad - SMALL;
                            }

                            // local vals
                            scalar betaLocal = std::max(
                                (p_max[N * ir + i] - p_phi[N * ir + i]) / rgrad,
                                (p_min[N * ir + i] - p_phi[N * ir + i]) /
                                    rgrad);
                            scalar beta12Local = std::max(
                                -(p_max[N * ir + i] - p_phi[N * ir + i]) /
                                    rgrad,
                                -(p_min[N * ir + i] - p_phi[N * ir + i]) /
                                    rgrad);

                            scalar y = std::min(betaLocal, beta12Local);

                            if (blendingFactorMax_ < 2.0 - SMALL)
                            {
                                y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                            }

                            if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                            {
                                beta[i] = std::min(beta[i], y);
                            }
                        }
                    }
                }
            }
        }
    }

    // Boundary ip
    {
        // get geometric fields
        STKScalarField* coordsSTKFieldPtr = metaData.template get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));

        // nodal fields to gather; gather everything other than what we
        // are assembling
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;

        // geometry related to populate
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;

        // fixed size arrays
        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        // get pointers
        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            // define vector of parent topos; should always be UNITY in
            // size
            std::vector<stk::topology> parentTopo;

            // setup for buckets; union parts and ask for locally owned
            stk::mesh::PartVector partVec =
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts();

            // define some common selectors
            stk::mesh::Selector selAllSides =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);
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
                    accel::MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                // extract master element specifics
                const label nodesPerSide = meFC->nodesPerElement_;
                const label numScsIp = meFC->numIntPoints_;
                const label* ipNodeMap = meFC->ipNodeMap();

                // algorithm related
                ws_phi.resize(nodesPerSide * N);
                ws_gradPhi.resize(nodesPerSide * N * SPATIAL_DIM);
                ws_min.resize(nodesPerSide * N);
                ws_max.resize(nodesPerSide * N);
                ws_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinate_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinates.resize(nodesPerSide * SPATIAL_DIM);

                // pointers
                scalar* p_phi = &ws_phi[0];
                scalar* p_gradPhi = &ws_gradPhi[0];
                scalar* p_min = &ws_min[0];
                scalar* p_max = &ws_max[0];
                scalar* p_coordinates = &ws_coordinates[0];
                scalar* p_shape_function = &ws_shape_function[0];
                scalar* p_coordinate_shape_function =
                    &ws_coordinate_shape_function[0];

                if (interpolationScheme_ ==
                    interpolationSchemeType::linearLinear)
                {
                    meFC->shifted_shape_fcn(&p_shape_function[0]);
                }
                else
                {
                    meFC->shape_fcn(&p_shape_function[0]);
                }

                // Coordinate shape function - always use trilinear (standard)
                // shape functions
                meFC->shape_fcn(&p_coordinate_shape_function[0]);

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    //===============================================
                    // gather nodal data; this is how we do it now..
                    //===============================================
                    stk::mesh::Entity const* sideNodeRels =
                        sideBucket.begin_nodes(iSide);
                    label numNodes = sideBucket.num_nodes(iSide);

                    // sanity check on num nodes
                    STK_ThrowAssert(numNodes == nodesPerSide);

                    for (label ni = 0; ni < numNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];

                        const scalar* phi =
                            stk::mesh::field_data(*phiSTKFieldPtr, node);
                        const scalar* min =
                            stk::mesh::field_data(*minValueSTKFieldPtr, node);
                        const scalar* max =
                            stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                        const scalar* gradPhi =
                            stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
                        const scalar* coords =
                            stk::mesh::field_data(*coordsSTKFieldPtr, node);

                        // gather vectors

                        const label offset = ni * N;
                        for (label i = 0; i < N; ++i)
                        {
                            p_phi[offset + i] = phi[i];
                            p_min[offset + i] = min[i];
                            p_max[offset + i] = max[i];

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_gradPhi[SPATIAL_DIM * offset +
                                          SPATIAL_DIM * i + j] =
                                    gradPhi[SPATIAL_DIM * i + j];
                            }
                        }

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                        }
                    }

                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        // nearest node
                        const label nn = ipNodeMap[ip];

                        stk::mesh::Entity node = sideNodeRels[nn];

                        // interpolate to scs ip
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            p_coordIp[j] = 0.0;
                        }
                        const label offset = ip * nodesPerSide;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r_coord =
                                p_coordinate_shape_function[offset + ic];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_coordIp[j] +=
                                    r_coord *
                                    p_coordinates[ic * SPATIAL_DIM + j];
                            }
                        }

                        // position vector from left node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[nn * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, node);

                        for (label i = 0; i < N; i++)
                        {
                            scalar rgrad = 0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                rgrad += p_gradPhi[nn * N * SPATIAL_DIM +
                                                   i * SPATIAL_DIM + j] *
                                         p_dxj[j];
                            }

                            // stabilisation
                            if (rgrad > 0)
                            {
                                rgrad = rgrad + SMALL;
                            }
                            else
                            {
                                rgrad = rgrad - SMALL;
                            }

                            // local vals
                            scalar betaLocal = std::max(
                                (p_max[N * nn + i] - p_phi[N * nn + i]) / rgrad,
                                (p_min[N * nn + i] - p_phi[N * nn + i]) /
                                    rgrad);
                            scalar beta12Local = std::max(
                                -(p_max[N * nn + i] - p_phi[N * nn + i]) /
                                    rgrad,
                                -(p_min[N * nn + i] - p_phi[N * nn + i]) /
                                    rgrad);

                            scalar y = std::min(betaLocal, beta12Local);

                            if (blendingFactorMax_ < 2.0 - SMALL)
                            {
                                y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                            }

                            if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                            {
                                beta[i] = std::min(beta[i], y);
                            }
                        }
                    }
                }
            }
        }
    }

    // Relax beta field
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().locally_owned_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());

    const auto& blendingFactorSTKFieldRefPrev =
        this->blendingFactorRef().prevIterRef().stkFieldRef();

    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        scalar* beta =
            stk::mesh::field_data(*blendingFactorSTKFieldPtr, nodeBucket);
        const scalar* betaPrev =
            stk::mesh::field_data(blendingFactorSTKFieldRefPrev, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            for (label i = 0; i < N; ++i)
            {
                beta[N * iNode + i] =
                    0.75 * betaPrev[N * iNode + i] + 0.25 * beta[N * iNode + i];
            }
        }
    }

    // synchronize in case of parallel
    this->blendingFactorRef().synchronizeGhostedEntities(iZone);
}

template <size_t N, size_t M>
void nodeField<N, M>::updateGradientField(label iZone)
{
    // ensure zone is active
    assert(this->isZoneSet(iZone));

    // this call assumes the following:
    // 1) all field nodes on the current zone are up-to-date, including all
    // ghosted nodes
    // 2) in the case of an existing interface on the current zone, and given
    // the fact that this interface connects another zone, then the other zone's
    // field nodes has to be update-to-date including the ghosted nodes

    const auto& bulkData = this->meshRef().bulkDataRef();
    const auto& metaData = this->meshRef().metaDataRef();

    // relax the existing gradient
    this->gradRef().relax(iZone, 1.0 - gradURF_);

    // use incremental gradient change?
    scalar incMult = incrementalGradientChange_ ? 1.0 : 0.0;

    // get fields
    const STKScalarField* phiSTKFieldPtr = this->stkFieldPtr();
    STKScalarField* gradPhiSTKFieldPtr = this->gradRef().stkFieldPtr();

    // Interior ip contribution
    {
        std::vector<scalar> ws_shape_function;

        // define some common selectors
        stk::mesh::Selector selAllElements =
            this->metaDataRef().universal_part() &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());

        // get geometric fields
        STKScalarField* coordinates = metaData.template get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));
        const STKScalarField* dualNodalVolumeSTKFieldPtr =
            metaData.template get_field<scalar>(
                stk::topology::NODE_RANK, this->getDualNodalVolumeID_(iZone));

        // fixed size containers
        std::vector<scalar> phiIp(N);

        // pointers to fixed size
        scalar* p_phiIp = &phiIp[0];

        // scratch fields
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_dualVolume;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scs_areav;

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            // extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;

            ws_shape_function.resize(numScsIp * nodesPerElement);
            scalar* p_shape_function = &ws_shape_function[0];

            if (this->interpolationScheme_ ==
                interpolationSchemeType::linearLinear)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            // set sizes
            ws_phi.resize(nodesPerElement * N);
            ws_dualVolume.resize(nodesPerElement);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_scs_areav.resize(numScsIp * SPATIAL_DIM);

            // pointers to ws
            scalar* p_phi = &ws_phi[0];
            scalar* p_dualVolume = &ws_dualVolume[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_scs_areav = &ws_scs_areav[0];

            const label* lrscv = meSCS->adjacentNodes();

            for (stk::mesh::Bucket::size_type k = 0; k < nElementsPerBucket;
                 ++k)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(k);
                const label numNodes = elementBucket.num_nodes(k);

                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    // pointers to real data
                    const scalar* phi =
                        stk::mesh::field_data(*phiSTKFieldPtr, node);
                    const scalar* coords =
                        stk::mesh::field_data(*coordinates, node);

                    // gather scalars
                    p_dualVolume[ni] = *stk::mesh::field_data(
                        *dualNodalVolumeSTKFieldPtr, node);

                    // gather vectors
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    }

                    // gather N-dim
                    for (label i = 0; i < N; ++i)
                    {
                        p_phi[ni * N + i] = phi[i];
                    }
                }

                // compute geometry
                scalar scs_error = 0.0;
                meSCS->determinant(
                    1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                // start ip loop
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    // pointer to fields to assemble
                    scalar* gradPhiL =
                        stk::mesh::field_data(*gradPhiSTKFieldPtr, nodeL);
                    scalar* gradPhiR =
                        stk::mesh::field_data(*gradPhiSTKFieldPtr, nodeR);

                    // interpolate to scs point; operate on saved off
                    // ws_field
                    for (label j = 0; j < N; ++j)
                    {
                        p_phiIp[j] = 0.0;
                    }
                    const label offset = ip * nodesPerElement;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        for (label j = 0; j < N; ++j)
                        {
                            p_phiIp[j] += p_shape_function[offset + ic] *
                                          p_phi[ic * N + j];
                        }
                    }

                    // left and right volume
                    scalar inv_volL = 1.0 / p_dualVolume[il];
                    scalar inv_volR = 1.0 / p_dualVolume[ir];

                    // assemble to il/ir
                    for (label i = 0; i < N; ++i)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            gradPhiL[i * N + j] +=
                                gradURF_ *
                                (p_phiIp[i] - incMult * p_phi[il * N + i]) *
                                p_scs_areav[ip * SPATIAL_DIM + j] * inv_volL;
                            gradPhiR[i * N + j] -=
                                gradURF_ *
                                (p_phiIp[i] - incMult * p_phi[ir * N + i]) *
                                p_scs_areav[ip * SPATIAL_DIM + j] * inv_volR;
                        }
                    }
                }
            }
        }
    }

    // Boundary ip contribution
    {
        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef =
            *metaData.template get_field<scalar>(
                metaData.side_rank(), this->getExposedAreaVectorID_(iZone));
        const STKScalarField* dualNodalVolumeSTKFieldPtr =
            metaData.template get_field<scalar>(
                stk::topology::NODE_RANK, this->getDualNodalVolumeID_(iZone));

        // nodal fields to gather; gather everything other than what we
        // are assembling
        std::vector<scalar> ws_phi;

        // geometry related to populate
        std::vector<scalar> ws_shape_function;

        // ip data
        std::vector<scalar> phiIp(N);

        // pointers ..
        scalar* p_phiIp = &phiIp[0];

        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            boundaryPhysicalType physicalType =
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();

            switch (physicalType)
            {
                case boundaryPhysicalType::symmetry:
                    break;

                default:
                    {
                        // setup for buckets; union parts and ask for locally
                        // owned
                        stk::mesh::PartVector partVec =
                            this->meshPtr()
                                ->zonePtr(iZone)
                                ->boundaryRef(iBoundary)
                                .parts();

                        // define some common selectors
                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(partVec);

                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
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
                            MasterElement* meFC = accel::MasterElementRepo::
                                get_surface_master_element(
                                    sideBucket.topology());

                            // extract master element specifics
                            const label nodesPerSide = meFC->nodesPerElement_;
                            const label numScsIp = meFC->numIntPoints_;
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // algorithm related
                            ws_phi.resize(nodesPerSide * N);
                            ws_shape_function.resize(numScsIp * nodesPerSide);

                            // pointers
                            scalar* p_phi = &ws_phi[0];
                            scalar* p_shape_function = &ws_shape_function[0];

                            if (this->interpolationScheme_ ==
                                interpolationSchemeType::linearLinear)
                            {
                                meFC->shifted_shape_fcn(&p_shape_function[0]);
                            }
                            else
                            {
                                meFC->shape_fcn(&p_shape_function[0]);
                            }

                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // face data
                                scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef,
                                    sideBucket,
                                    iSide);

                                //===============================================
                                // gather nodal data; this is how we do it now..
                                //===============================================
                                stk::mesh::Entity const* sideNodeRels =
                                    sideBucket.begin_nodes(iSide);
                                label numNodes = sideBucket.num_nodes(iSide);

                                // sanity check on num nodes
                                STK_ThrowAssert(numNodes == nodesPerSide);

                                for (label ni = 0; ni < numNodes; ++ni)
                                {
                                    stk::mesh::Entity node = sideNodeRels[ni];

                                    // pointers to real data
                                    scalar* phi = stk::mesh::field_data(
                                        *phiSTKFieldPtr, node);

                                    // gather vectors
                                    const label offset = ni * N;
                                    for (label j = 0; j < N; ++j)
                                    {
                                        p_phi[offset + j] = phi[j];
                                    }
                                }

                                // start assembly
                                for (label ip = 0; ip < numScsIp; ++ip)
                                {
                                    // nearest node
                                    const label nn = ipNodeMap[ip];

                                    stk::mesh::Entity node = sideNodeRels[nn];

                                    // pointer to fields to assemble
                                    scalar* gradPhi = stk::mesh::field_data(
                                        *gradPhiSTKFieldPtr, node);

                                    // suplemental
                                    scalar vol = *stk::mesh::field_data(
                                        *dualNodalVolumeSTKFieldPtr, node);

                                    // interpolate to scs point; operate on
                                    // saved off ws_field
                                    for (label j = 0; j < N; ++j)
                                    {
                                        p_phiIp[j] = 0.0;
                                    }

                                    const label offset = ip * nodesPerSide;
                                    for (label ic = 0; ic < nodesPerSide; ++ic)
                                    {
                                        const scalar r =
                                            p_shape_function[offset + ic];
                                        for (label j = 0; j < N; ++j)
                                        {
                                            p_phiIp[j] += r * p_phi[ic * N + j];
                                        }
                                    }

                                    // nearest node volume
                                    scalar inv_vol = 1.0 / vol;

                                    // assemble to nearest node
                                    for (label i = 0; i < N; ++i)
                                    {
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            gradPhi[i * N + j] +=
                                                gradURF_ *
                                                (p_phiIp[i] -
                                                 incMult * p_phi[nn * N + i]) *
                                                areaVec[ip * SPATIAL_DIM + j] *
                                                inv_vol;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;
            }
        }
    }

    if (limitGradient_)
    {
        limitGradientField_(iZone);
    }

    // correct for symmetry planes
    if (correctGradient_)
    {
        correctGradientField_(iZone);
    }

    // synchronize in case of parallel
    this->gradRef().synchronizeGhostedEntities(iZone);
}

template <size_t N, size_t M>
void nodeField<N, M>::limitGradientField_(label iZone)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // required for calculation of limiter
    this->updateMinMaxFields(iZone);

    // get fields
    STKScalarField* gradPhiSTKFieldPtr = this->gradRef().stkFieldPtr();
    STKScalarField* gradLimiterSTKFieldPtr =
        this->gradientLimiterRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtr = this->stkFieldPtr();
    const STKScalarField* minValueSTKFieldPtr =
        this->minValueRef().stkFieldPtr();
    const STKScalarField* maxValueSTKFieldPtr =
        this->maxValueRef().stkFieldPtr();

    // Extract coordinate field
    STKScalarField* coordsSTKFieldPtr = metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));

    // Set the upper limit for the limiter field
    this->gradientLimiterRef().setToValue(
        std::vector<scalar>(N, 1.0),
        this->meshPtr()->zonePtr(iZone)->interiorParts());

    // Process interior elements
    {
        // Declare all temporary arrays at the top
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scs_areav;

        stk::mesh::Selector selAllElements =
            this->metaDataRef().universal_part() &
            stk::mesh::selectField(this->stkFieldRef()) &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

        // fixed size arrays
        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        // get pointers
        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            // Extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;

            // resize arrays
            ws_phi.resize(nodesPerElement * N);
            ws_gradPhi.resize(nodesPerElement * N * SPATIAL_DIM);
            ws_min.resize(nodesPerElement * N);
            ws_max.resize(nodesPerElement * N);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_shape_function.resize(numScsIp * nodesPerElement);
            ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);
            ws_scs_areav.resize(numScsIp * SPATIAL_DIM);

            // pointers
            scalar* p_phi = &ws_phi[0];
            scalar* p_gradPhi = &ws_gradPhi[0];
            scalar* p_min = &ws_min[0];
            scalar* p_max = &ws_max[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_shape_function = &ws_shape_function[0];
            scalar* p_coordinate_shape_function =
                &ws_coordinate_shape_function[0];
            scalar* p_scs_areav = &ws_scs_areav[0];

            if (interpolationScheme_ == interpolationSchemeType::linearLinear)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            // Coordinate shape function - always use trilinear (standard) shape
            // functions
            meSCS->shape_fcn(&p_coordinate_shape_function[0]);

            const label* lrscv = meSCS->adjacentNodes();

            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElementsPerBucket;
                 ++iElement)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);
                const label nNodesPerElement =
                    elementBucket.num_nodes(iElement);

                // Gather nodal data
                for (label iNode = 0; iNode < nNodesPerElement; ++iNode)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    // pointers to real fields
                    const scalar* phi =
                        stk::mesh::field_data(*phiSTKFieldPtr, node);
                    const scalar* min =
                        stk::mesh::field_data(*minValueSTKFieldPtr, node);
                    const scalar* max =
                        stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                    const scalar* gradPhi =
                        stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr, node);

                    // gather vectors
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordinates[iNode * SPATIAL_DIM + i] = coords[i];
                    }

                    // gather N-dim
                    const label offset = iNode * N;
                    for (label i = 0; i < N; ++i)
                    {
                        p_phi[offset + i] = phi[i];
                        p_min[offset + i] = min[i];
                        p_max[offset + i] = max[i];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_gradPhi[SPATIAL_DIM * offset + SPATIAL_DIM * i +
                                      j] = gradPhi[SPATIAL_DIM * i + j];
                        }
                    }
                }

                // Compute geometry
                scalar scs_error = 0.0;
                meSCS->determinant(
                    1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

                // Loop over integration points
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // Left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    // Interpolate coordinates to integration point
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        p_coordIp[j] = 0.0;
                    }
                    const label offset = ip * nodesPerElement;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r_coord =
                            p_coordinate_shape_function[offset + ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_coordIp[j] +=
                                r_coord * p_coordinates[ic * SPATIAL_DIM + j];
                        }
                    }

                    // Process left node
                    {
                        // Distance vector from left node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[il * SPATIAL_DIM + j];
                        }

                        scalar* limiter = stk::mesh::field_data(
                            *gradLimiterSTKFieldPtr, nodeL);

                        for (label i = 0; i < N; i++)
                        {
                            // Compute gradient projection: gradi·dx
                            scalar gradiDx = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                gradiDx += p_gradPhi[il * N * SPATIAL_DIM +
                                                     i * SPATIAL_DIM + j] *
                                           p_dxj[j];
                            }

                            scalar deltaMax =
                                p_max[il * N + i] - p_phi[il * N + i];
                            scalar deltaMin =
                                p_min[il * N + i] - p_phi[il * N + i];

                            if (gradiDx > deltaMax + SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMax / gradiDx);
                            }
                            else if (gradiDx < deltaMin - SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMin / gradiDx);
                            }
                        }
                    }

                    // Process right node
                    {
                        // Distance vector from right node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[ir * SPATIAL_DIM + j];
                        }

                        scalar* limiter = stk::mesh::field_data(
                            *gradLimiterSTKFieldPtr, nodeR);

                        for (label i = 0; i < N; i++)
                        {
                            // Compute gradient projection: gradi·dx
                            scalar gradiDx = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                gradiDx += p_gradPhi[ir * N * SPATIAL_DIM +
                                                     i * SPATIAL_DIM + j] *
                                           p_dxj[j];
                            }

                            scalar deltaMax =
                                p_max[ir * N + i] - p_phi[ir * N + i];
                            scalar deltaMin =
                                p_min[ir * N + i] - p_phi[ir * N + i];

                            if (gradiDx > deltaMax + SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMax / gradiDx);
                            }
                            else if (gradiDx < deltaMin - SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMin / gradiDx);
                            }
                        }
                    }
                }
            }
        }
    }

    // Process boundary sides
    {
        // Declare temporary arrays at the top
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;

        // fixed size arrays
        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        // get pointers
        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            // Setup for buckets; union parts and ask for locally owned
            stk::mesh::PartVector partVec =
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts();

            // Define some common selectors
            stk::mesh::Selector selAllSides =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);

            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                const stk::mesh::Bucket::size_type nSidesPerBucket =
                    sideBucket.size();

                // Extract master element
                MasterElement* meFC =
                    accel::MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                // Extract master element specifics
                const label nodesPerSide = meFC->nodesPerElement_;
                const label numScsIp = meFC->numIntPoints_;
                const label* ipNodeMap = meFC->ipNodeMap();

                // Resize arrays
                ws_phi.resize(nodesPerSide * N);
                ws_gradPhi.resize(nodesPerSide * N * SPATIAL_DIM);
                ws_min.resize(nodesPerSide * N);
                ws_max.resize(nodesPerSide * N);
                ws_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinate_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinates.resize(nodesPerSide * SPATIAL_DIM);

                // pointers
                scalar* p_phi = &ws_phi[0];
                scalar* p_gradPhi = &ws_gradPhi[0];
                scalar* p_min = &ws_min[0];
                scalar* p_max = &ws_max[0];
                scalar* p_coordinates = &ws_coordinates[0];
                scalar* p_shape_function = &ws_shape_function[0];
                scalar* p_coordinate_shape_function =
                    &ws_coordinate_shape_function[0];

                if (interpolationScheme_ ==
                    interpolationSchemeType::linearLinear)
                {
                    meFC->shifted_shape_fcn(&p_shape_function[0]);
                }
                else
                {
                    meFC->shape_fcn(&p_shape_function[0]);
                }

                // Coordinate shape function - always use trilinear (standard)
                // shape functions
                meFC->shape_fcn(&p_coordinate_shape_function[0]);

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    // Gather nodal data
                    stk::mesh::Entity const* sideNodeRels =
                        sideBucket.begin_nodes(iSide);
                    label numNodes = sideBucket.num_nodes(iSide);

                    // Sanity check on num nodes
                    STK_ThrowAssert(numNodes == nodesPerSide);

                    for (label ni = 0; ni < numNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];

                        // pointers to real data
                        const scalar* coords =
                            stk::mesh::field_data(*coordsSTKFieldPtr, node);
                        const scalar* phi =
                            stk::mesh::field_data(*phiSTKFieldPtr, node);
                        const scalar* min =
                            stk::mesh::field_data(*minValueSTKFieldPtr, node);
                        const scalar* max =
                            stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                        const scalar* gradPhi =
                            stk::mesh::field_data(*gradPhiSTKFieldPtr, node);

                        // gather vectors
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                        }

                        // gather N-dim fields
                        for (label i = 0; i < N; ++i)
                        {
                            p_phi[ni * N + i] = phi[i];
                            p_min[ni * N + i] = min[i];
                            p_max[ni * N + i] = max[i];

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_gradPhi[ni * N * SPATIAL_DIM +
                                          i * SPATIAL_DIM + j] =
                                    gradPhi[i * SPATIAL_DIM + j];
                            }
                        }
                    }

                    // Loop over integration points
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        // Nearest node
                        const label nn = ipNodeMap[ip];

                        stk::mesh::Entity node = sideNodeRels[nn];

                        // Interpolate coordinates to integration point
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            p_coordIp[j] = 0.0;
                        }
                        const label offset = ip * nodesPerSide;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r_coord =
                                p_coordinate_shape_function[offset + ic];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_coordIp[j] +=
                                    r_coord *
                                    p_coordinates[ic * SPATIAL_DIM + j];
                            }
                        }

                        // Distance vector from node to ip
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[nn * SPATIAL_DIM + j];
                        }

                        scalar* limiter = stk::mesh::field_data(
                            *gradLimiterSTKFieldPtr, node);

                        for (label i = 0; i < N; i++)
                        {
                            // Compute gradient projection: gradi·dx
                            scalar gradiDx = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                gradiDx += p_gradPhi[nn * N * SPATIAL_DIM +
                                                     i * SPATIAL_DIM + j] *
                                           p_dxj[j];
                            }

                            scalar deltaMax =
                                p_max[nn * N + i] - p_phi[nn * N + i];
                            scalar deltaMin =
                                p_min[nn * N + i] - p_phi[nn * N + i];

                            if (gradiDx > deltaMax + SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMax / gradiDx);
                            }
                            else if (gradiDx < deltaMin - SMALL)
                            {
                                limiter[i] =
                                    std::min(limiter[i], deltaMin / gradiDx);
                            }
                        }
                    }
                }
            }
        }
    }

    // Apply limiter
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().locally_owned_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        const scalar* limiterb =
            stk::mesh::field_data(*gradLimiterSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = nodeBucket[iNode];
            scalar* gradPhi = stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
            for (label i = 0; i < N; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    gradPhi[i * N + j] *= limiterb[N * iNode + i];
                }
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::updateScale()
{
    scalar maxValue = this->max();
    scalar minValue = this->min();

    if (messager::parallel())
    {
        messager::maxReduce(maxValue);
        messager::minReduce(minValue);
    }

    scalar scaleValue = std::max(maxValue, 0.0) - std::min(minValue, 0.0);
    if (scaleValue < 1.0e-9)
    {
        scaleValue = 1.0;
    }
    scaleValue = std::max(scaleValue, maxValue);

    this->setScale(scaleValue);
    this->setMaxScale(maxValue);
    this->setMinScale(minValue);
}

template <size_t N, size_t M>
void nodeField<N, M>::synchronize(label iZone)
{
    assert(this->isZoneSet(iZone));

    this->synchronizeGhostedEntities(iZone);

    if (nodeSideFieldPtr_ != nullptr)
    {
        this->nodeSideFieldRef().synchronizeGhostedEntities(iZone);
    }

    if (sideFieldPtr_ != nullptr)
    {
        this->sideFieldRef().synchronizeGhostedEntities(iZone);
    }

    if (sideFluxFieldPtr_ != nullptr)
    {
        this->sideFluxFieldRef().synchronizeGhostedEntities(iZone);
    }
}

template <size_t N, size_t M>
scalar nodeField<N, M>::max()
{
    scalar maxValue = -BIG;

    for (label iZone = 0; iZone < this->meshPtr()->nZones(); iZone++)
    {
        if (this->isZoneSet(iZone))
        {
            auto& stkFieldRef = this->stkFieldRef();

            // Interior
            stk::mesh::Selector selAllNodes =
                this->metaDataRef().universal_part() &
                stk::mesh::selectField(this->stkFieldRef()) &
                stk::mesh::selectUnion(
                    this->meshPtr()->zonePtr(iZone)->interiorParts());

            stk::mesh::BucketVector const& nodeBuckets =
                this->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                selAllNodes);
            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;
                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();
                const scalar* valueb =
                    stk::mesh::field_data(stkFieldRef, nodeBucket);
                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    if constexpr (N == 1)
                    {
                        // Signed scalar field: use the actual value
                        maxValue = std::max(maxValue, valueb[iNode]);
                    }
                    else
                    {
                        // vector/tensor-like: compare by magnitude
                        scalar magValue = 0.0;
                        for (label i = 0; i < N; i++)
                        {
                            magValue +=
                                valueb[N * iNode + i] * valueb[N * iNode + i];
                        }
                        magValue = std::sqrt(magValue);

                        // store max value
                        maxValue = std::max(maxValue, magValue);
                    }
                }
            }
        }
    }

    return maxValue;
}

template <size_t N, size_t M>
scalar nodeField<N, M>::min()
{
    scalar minValue = BIG;

    for (label iZone = 0; iZone < this->meshPtr()->nZones(); iZone++)
    {
        if (this->isZoneSet(iZone))
        {
            auto& stkFieldRef = this->stkFieldRef();

            // Interior
            stk::mesh::Selector selAllNodes =
                this->metaDataRef().universal_part() &
                stk::mesh::selectField(this->stkFieldRef()) &
                stk::mesh::selectUnion(
                    this->meshPtr()->zonePtr(iZone)->interiorParts());

            stk::mesh::BucketVector const& nodeBuckets =
                this->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                selAllNodes);
            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;
                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();
                const scalar* valueb =
                    stk::mesh::field_data(stkFieldRef, nodeBucket);
                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    if constexpr (N == 1)
                    {
                        // Signed scalar field: use the actual value
                        minValue = std::min(minValue, valueb[iNode]);
                    }
                    else
                    {
                        // vector/tensor-like: compare by magnitude
                        scalar magValue = 0.0;
                        for (label i = 0; i < N; i++)
                        {
                            magValue +=
                                valueb[N * iNode + i] * valueb[N * iNode + i];
                        }
                        magValue = std::sqrt(magValue);

                        // store min value
                        minValue = std::min(minValue, magValue);
                    }
                }
            }
        }
    }

    return minValue;
}

template <size_t N, size_t M>
void nodeField<N, M>::correctBoundaryNodes(label iZone, label iBoundary)
{
    assert(nodeSideFieldPtr_);
    assert(correctedBoundaryNodeValues_);

    // Get fields
    auto& stkFieldRef = this->stkFieldRef();
    const auto& nodeSideSTKFieldRef = this->nodeSideFieldRef().stkFieldRef();

    // define some common selectors; select All nodes
    stk::mesh::Selector selAllSideNodes =
        this->metaDataRef().universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts());

    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                        selAllSideNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket =
            sideNodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = sideNodeBucket[iNode];

            scalar* value = stk::mesh::field_data(stkFieldRef, node);
            const scalar* snvalue =
                stk::mesh::field_data(nodeSideSTKFieldRef, node);

            for (label i = 0; i < N; i++)
            {
                value[i] = snvalue[i];
            }
        }
    }
}

template <size_t N, size_t M>
void nodeField<N, M>::relax(label iZone, const scalar urf)
{
    // Get fields
    auto& stkFieldRef = this->stkFieldRef();

    // define some common selectors; select All nodes
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().universal_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = nodeBucket[iNode];

            scalar* value = stk::mesh::field_data(stkFieldRef, node);

            for (label i = 0; i < N; i++)
            {
                value[i] *= urf;
            }
        }
    }
}

// Access

template <size_t N, size_t M>
nodeField<N, M>& nodeField<N, M>::prevIterRef()
{
    return *prevIterFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N, M>& nodeField<N, M>::prevIterRef() const
{
    return *prevIterFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N, M>& nodeField<N, M>::prevTimeRef()
{
    return *prevTimeFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N, M>& nodeField<N, M>::prevTimeRef() const
{
    return *prevTimeFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N>& nodeField<N, M>::blendingFactorRef()
{
    return *blendingFactorFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N>& nodeField<N, M>::blendingFactorRef() const
{
    return *blendingFactorFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N>* nodeField<N, M>::blendingFactorPtr()
{
    return blendingFactorFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N>* nodeField<N, M>::blendingFactorPtr() const
{
    return blendingFactorFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N>& nodeField<N, M>::gradientLimiterRef()
{
    return *gradientLimiterFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N>& nodeField<N, M>::gradientLimiterRef() const
{
    return *gradientLimiterFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N>* nodeField<N, M>::gradientLimiterPtr()
{
    return gradientLimiterFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N>* nodeField<N, M>::gradientLimiterPtr() const
{
    return gradientLimiterFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N, M>& nodeField<N, M>::maxValueRef()
{
    return *maxValueFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N, M>& nodeField<N, M>::maxValueRef() const
{
    return *maxValueFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<N, M>& nodeField<N, M>::minValueRef()
{
    return *minValueFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<N, M>& nodeField<N, M>::minValueRef() const
{
    return *minValueFieldPtr_.get();
}

template <size_t N, size_t M>
nodeField<M>& nodeField<N, M>::gradRef()
{
    return *gradFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeField<M>& nodeField<N, M>::gradRef() const
{
    return *gradFieldPtr_.get();
}

template <size_t N, size_t M>
nodeSideField<scalar, N>& nodeField<N, M>::nodeSideFieldRef()
{
    return *nodeSideFieldPtr_.get();
}

template <size_t N, size_t M>
const nodeSideField<scalar, N>& nodeField<N, M>::nodeSideFieldRef() const
{
    return *nodeSideFieldPtr_.get();
}

template <size_t N, size_t M>
sideField<scalar, N>& nodeField<N, M>::sideFieldRef()
{
    return *sideFieldPtr_.get();
}

template <size_t N, size_t M>
const sideField<scalar, N>& nodeField<N, M>::sideFieldRef() const
{
    return *sideFieldPtr_.get();
}

template <size_t N, size_t M>
sideField<scalar, N>& nodeField<N, M>::sideFluxFieldRef()
{
    return *sideFluxFieldPtr_.get();
}

template <size_t N, size_t M>
const sideField<scalar, N>& nodeField<N, M>::sideFluxFieldRef() const
{
    return *sideFluxFieldPtr_.get();
}

// Out-of-line definitions

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& os, const nodeField<N, M>& field)
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

            for (label iZone = 0; iZone < field.meshPtr()->nZones(); iZone++)
            {
                if (field.isZoneSet(iZone))
                {
                    // select all locally owned nodes relevant to the field
                    stk::mesh::Selector selAllNodes =
                        field.metaDataRef().universal_part() &
                        stk::mesh::selectField(field.stkFieldRef()) &
                        stk::mesh::selectUnion(
                            field.meshPtr()->zonePtr(iZone)->interiorParts());

                    stk::mesh::BucketVector const& nodeBuckets =
                        field.bulkDataRef().get_buckets(
                            stk::topology::NODE_RANK, selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             nodeBuckets.begin();
                         ib != nodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& nodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nNodesPerBucket =
                            nodeBucket.size();
                        scalar* valueb = stk::mesh::field_data(
                            field.stkFieldRef(), nodeBucket);
                        for (stk::mesh::Bucket::size_type iNode = 0;
                             iNode < nNodesPerBucket;
                             ++iNode)
                        {
                            for (label i = 0; i < N; i++)
                            {
                                os << std::scientific << std::setprecision(14)
                                   << valueb[N * iNode + i] << " ";
                            }

                            os << std::endl;
                        }
                    }
                }
            }

            os << "}" << std::endl;
        }
        messager::barrier();
    }

    return os;
}

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& os,
                         const boundaryConditionDictionary& dict)
{
    if (messager::master())
    {
        os << "index: " << dict.index() << std::endl;
        os << "name : " << dict.name() << std::endl;
        os << "type : " << toString(dict.type()) << std::endl;
    }

    return os;
}

} // namespace accel
