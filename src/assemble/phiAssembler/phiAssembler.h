// File : phiAssembler.h
// Created : Thu Feb 22 2024 13:38:51 (+0100)
// Author : Fabian Wermelinger
// Description: Assembler for generic systems represented by a single physical
// variable `phi` (requires phi field to be specified via
// interface). Coupled assemblers should NOT inherit from this
// assembler since they require knowledge of multiple physical
// variables and not just one (`phi` in general). Knowledge about
// the transport fields is sufficient for this assembler type.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PHIASSEMBLER_H
#define PHIASSEMBLER_H

// code
#include "assembler.h"
#include "elementField.h"
#include "macros.h"
#include "mesh.h"
#include "nodeField.h"
#include "types.h"

namespace accel
{

enum transportMode
{
    null = 0,
    advection,
    diffusion,
    advectionDiffusion
};

template <size_t N>
class phiAssembler : public assembler<N>
{
public:
    using Base = assembler<N>;
    using Base::BLOCKSIZE;
    using Context = typename Base::Context;
    using FieldType = nodeField<N, N * SPATIAL_DIM>;

    using Bucket = stk::mesh::Bucket;
    using BucketVector = stk::mesh::BucketVector;
    using GammaFunction =
        std::function<void(const domain* domain, STKScalarField& gamma)>;

protected:
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;
    using Base::field_broker_;

    enum class GammaType
    {
        Constant,
        External,
        Function,
        NotSet
    };

public:
    phiAssembler(fieldBroker* field_broker)
        : Base(field_broker), phi_(nullptr), GammaSTKFieldPtr_(nullptr),
          transportMode_(null), GammaFunc_(nullptr),
          GammaType_(GammaType::NotSet)
    {
    }

    void setup(FieldType* phi,
               transportMode transportMode,
               const std::vector<std::shared_ptr<domain>>& domains,
               const scalar gamma = 0.0)
    {
        phi_ = phi;
        transportMode_ = transportMode;

        switch (transportMode_)
        {
            case advection:
                {
                    GammaType_ = GammaType::Constant;

                    // force gamma to 0.0
                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), 0.0);
                    }

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case diffusion:
                {
                    GammaType_ = GammaType::Constant;

                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), gamma);
                    }

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case advectionDiffusion:
                {
                    GammaType_ = GammaType::Constant;

                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), gamma);
                    }

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            default:
                break;
        }
    }

    void setup(FieldType* phi,
               transportMode transportMode,
               const std::vector<std::shared_ptr<domain>>& domains,
               STKScalarField* gamma)
    {
        phi_ = phi;
        transportMode_ = transportMode;

        switch (transportMode_)
        {
            case advection:
                {
                    GammaType_ = GammaType::Constant;

                    // force gamma to 0.0
                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), 0.0);
                    }

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case diffusion:
                {
                    GammaType_ = GammaType::External;
                    STK_ThrowRequireMsg(
                        gamma->entity_rank() == stk::topology::NODE_RANK,
                        "phiAssembler::setup: field `"
                            << gamma->name() << "` has entity-rank "
                            << gamma->entity_rank()
                            << " but should have NODE_RANK");
                    GammaSTKFieldPtr_ = gamma;

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case advectionDiffusion:
                {
                    GammaType_ = GammaType::External;
                    STK_ThrowRequireMsg(
                        gamma->entity_rank() == stk::topology::NODE_RANK,
                        "phiAssembler::setup: field `"
                            << gamma->name() << "` has entity-rank "
                            << gamma->entity_rank()
                            << " but should have NODE_RANK");
                    GammaSTKFieldPtr_ = gamma;

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            default:
                break;
        }
    }

    void setup(FieldType* phi,
               transportMode transportMode,
               const std::vector<std::shared_ptr<domain>>& domains,
               GammaFunction f)
    {
        phi_ = phi;
        transportMode_ = transportMode;

        switch (transportMode_)
        {
            case advection:
                {
                    GammaType_ = GammaType::Constant;

                    // force gamma to 0.0
                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), 0.0);
                    }

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case diffusion:
                {
                    GammaType_ = GammaType::Function;
                    GammaFunc_ = f;

                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), 0.0);
                    }

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            case advectionDiffusion:
                {
                    GammaType_ = GammaType::Function;
                    GammaFunc_ = f;

                    for (const auto& domain : domains)
                    {
                        initializeGamma_(domain.get(), 0.0);
                    }

                    // assign mDot pointers
                    mDotSTKFieldPtr_ = this->mDotRef().stkFieldPtr();
                    mDotSideSTKFieldPtr_ =
                        this->mDotRef().sideFieldPtr()
                            ? this->mDotRef().sideFieldRef().stkFieldPtr()
                            : nullptr;
                    divUSTKFieldPtr_ = this->divRef().stkFieldPtr();
                    assert(mDotSTKFieldPtr_);
                    assert(divUSTKFieldPtr_);

                    // force creation of 0-dummy beta field in case a
                    // high-resolution is not required
                    if (phi_->blendingFactorPtr() == nullptr)
                    {
                        phi_->setupBlendingFactorField(
                            /*enable only a dummy field*/ true);
                    }

                    assert(phi_);
                    assert(GammaSTKFieldPtr_);
                }
                break;

            default:
                break;
        }
    }

protected:
    FieldType* phi_;
    STKScalarField* GammaSTKFieldPtr_;
    STKScalarField* mDotSTKFieldPtr_;
    STKScalarField* mDotSideSTKFieldPtr_;
    STKScalarField* divUSTKFieldPtr_;
    transportMode transportMode_;

    void preAssemble_(const domain* domain, Context*) override
    {
        computeGamma_(domain);
    }

    virtual void postAssemble_(const domain* domain, Context* ctx) override
    {
        assert(phi_);
        assembleRelaxation_(domain, ctx->getAMatrix(), phi_->urf());
    }

    void assemble_(const domain* domain, Context* ctx) override
    {
        assert(phi_);
        assembleNodeTerms_(domain, ctx); // pointwise kernels
        assembleElemTerms_(domain, ctx); // stencil kernels
    }

    void computeGamma_(const domain* domain);
    virtual void
    assembleRelaxation_(const domain* domain, Matrix& A, const scalar urf);
    virtual void assembleBoundaryRelaxation_(const domain* domain,
                                             Vector& b,
                                             const scalar urf);

    // kernel drivers
    virtual void assembleNodeTermsFused_(const domain* domain, Context* ctx);
    virtual void assembleNodeTermsFusedSteady_(const domain* domain,
                                               Context* ctx);
    virtual void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                           Context* ctx);
    virtual void
    assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                               Context* ctx);
    virtual void assembleElemTermsInterior_(const domain* domain, Context* ctx);

    virtual void assembleElemTermsBoundary_(const domain* domain, Context* ctx);

    // boundary conditions
    virtual void
    assembleElemTermsBoundaryWallFixedValue_(const domain* domain,
                                             const boundary* boundary,
                                             Context* ctx);
    virtual void
    assembleElemTermsBoundaryWallZeroGradient_(const domain* domain,
                                               const boundary* boundary,
                                               Context* ctx);
    virtual void
    assembleElemTermsBoundaryWallSpecifiedFlux_(const domain* domain,
                                                const boundary* boundary,
                                                Context* ctx);
    virtual void assembleElemTermsBoundaryWallMixed_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx);
    virtual void
    assembleElemTermsBoundaryInletFixedValue_(const domain* domain,
                                              const boundary* boundary,
                                              Context* ctx);
    virtual void
    assembleElemTermsBoundaryOutletZeroGradient_(const domain* domain,
                                                 const boundary* boundary,
                                                 Context* ctx);
    virtual void assembleElemTermsBoundaryOpening_(const domain* domain,
                                                   const boundary* boundary,
                                                   Context* ctx);

    // Auxiliary field access

    virtual nodeField<1, SPATIAL_DIM>& rhoRef()
    {
        return field_broker_->rhoRef();
    }

    virtual const nodeField<1, SPATIAL_DIM>& rhoRef() const
    {
        return field_broker_->rhoRef();
    }

    virtual elementField<scalar, 1>& mDotRef()
    {
        return field_broker_->mDotRef();
    }

    virtual const elementField<scalar, 1>& mDotRef() const
    {
        return field_broker_->mDotRef();
    }

    virtual nodeField<1>& divRef()
    {
        return field_broker_->mDotRef().divRef();
    }

    const virtual nodeField<1>& divRef() const
    {
        return field_broker_->mDotRef().divRef();
    }

    virtual std::string getCoordinatesID_(const domain* /*domain*/) const
    {
        return mesh::coordinates_ID;
    }

    virtual std::string getDualNodalVolumeID_(const domain* /*domain*/) const
    {
        return mesh::dual_nodal_volume_ID;
    }

    virtual std::string getExposedAreaVectorID_(const domain* /*domain*/) const
    {
        return mesh::exposed_area_vector_ID;
    }

private:
    GammaFunction GammaFunc_;
    GammaType GammaType_;

    void initializeGamma_(const domain* domain, const scalar v);

    void assembleNodeTerms_(const domain* domain, Context* ctx)
    {
        assembleNodeTermsFused_(domain, ctx);
    }

    void assembleElemTerms_(const domain* domain, Context* ctx)
    {
        assembleElemTermsInterior_(domain, ctx);
        assembleElemTermsBoundary_(domain, ctx);
    }
};

template <size_t N>
void phiAssembler<N>::initializeGamma_(const domain* domain, const scalar v)
{
    assert(phi_);

    auto& mesh = field_broker_->meshRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // For next iterations, STK will only retrieve Gamma, not create it
    GammaSTKFieldPtr_ = &metaData.declare_field<scalar>(
        stk::topology::NODE_RANK, "Gamma_" + phi_->name());

    // Put the stk field on interior mesh parts for every active zone

    // put domain parts on mesh
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    for (const stk::mesh::Part* part : partVec)
    {
        // check if already defined from a previous pass
        if (!GammaSTKFieldPtr_->defined_on(*part))
        {
            stk::mesh::put_field_on_mesh(
                *GammaSTKFieldPtr_, *part, 1, nullptr); // scalar field
        }
    }

    // initialize field
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // define some common selectors; select all nodes
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    const BucketVec& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (Bucket::size_type ib = 0; ib < nodeBuckets.size(); ib++)
    {
        Bucket& nodeBucket = *nodeBuckets[ib];
        const Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // field chunks in bucket
        scalar* Gammab = stk::mesh::field_data(*GammaSTKFieldPtr_, nodeBucket);

        for (Bucket::size_type iNode = 0; iNode < nNodesPerBucket; iNode++)
        {
            Gammab[iNode] = v;
        }
    }
}

template <size_t N>
void phiAssembler<N>::computeGamma_(const domain* domain)
{
    using BucketVec = stk::mesh::BucketVector;

    if (GammaType_ == GammaType::Function)
    {
        assert(GammaSTKFieldPtr_);
        GammaFunc_(domain, *GammaSTKFieldPtr_);
    }
    else if (GammaType_ == GammaType::NotSet)
    {
        errorMsg(
            "phiAssembler::computeGamma_: type of Gamma evaluation is not set");
    }
}

template <size_t N>
void phiAssembler<N>::assembleRelaxation_(const domain* domain,
                                          Matrix& A,
                                          const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    // select all locally owned nodes for this domain
    const auto& mesh = field_broker_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    const BucketVec& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

    // loop over local nodes and relax associated matrix rows
    const scalar urf_inv = 1.0 / urf;
    for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
    {
        const Bucket& bucket = *nodeBuckets[ib];
        const Bucket::size_type n_entities = bucket.size();
        for (Bucket::size_type i = 0; i < n_entities; ++i)
        {
            const stk::mesh::Entity entity = bucket[i];
            const auto lid = bulkData.local_id(entity);

            // relax diagonal
            assert(lid < A.nRows());
            scalar* diag = A.diag(lid);
            for (label k = 0; k < BLOCKSIZE; k++)
            {
                diag[BLOCKSIZE * k + k] *= urf_inv;
            }
        }
    }
}

template <size_t N>
void phiAssembler<N>::assembleBoundaryRelaxation_(const domain* domain,
                                                  Vector& b,
                                                  const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    // select all locally owned nodes for this domain
    const auto& mesh = field_broker_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const zone* zonePtr = domain->zonePtr();

    // Collect all boundary parts: It is ensured in this way that no
    // node duplications are made for adjacent boundary parts
    stk::mesh::PartVector partVec;

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
        const stk::mesh::PartVector& parts = boundaryRef.parts();

        boundaryPhysicalType type = boundaryRef.type();
        switch (type)
        {
            case boundaryPhysicalType::wall:
                {
                    if (domain->type() == domainType::fluid)
                    {
                        const boundaryConditionType UBCType =
                            field_broker_->URef()
                                .boundaryConditionRef(domain->index(),
                                                      iBoundary)
                                .type();

                        switch (UBCType)
                        {
                            case boundaryConditionType::noSlip:
                                {
                                    for (auto part : parts)
                                    {
                                        partVec.push_back(part);
                                    }
                                }
                                break;

                            case boundaryConditionType::slip:
                                {
                                    if (mesh.controlsRef().isTransient())
                                    {
                                        // only no-slip walls in case of
                                        // transient are relaxed
                                    }
                                    else
                                    {
                                        for (auto part : parts)
                                        {
                                            partVec.push_back(part);
                                        }
                                    }
                                }
                                break;

                            default:
                                errorMsg("invalid velocity boundary "
                                         "condition at wall");
                        }
                    }
                    else
                    {
                        for (auto part : parts)
                        {
                            partVec.push_back(part);
                        }
                    }
                }
                break;

            default:
                {
                    if (mesh.controlsRef().isTransient())
                    {
                        // only no-slip walls in case of transient are relaxed
                    }
                    else
                    {
                        for (auto part : parts)
                        {
                            partVec.push_back(part);
                        }
                    }
                }
                break;
        }
    }

    // Apply relaxation
    {
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

        const BucketVec& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        // loop over local nodes and relax associated matrix rows
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& bucket = *nodeBuckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            for (Bucket::size_type i = 0; i < n_entities; ++i)
            {
                const stk::mesh::Entity entity = bucket[i];
                const auto lid = bulkData.local_id(entity);

                // relax diagonal
                assert(lid < b.size() / BLOCKSIZE);
                scalar* rhs_val = &b[BLOCKSIZE * lid];
                for (label k = 0; k < BLOCKSIZE; k++)
                {
                    rhs_val[k] *= urf;
                }
            }
        }
    }
}

} /* namespace accel */

// kernel implementations
#include "phiAssemblerElemBoundaryConditions.hpp"
#include "phiAssemblerElemTerms.hpp"
#include "phiAssemblerNodeTerms.hpp"

#endif // PHIASSEMBLER_H
