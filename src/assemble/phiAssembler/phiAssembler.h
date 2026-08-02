// File       : phiAssembler.h
// Created    : Thu Feb 22 2024 13:38:51 (+0100)
// Author     : Fabian Wermelinger
// Description: Assembler for generic systems represented by a single physical
//              variable `phi` (requires phi field to be specified via
//              interface). Coupled assemblers should NOT inherit from this
//              assembler since they require knowledge of multiple physical
//              variables and not just one (`phi` in general). Knowledge about
//              the transport fields is sufficient for this assembler type.
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef PHIASSEMBLER_H
#define PHIASSEMBLER_H

// code
#include "assembler.h"
#include "elementField.h"
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

    void preAssemble(const domain* domain, Context* ctx) override
    {
        computeGamma_(domain);
    }

    void assemble(const domain* domain, Context* ctx) override
    {
        assert(phi_);
        assembleNodeTerms_(domain, ctx); // pointwise kernels
        assembleElemTerms_(domain, ctx); // stencil kernels
    }

    virtual void postAssemble(const domain* domain, Context* ctx) override
    {
        assert(phi_);
        applyConstraints(domain, ctx);
        assembleRelaxation_(domain, ctx->getAMatrix(), phi_->urf());
    }

    void applyConstraints(const domain* domain, Context* ctx)
    {
        // Interface constraints applied under either of two conditions:
        // 1) conformal interfaces (in case non-conformal treatment forces to
        // false)
        if (field_broker_->meshRef().hasInterfaces())
        {
            // select all locally owned nodes for this domain
            const auto& mesh = field_broker_->meshRef();
            const stk::mesh::MetaData& metaData = mesh.metaDataRef();
            const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

            // get fields
            const auto& phiSTKFieldRef = phi_->stkFieldRef();

            auto& A = ctx->getAMatrix();
            auto& b = ctx->getBVector();

            for (const interface* interf : domain->zonePtr()->interfacesRef())
            {
                // constraints on both sides are done in a single pass: during
                // master visit
                if (!interf->isMasterZone(domain->index()))
                    continue;

                if (!interf->isConformalTreatment())
                {
                    // nothing to do
                }
                else
                {
                    // scale slave row by row-1 diagonal (mirrors coupled) so
                    // computeDUCoefficients reads the merged diag, not a bare 1
                    const bool scaledConstraints = true;

                    const auto* graph = A.getGraph();

                    if constexpr (N == 1)
                    {
                        // conformal row-to-row mapping
                        const auto& matchingNodePairConnectivityMap =
                            interf->conformalRowToRowMap(graph);

                        // get pairs
                        const auto& nodePairs =
                            interf->matchingNodePairVector();

                        // matrix connection data
                        const auto& diagOffsets = A.diagOffsetRef();

                        label iPair = 0;
                        for (const auto& nodePair : nodePairs)
                        {
                            // get required local data for the matching pair

                            // data for stencils of node 1 and node 2
                            const auto& node1 = nodePair.first;
                            const auto& node2 = nodePair.second;
                            const label lid1 = static_cast<label>(
                                graph->localToRow(bulkData.local_id(node1)));
                            const label lid2 = static_cast<label>(
                                graph->localToRow(bulkData.local_id(node2)));

                            // pair not in this (subset) graph
                            if (lid1 < 0 || lid2 < 0)
                            {
                                iPair++;
                                continue;
                            }

                            // split pair (a row is off-rank): done by the
                            // cross-rank hook, not the local path
                            if (lid1 >= static_cast<label>(A.nRows()) ||
                                lid2 >= static_cast<label>(A.nRows()))
                            {
                                iPair++;
                                continue;
                            }

                            auto vals1 = A.rowVals(lid1);
                            const label diagOffset1 = diagOffsets[lid1];
                            const scalar phi1 =
                                *stk::mesh::field_data(phiSTKFieldRef, node1);

                            auto vals2 = A.rowVals(lid2);
                            const label diagOffset2 = diagOffsets[lid2];
                            const scalar phi2 =
                                *stk::mesh::field_data(phiSTKFieldRef, node2);

                            // define a mapper from the stencil of node 2 to the
                            // stencil of node 1
                            const std::vector<label>& mapper =
                                matchingNodePairConnectivityMap.map[iPair];

                            // add vals2 to vals1
                            for (label iCol = 0; iCol < mapper.size(); iCol++)
                            {
                                vals1[mapper[iCol]] += vals2[iCol];
                            }

                            // add rhs2 to rhs1
                            b[lid1] += b[lid2];

                            // Force value at node 2 to be equal to that at node
                            // 1

                            // zero row of node 2
                            for (label i = 0; i < vals2.size(); i++)
                            {
                                vals2[i] = 0;
                            }

                            // zero rhs of node 2
                            b[lid2] = 0;

                            // set diagonal of row 2
                            vals2[diagOffset2] = 1;

                            // set off-diagonal of row 2
                            vals2[diagOffset1] = -1;

                            // set rhs for res
                            b[lid2] -= phi2 - phi1;

                            if (scaledConstraints)
                            {
                                // multiply row 2 by diagonal scalar from row 1
                                scalar diag1 = vals1[diagOffset1];

                                // multiply diagonal
                                vals2[diagOffset2] *= diag1;

                                // multiply off-diagonal
                                vals2[diagOffset1] *= diag1;

                                // multiply rhs
                                b[lid2] *= diag1;
                            }

                            // increment
                            iPair++;
                        }
                    }
                    else if constexpr (N == SPATIAL_DIM)
                    {
                        // Get rotation tensor (identity in case of translation
                        // periodicity or general connection)
                        const utils::matrix& rotMat =
                            interf->interfaceSideInfoPtr(domain->index())
                                ->rotationMatrix_;

                        // conformal row-to-row mapping
                        const auto& matchingNodePairConnectivityMap =
                            interf->conformalRowToRowMap(graph);

                        // get pairs
                        const auto& nodePairs =
                            interf->matchingNodePairVector();

                        // matrix connection data
                        const auto& diagOffsets = A.diagOffsetRef();

                        // reusable snapshot buffers for the scaledConstraints
                        // scaling (avoid in-place aliasing when R != I)
                        std::vector<scalar> ws_srcDiag(
                            SPATIAL_DIM * SPATIAL_DIM, 0.0);
                        std::vector<scalar> ws_srcOff(SPATIAL_DIM * SPATIAL_DIM,
                                                      0.0);
                        std::vector<scalar> ws_srcRhs(SPATIAL_DIM, 0.0);

                        // column block, row-major DIM x DIM; fully overwritten
                        // for every column below, so it never needs clearing
                        scalar blk[SPATIAL_DIM * SPATIAL_DIM];

                        label iPair = 0;
                        for (const auto& nodePair : nodePairs)
                        {
                            // get required local data for the matching pair

                            // data for stencils of node 1 and node 2
                            const auto& node1 = nodePair.first;
                            const auto& node2 = nodePair.second;
                            const label lid1 = static_cast<label>(
                                graph->localToRow(bulkData.local_id(node1)));
                            const label lid2 = static_cast<label>(
                                graph->localToRow(bulkData.local_id(node2)));

                            // pair not in this (subset) graph
                            if (lid1 < 0 || lid2 < 0)
                            {
                                iPair++;
                                continue;
                            }

                            // split pair (a row is off-rank): done by the
                            // cross-rank hook, not the local path
                            if (lid1 >= static_cast<label>(A.nRows()) ||
                                lid2 >= static_cast<label>(A.nRows()))
                            {
                                iPair++;
                                continue;
                            }

                            auto vals1 = A.rowVals(lid1);
                            const label diagOffset1 = diagOffsets[lid1];
                            const scalar* phi1 =
                                stk::mesh::field_data(phiSTKFieldRef, node1);

                            auto vals2 = A.rowVals(lid2);
                            const label diagOffset2 = diagOffsets[lid2];
                            const scalar* phi2 =
                                stk::mesh::field_data(phiSTKFieldRef, node2);

                            // define a mapper from the stencil of node 2 to the
                            // stencil of node 1
                            const std::vector<label>& mapper =
                                matchingNodePairConnectivityMap.map[iPair];
                            const std::vector<uint8_t>& isRemapped =
                                matchingNodePairConnectivityMap.remapped[iPair];

                            // add vals2 to vals1 (rotate vals2 first)
                            for (label iCol = 0; iCol < mapper.size(); iCol++)
                            {
                                const bool remappedColumn =
                                    isRemapped[iCol] != 0;

                                const label o2 =
                                    iCol * SPATIAL_DIM * SPATIAL_DIM;
                                const label o1 =
                                    mapper[iCol] * SPATIAL_DIM * SPATIAL_DIM;

                                // blk = A22 * R^T (remapped column) or A22
                                // (interior); split to stay at O(DIM^3)
                                if (remappedColumn)
                                {
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            scalar s = 0.0;
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                                s +=
                                                    vals2[o2 + k * SPATIAL_DIM +
                                                          l] *
                                                    rotMat(j, l);
                                            blk[k * SPATIAL_DIM + j] = s;
                                        }
                                }
                                else
                                {
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                            blk[k * SPATIAL_DIM + j] =
                                                vals2[o2 + k * SPATIAL_DIM + j];
                                }

                                // apply for block (SPATIAL_DIM x SPATIAL_DIM)
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        scalar sum = 0.0;
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                            sum += rotMat(i, k) *
                                                   blk[k * SPATIAL_DIM + j];
                                        vals1[o1 + i * SPATIAL_DIM + j] += sum;
                                    }
                                }
                            }

                            // add rhs2 to rhs1 (rotate rhs2 first)
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    b[lid1 * SPATIAL_DIM + i] +=
                                        rotMat(i, j) *
                                        b[lid2 * SPATIAL_DIM + j];
                                }
                            }

                            // Force value at node 2 to be equal to that at node
                            // 1

                            // zero row of node 2
                            for (label i = 0; i < vals2.size(); i++)
                            {
                                vals2[i] = 0;
                            }

                            // zero rhs of node 2
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                b[lid2 * SPATIAL_DIM + i] = 0.0;
                            }

                            // set diagonal of row 2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                vals2[SPATIAL_DIM * SPATIAL_DIM * diagOffset2 +
                                      i * SPATIAL_DIM + i] = 1;
                            }

                            // set off-diagonal of row 2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    vals2[SPATIAL_DIM * SPATIAL_DIM *
                                              diagOffset1 +
                                          i * SPATIAL_DIM + j] = -rotMat(j, i);
                                }
                            }

                            // set rhs for res
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                b[lid2 * SPATIAL_DIM + i] -= phi2[i];
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    b[lid2 * SPATIAL_DIM + i] +=
                                        (rotMat(j, i) * phi1[j]);
                                }
                            }

                            if (scaledConstraints)
                            {
                                // Rotate diagonal tensor from row 1 (master) to
                                // row 2 (slave) frame Since R rotates from
                                // slave to master: v_master = R * v_slave To
                                // rotate tensor from master to slave: D1_rot =
                                // R^T * D1 * R D1_rot[i,j] = sum_k sum_l R[k,i]
                                // * D1[k,l] * R[l,j]

                                // snapshot source to avoid in-place aliasing
                                // (full D1_rot when R != I)
                                scalar* srcDiag = ws_srcDiag.data();
                                for (label p = 0; p < SPATIAL_DIM * SPATIAL_DIM;
                                     ++p)
                                    srcDiag[p] =
                                        vals2[SPATIAL_DIM * SPATIAL_DIM *
                                                  diagOffset2 +
                                              p];
                                // multiply diagonal block by rotated diagonal
                                // of row 1
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; ++i)
                                    {
                                        scalar sum = 0.0;
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            // rotated diagonal (R^T * D1 *
                                            // R)[i,k]
                                            scalar d_rot_ik = 0.0;
                                            for (label m = 0; m < SPATIAL_DIM;
                                                 ++m)
                                            {
                                                for (label n = 0;
                                                     n < SPATIAL_DIM;
                                                     ++n)
                                                {
                                                    d_rot_ik +=
                                                        rotMat(m, i) *
                                                        vals1[SPATIAL_DIM *
                                                                  SPATIAL_DIM *
                                                                  diagOffset1 +
                                                              m * SPATIAL_DIM +
                                                              n] *
                                                        rotMat(n, k);
                                                }
                                            }
                                            sum += d_rot_ik *
                                                   srcDiag[k * SPATIAL_DIM + j];
                                        }
                                        vals2[SPATIAL_DIM * SPATIAL_DIM *
                                                  diagOffset2 +
                                              i * SPATIAL_DIM + j] = sum;
                                    }
                                }

                                // snapshot source to avoid in-place aliasing
                                // (full D1_rot when R != I)
                                scalar* srcOff = ws_srcOff.data();
                                for (label p = 0; p < SPATIAL_DIM * SPATIAL_DIM;
                                     ++p)
                                    srcOff[p] =
                                        vals2[SPATIAL_DIM * SPATIAL_DIM *
                                                  diagOffset1 +
                                              p];
                                // multiply off-diagonal block by rotated
                                // diagonal of row 1
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; ++i)
                                    {
                                        scalar sum = 0.0;
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            // rotated diagonal (R^T * D1 *
                                            // R)[i,k]
                                            scalar d_rot_ik = 0.0;
                                            for (label m = 0; m < SPATIAL_DIM;
                                                 ++m)
                                            {
                                                for (label n = 0;
                                                     n < SPATIAL_DIM;
                                                     ++n)
                                                {
                                                    d_rot_ik +=
                                                        rotMat(m, i) *
                                                        vals1[SPATIAL_DIM *
                                                                  SPATIAL_DIM *
                                                                  diagOffset1 +
                                                              m * SPATIAL_DIM +
                                                              n] *
                                                        rotMat(n, k);
                                                }
                                            }
                                            sum += d_rot_ik *
                                                   srcOff[k * SPATIAL_DIM + j];
                                        }
                                        vals2[SPATIAL_DIM * SPATIAL_DIM *
                                                  diagOffset1 +
                                              i * SPATIAL_DIM + j] = sum;
                                    }
                                }

                                // snapshot rhs to avoid in-place aliasing
                                scalar* srcRhs = ws_srcRhs.data();
                                for (label p = 0; p < SPATIAL_DIM; ++p)
                                    srcRhs[p] = b[lid2 * SPATIAL_DIM + p];
                                // multiply rhs by rotated diagonal of row 1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    scalar sum = 0.0;
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        // rotated diagonal (R^T * D1 * R)[i,k]
                                        scalar d_rot_ik = 0.0;
                                        for (label m = 0; m < SPATIAL_DIM; ++m)
                                        {
                                            for (label n = 0; n < SPATIAL_DIM;
                                                 ++n)
                                            {
                                                d_rot_ik +=
                                                    rotMat(m, i) *
                                                    vals1[SPATIAL_DIM *
                                                              SPATIAL_DIM *
                                                              diagOffset1 +
                                                          m * SPATIAL_DIM + n] *
                                                    rotMat(n, k);
                                            }
                                        }
                                        sum += d_rot_ik * srcRhs[k];
                                    }
                                    b[lid2 * SPATIAL_DIM + i] = sum;
                                }
                            }

                            // increment
                            iPair++;
                        }
                    }
                    else
                    {
                        errorMsg("1:1 interface treatment is only available "
                                 "for scalar "
                                 "and vectorial transport equations");
                    }

                    // cross-rank fold for split pairs; no-op in base,
                    // overridden in the flow assemblers (developed in .cpp)
                    applyConstraintsCrossRank_(interf, domain, ctx);
                }
            }
        }
    }

    // Cross-rank conformal fold for split pairs: rotation T for vectors,
    // identity for scalars (coupled assemblers call applyCrossRankFold direct).
    virtual void applyConstraintsCrossRank_(const interface* interf,
                                            const domain* domain,
                                            Context* ctx)
    {
        if (!messager::parallel())
            return;

        std::array<scalar, N * N> T{};
        if constexpr (N == 1)
        {
            T[0] = 1.0;
        }
        else
        {
            const utils::matrix& R =
                interf->interfaceSideInfoPtr(domain->index())->rotationMatrix_;
            for (label i = 0; i < static_cast<label>(N); ++i)
                for (label j = 0; j < static_cast<label>(N); ++j)
                    T[i * N + j] = R(i, j);
        }

        const auto& field = phi_->stkFieldRef();
        auto gather = [&](stk::mesh::Entity n, scalar* out)
        {
            const scalar* v = stk::mesh::field_data(field, n);
            for (label i = 0; i < static_cast<label>(N); ++i)
                out[i] = v[i];
        };

        this->applyCrossRankFold(*ctx,
                                 interf->matchingNodePairVector(),
                                 T,
                                 0,
                                 static_cast<int>(N),
                                 gather);
    }

protected:
    FieldType* phi_;
    STKScalarField* GammaSTKFieldPtr_;
    STKScalarField* mDotSTKFieldPtr_;
    STKScalarField* mDotSideSTKFieldPtr_;
    STKScalarField* divUSTKFieldPtr_;
    transportMode transportMode_;

    // Size the Schur-augmented storage on this assembler's
    // linearSolverContext->coefficients to the current globally-
    // consistent CS count.

    void computeGamma_(const domain* domain);
    virtual void
    assembleRelaxation_(const domain* domain, Matrix& A, const scalar urf);
    virtual void assembleBoundaryRelaxation_(const domain* domain,
                                             Context* ctx,
                                             const scalar urf);

    virtual void applySymmetryConditions_(const domain* domain, Context* ctx)
    {
    }

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

    virtual void assembleElemTermsInterfaces_(const domain* domain,
                                              Context* ctx);
    virtual void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);

    // Optional per-IP hook called once per non-exposed interface fragment
    // from inside assembleElemTermsInterfaceSide_. Only invoked for N==1.
    virtual void recordInterfaceFlux_(stk::mesh::Entity /*face*/,
                                      label /*currentGaussPointId*/,
                                      scalar /*fluxValue*/)
    {
    }

    virtual void assembleElemTermsBoundary_(const domain* domain, Context* ctx);

    // boundary conditions
    virtual void assembleElemTermsBoundarySymmetry_(const domain* domain,
                                                    const boundary* boundary,
                                                    Context* ctx)
    {
    }

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
        assembleElemTermsInterfaces_(domain, ctx);
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
    const auto* graph = A.getGraph();
    const scalar urf_inv = 1.0 / urf;
    for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
    {
        const Bucket& bucket = *nodeBuckets[ib];
        const Bucket::size_type n_entities = bucket.size();
        for (Bucket::size_type i = 0; i < n_entities; ++i)
        {
            const stk::mesh::Entity entity = bucket[i];
            const int64_t row = graph->localToRow(bulkData.local_id(entity));
            if (row < 0) // node not part of this (subset) system
                continue;

            // relax diagonal
            assert(row < static_cast<int64_t>(A.nRows()));
            scalar* diag = A.diag(row);
            for (label k = 0; k < BLOCKSIZE; k++)
            {
                diag[BLOCKSIZE * k + k] *= urf_inv;
            }
        }
    }
}

template <size_t N>
void phiAssembler<N>::assembleBoundaryRelaxation_(const domain* domain,
                                                  Context* ctx,
                                                  const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    Vector& b = ctx->getBVector();
    const auto* graph = ctx->getAMatrix().getGraph();

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

    // add interface parts for relaxation under certain circumstances:
    // 1) if the interface is a fluid-solid interface
    // 2) if the interface whether inter-domain are connecting multiple
    //    domains, we only consider nodes nearest to exposed ip's
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            for (auto part :
                 interf->interfaceSideInfoPtr(domain->index())->currentPartVec_)
            {
                partVec.push_back(part);
            }
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
                const int64_t row =
                    graph->localToRow(bulkData.local_id(entity));
                if (row < 0) // node not part of this (subset) system
                    continue;

                // relax diagonal
                assert(row < static_cast<int64_t>(b.size() / BLOCKSIZE));
                scalar* rhs_val = &b[BLOCKSIZE * row];
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
#include "phiAssemblerElemInterfaceConditions.hpp"
#include "phiAssemblerElemTerms.hpp"
#include "phiAssemblerNodeTerms.hpp"

#endif // PHIASSEMBLER_H
