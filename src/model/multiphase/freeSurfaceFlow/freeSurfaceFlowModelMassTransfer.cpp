// File       : freeSurfaceFlowModelMassTransfer.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Interphase mass transfer (cavitation) for the VOF free surface
//              flow model: setup of the models/fields and the rate update.
//              The rate enters the volume fraction equation and the pressure
//              (continuity) equation through the assemblers, see
//              volumeFractionAssemblerNodeTerms.cpp and
//              bulkPressureCorrectionAssemblerNodeTerms.cpp.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "freeSurfaceFlowModel.h"
#include "simulation.h"

namespace accel
{

void freeSurfaceFlowModel::setupMassTransfer_(realm* realm)
{
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    auto isPhase = [&](label index) {
        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (phaseRef(iPhase).index_ == index)
                return true;
        }
        return false;
    };

    label pairCount = 0;
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        pairCount = 0;
        for (const auto& fpm : domain->fluidPairModels_)
        {
            const std::string path = "fluid_pair_models[" +
                                     std::to_string(pairCount++) +
                                     "].mass_transfer";

            const auto& mt = fpm.massTransfer_;
            if (mt.option_ == massTransferModelOption::none)
                continue;

            if (!isPhase(mt.liquidIndex_) || !isPhase(mt.vaporIndex_))
            {
                errorMsg(path + ": liquid and vapor phases must be phases of "
                                "the free surface flow model");
            }

            massTransferPair pair;
            pair.domainIndex_ = domain->index();
            pair.name_ = fpm.materialA_ + "_" + fpm.materialB_;
            pair.model_ = createMassTransferModel(mt);
            pair.liquidIndex_ = mt.liquidIndex_;
            pair.vaporIndex_ = mt.vaporIndex_;

            // nodal rate field: declared once per pair, shared across domains
            const std::string fieldName = "mass_transfer_rate." + pair.name_;
            auto it = mdotMassTransferSTKFieldPtrs_.find(fieldName);
            if (it == mdotMassTransferSTKFieldPtrs_.end())
            {
                auto* fieldPtr = &metaData.declare_field<scalar>(
                    stk::topology::NODE_RANK, fieldName);
                stk::io::set_field_output_type(*fieldPtr, fieldType[1]);
                it = mdotMassTransferSTKFieldPtrs_
                         .emplace(fieldName, fieldPtr)
                         .first;
            }
            pair.mdotSTKFieldPtr_ = it->second;

            const stk::mesh::PartVector& partVec =
                domain->zonePtr()->interiorParts();
            for (const stk::mesh::Part* part : partVec)
            {
                if (!pair.mdotSTKFieldPtr_->defined_on(*part))
                {
                    stk::mesh::put_field_on_mesh(
                        *pair.mdotSTKFieldPtr_, *part, nullptr);
                }
            }

            massTransferPairs_.push_back(std::move(pair));
        }
    }
}

std::vector<const freeSurfaceFlowModel::massTransferPair*>
freeSurfaceFlowModel::massTransferPairs(const domain* domain) const
{
    std::vector<const massTransferPair*> pairs;
    for (const auto& pair : massTransferPairs_)
    {
        if (pair.domainIndex_ == domain->index())
            pairs.push_back(&pair);
    }
    return pairs;
}

void freeSurfaceFlowModel::initializeMassTransferRate(
    const std::shared_ptr<domain> domain)
{
    if (controlsRef().solverRef().restartControl_.isRestart_)
        return;

    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const massTransferPair* pair : massTransferPairs(domain.get()))
    {
        for (stk::mesh::Bucket* bucket : nodeBuckets)
        {
            scalar* mdotb =
                stk::mesh::field_data(*pair->mdotSTKFieldPtr_, *bucket);
            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < bucket->size();
                 ++iNode)
            {
                mdotb[iNode] = 0.0;
            }
        }
    }
}

void freeSurfaceFlowModel::updateMassTransferRate(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    // pressure is stored relative to the domain reference pressure; phase
    // change depends on the absolute pressure
    const scalar pLevel = domain->referencePressure();
    const STKScalarField* pSTKFieldPtr = pRef().stkFieldPtr();

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const massTransferPair* pair : massTransferPairs(domain.get()))
    {
        const STKScalarField* alphaVSTKFieldPtr =
            this->alphaRef(pair->vaporIndex_).stkFieldPtr();
        const STKScalarField* rhoLSTKFieldPtr =
            this->rhoRef(pair->liquidIndex_).stkFieldPtr();
        const STKScalarField* rhoVSTKFieldPtr =
            this->rhoRef(pair->vaporIndex_).stkFieldPtr();

        for (stk::mesh::Bucket* bucket : nodeBuckets)
        {
            const scalar* pb = stk::mesh::field_data(*pSTKFieldPtr, *bucket);
            const scalar* alphaVb =
                stk::mesh::field_data(*alphaVSTKFieldPtr, *bucket);
            const scalar* rhoLb =
                stk::mesh::field_data(*rhoLSTKFieldPtr, *bucket);
            const scalar* rhoVb =
                stk::mesh::field_data(*rhoVSTKFieldPtr, *bucket);
            scalar* mdotb =
                stk::mesh::field_data(*pair->mdotSTKFieldPtr_, *bucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < bucket->size();
                 ++iNode)
            {
                const double raw = pair->model_->rate({pb[iNode] + pLevel,
                                                       alphaVb[iNode],
                                                       rhoLb[iNode],
                                                       rhoVb[iNode]});

                // CFX-like under-relaxation against the previous (relaxed)
                // value that is still stored in the field
                mdotb[iNode] = pair->model_->relax(raw, mdotb[iNode]);
            }
        }
    }
}

} /* namespace accel */
