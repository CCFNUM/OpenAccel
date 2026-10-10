// File       : multiphaseModel.cpp
// Created    : Sun Jan 26 2025 22:02:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "multiphaseModel.h"
#include "initialConditions.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

multiphaseModel::multiphaseModel(realm* realm) : flowModel(realm)
{
}

void multiphaseModel::setupMassTransfer_(realm* realm)
{
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    auto isPhase = [&](label index)
    {
        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (phaseRef(iPhase).index_ == index)
                return true;
        }
        return false;
    };

    // nodal field of a pair: declared once, shared across domains
    auto declareField =
        [&](std::map<std::pair<label, label>, STKScalarField*>& fieldPtrs,
            const std::string& fieldName,
            const fluidPairModel& fpm,
            const stk::mesh::PartVector& partVec)
    {
        const auto& mt = fpm.massTransfer_;
        const std::pair<label, label> key = {mt.liquidIndex_, mt.vaporIndex_};
        if (fieldPtrs.find(key) == fieldPtrs.end())
        {
            fieldPtrs[key] = &metaData.declare_field<scalar>(
                stk::topology::NODE_RANK,
                fieldName + "." + mt.liquidPhase_ + "_" + mt.vaporPhase_);
        }
        STKScalarField* fieldPtr = fieldPtrs[key];
        const scalar initialValue = 0.0;
        for (const stk::mesh::Part* part : partVec)
        {
            if (!fieldPtr->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *fieldPtr, *part, 1, &initialValue);
            }
        }
    };

    for (const auto& domain : realm->simulationRef().domainVector())
    {
        for (const auto& fpm : domain->fluidPairModels_)
        {
            const auto& mt = fpm.massTransfer_;
            if (mt.option_ != massTransferModelOption::none)
            {
                const std::string pairName =
                    mt.liquidPhase_ + "_" + mt.vaporPhase_;

                if (!isPhase(mt.liquidIndex_) || !isPhase(mt.vaporIndex_))
                {
                    errorMsg("mass_transfer of pair " + pairName +
                             ": liquid and vapor phases must be phases of "
                             "the multiphase model");
                }

                // continuity source assumes constant phase densities
                if (domain->isMaterialCompressible(
                        domain->globalToLocalMaterialIndex(mt.liquidIndex_)) ||
                    domain->isMaterialCompressible(
                        domain->globalToLocalMaterialIndex(mt.vaporIndex_)))
                {
                    errorMsg("mass_transfer of pair " + pairName +
                             ": compressible phases are not supported");
                }

                if (domain->multiphase_.freeSurfaceModel_
                        .fluxCorrectedTransport_)
                {
                    warningMsg("mass_transfer of pair " + pairName +
                               ": flux corrected transport sharpens the "
                               "volume fraction of the dispersed vapor");
                }

                const std::pair<label, label> key = {mt.liquidIndex_,
                                                     mt.vaporIndex_};
                const bool declared =
                    mdotSTKFieldPtrs_.find(key) != mdotSTKFieldPtrs_.end();

                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();

                declareField(
                    mdotSTKFieldPtrs_, "mass_transfer_rate", fpm, partVec);
                declareField(dmdotdalphaSTKFieldPtrs_,
                             "mass_transfer_rate_alpha_coeff",
                             fpm,
                             partVec);
                declareField(dmdotdpSTKFieldPtrs_,
                             "mass_transfer_rate_pressure_coeff",
                             fpm,
                             partVec);

                // the rate is the state of the relaxation
                if (!declared)
                {
                    STKScalarField* mdotSTKFieldPtr = mdotSTKFieldPtrs_[key];
                    stk::io::set_field_output_type(*mdotSTKFieldPtr,
                                                   fieldType[1]);
                    realm->registerRestartField(mdotSTKFieldPtr->name());
                }
            }
        }
    }
}

void multiphaseModel::initializeMassTransferRate(
    const std::shared_ptr<domain> domain)
{
    // fields are zero by default: restore on restart only
    const auto& restart_ctrl = controlsRef().solverRef().restartControl_;
    if (!restart_ctrl.isRestart_)
        return;

    for (const auto& fpm : domain->fluidPairModels_)
    {
        if (fpm.massTransfer_.option_ != massTransferModelOption::none)
        {
            STKScalarField* mdotSTKFieldPtr =
                mdotSTKFieldPtrs_.at({fpm.massTransfer_.liquidIndex_,
                                      fpm.massTransfer_.vaporIndex_});

            stk::io::MeshField mf(mdotSTKFieldPtr,
                                  std::to_string(std::hash<std::string>{}(
                                      mdotSTKFieldPtr->name())),
                                  restart_ctrl.timeMatchOption_);

            scalar restart_time = restart_ctrl.restartTime_;
            if (restart_time == 0.0)
            {
                restart_time = this->meshRef().ioBrokerRef().get_max_time();
            }
            mf.set_read_time(restart_time);
            mf.set_single_state(false);

            for (const stk::mesh::Part* part :
                 domain->zonePtr()->interiorParts())
            {
                mf.add_subset(*part);
            }

            this->meshRef().ioBrokerRef().read_input_field(mf);

            if (mf.field_restored())
            {
                stk::mesh::communicate_field_data(this->meshRef().bulkDataRef(),
                                                  {mdotSTKFieldPtr});
            }
            else
            {
                warningMsg("Field " + mdotSTKFieldPtr->name() +
                           " not found in the restart file: starting from "
                           "zero");
            }
        }
    }
}

void multiphaseModel::updateMassTransferRate(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    // phase change depends on the absolute pressure
    const scalar pLevel = domain->referencePressure();
    const STKScalarField* pSTKFieldPtr = pRef().stkFieldPtr();

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const auto& fpm : domain->fluidPairModels_)
    {
        const auto& mt = fpm.massTransfer_;
        if (mt.option_ != massTransferModelOption::none)
        {
            const std::pair<label, label> key = {mt.liquidIndex_,
                                                 mt.vaporIndex_};

            // Rayleigh-Plesset: mdot = C (1-a_v) rho_v sqrt(2/3 (p_v-p)/rho_l)
            const scalar pv = mt.saturationPressure_;
            const scalar pMin = mt.pressureClippingForRate_ ? 0.0 : -BIG;
            const scalar Cvap = mt.vaporizationCoefficient_ * 3.0 *
                                mt.nucleationSiteVolumeFraction_ /
                                mt.nucleationSiteRadius_;
            const scalar Ccond =
                mt.condensationCoefficient_ * 3.0 / mt.nucleationSiteRadius_;
            const scalar omega = mt.underRelaxation_;

            // sqrt(|dp|) has an unbounded slope at p_v: linear in dp below
            // dpMin
            const scalar dpMin = 1.0e-3 * pv;

            const STKScalarField* alphaVSTKFieldPtr =
                this->alphaRef(mt.vaporIndex_).stkFieldPtr();
            const STKScalarField* rhoLSTKFieldPtr =
                this->rhoRef(mt.liquidIndex_).stkFieldPtr();
            const STKScalarField* rhoVSTKFieldPtr =
                this->rhoRef(mt.vaporIndex_).stkFieldPtr();
            STKScalarField* mdotSTKFieldPtr = mdotSTKFieldPtrs_.at(key);
            STKScalarField* dmdotdalphaSTKFieldPtr =
                dmdotdalphaSTKFieldPtrs_.at(key);
            STKScalarField* dmdotdpSTKFieldPtr = dmdotdpSTKFieldPtrs_.at(key);

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* pb =
                    stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
                const scalar* alphaVb =
                    stk::mesh::field_data(*alphaVSTKFieldPtr, nodeBucket);
                const scalar* rhoLb =
                    stk::mesh::field_data(*rhoLSTKFieldPtr, nodeBucket);
                const scalar* rhoVb =
                    stk::mesh::field_data(*rhoVSTKFieldPtr, nodeBucket);
                scalar* mdotb =
                    stk::mesh::field_data(*mdotSTKFieldPtr, nodeBucket);
                scalar* dmdotdalphab =
                    stk::mesh::field_data(*dmdotdalphaSTKFieldPtr, nodeBucket);
                scalar* dmdotdpb =
                    stk::mesh::field_data(*dmdotdpSTKFieldPtr, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    const scalar dp = pv - std::max(pb[iNode] + pLevel, pMin);
                    const scalar alphaV =
                        std::min(std::max(alphaVb[iNode], 0.0), 1.0);

                    // growth/collapse mass flux per unit pressure difference
                    const scalar w =
                        rhoVb[iNode] * std::sqrt(2.0 / 3.0 / rhoLb[iNode] /
                                                 std::max(std::abs(dp), dpMin));

                    const bool vaporization = dp > 0.0;
                    const scalar K =
                        (vaporization ? Cvap : Ccond) * w * std::abs(dp);
                    const scalar mdot =
                        vaporization ? K * (1.0 - alphaV) : -K * alphaV;

                    // under-relaxed rate and its (secant in p) sensitivities
                    mdotb[iNode] = omega * mdot + (1.0 - omega) * mdotb[iNode];
                    dmdotdalphab[iNode] = -omega * K;
                    dmdotdpb[iNode] = -omega *
                                      (vaporization ? Cvap * (1.0 - alphaV)
                                                    : Ccond * alphaV) *
                                      w;
                }
            }
        }
    }
}

void multiphaseModel::correctMassTransferRate(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    const STKScalarField* pSTKFieldPtr = pRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtrPrevIter =
        pRef().prevIterRef().stkFieldPtr();

    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (const auto& fpm : domain->fluidPairModels_)
    {
        const auto& mt = fpm.massTransfer_;
        if (mt.option_ != massTransferModelOption::none)
        {
            const std::pair<label, label> key = {mt.liquidIndex_,
                                                 mt.vaporIndex_};

            STKScalarField* mdotSTKFieldPtr = mdotSTKFieldPtrs_.at(key);
            const STKScalarField* dmdotdpSTKFieldPtr =
                dmdotdpSTKFieldPtrs_.at(key);

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const label nNodesPerBucket = nodeBucket.size();

                // field chunks in bucket
                const scalar* pb =
                    stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
                const scalar* pbPrevIter =
                    stk::mesh::field_data(*pSTKFieldPtrPrevIter, nodeBucket);
                const scalar* dmdotdpb =
                    stk::mesh::field_data(*dmdotdpSTKFieldPtr, nodeBucket);
                scalar* mdotb =
                    stk::mesh::field_data(*mdotSTKFieldPtr, nodeBucket);

                for (label iNode = 0; iNode < nNodesPerBucket; ++iNode)
                {
                    mdotb[iNode] +=
                        dmdotdpb[iNode] * (pb[iNode] - pbPrevIter[iNode]);
                }
            }
        }
    }
}

} /* namespace accel */
