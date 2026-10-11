// File       : fsiConvergence.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: FSI convergence checks
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "fsiConvergence.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#include "messager.h"
#include "pressureCorrectionEquation.h"
#include "segregatedFlowEquations.h"
#include "segregatedFreeSurfaceFlowEquations.h"
#include "simulation.h"
#include "solidDisplacementEquation.h"
#include "version.h"

namespace accel
{

fsiConvergence::fsiConvergence(simulation& sim) : sim_(sim)
{
    for (const auto& domainPtr : sim_.domainVector_)
    {
        for (const interface* interf : domainPtr->interfacesRef())
        {
            if (interf->isFluidSolidType() && interf->isRobinNeumann())
            {
                isRobinNeumann_ = true;
            }
            if (interf->isFluidSolidType() && interf->writeForceBalance())
            {
                writeForceBalance_ = true;
            }
        }
    }

    if (enabled() && !hasCriteria() && messager::master())
    {
        warningMsg("fsi_convergence is enabled without criteria: the outer "
                   "loop convergence uses the equation residuals; FSI "
                   "diagnostics/force balance are still written");
    }

    // RN parameter validity depends on the solid material only, hence is
    // fixed for the run (collective call: constructed on all ranks)
    if (isRobinNeumann_)
    {
        rnFallback_ =
            !computeRobinNeumannParams(
                 sim_, sim_.domainVector_, sim_.controlsRef().getTimestep())
                 .valid;
    }
}

bool fsiConvergence::hasCriteria() const
{
    return !sim_.controlsRef()
                .solverRef()
                .solverControl_.basicSettings_.convergenceCriteria_
                .fsiConvergence_.criteria_.empty();
}

bool fsiConvergence::enabled() const
{
    const auto& fsiConv =
        sim_.controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceCriteria_.fsiConvergence_;
    return fsiConv.enabled_;
}

void fsiConvergence::resetForTimeStep()
{
    if (!enabled())
    {
        return;
    }

    fsiInterfaceResidualNormMax_.clear();
    fsiInterfaceMaxTotalDisplNorm_.clear();
    fsiInterfaceResidualNorms_.clear();
    fsiInterfaceResidualNorm_ = 0.0;

    fsiInterfaceTractionPrev_.clear();
    fsiInterfaceTractionResidualNormMax_.clear();
    fsiForceResidualNorms_.clear();
    fsiForceResidualNorm_ = 0.0;

    fsiInterfacePCorrPrev_.clear();
    interfacePCorrChangeNorm_ = 0.0;
    interfacePCorrChangeNormRef_ = 1.0;
    robinAccelTermRms_ = 0.0;
    robinPCorrFluxRms_ = 0.0;
}

void fsiConvergence::update()
{
    if (!enabled())
    {
        return;
    }

    const auto& fsiConv =
        sim_.controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceCriteria_.fsiConvergence_;

    fsiInterfaceResidualNorm_ = 0.0;
    fsiForceResidualNorm_ = 0.0;

    for (const auto& criterion : fsiConv.criteria_)
    {
        switch (criterion)
        {
            case fsiConvergenceType::fsiInterfaceResidual:
                updateFsiInterfaceResidual_(fsiConv.printResiduals_);
                break;
            case fsiConvergenceType::fsiForceResidual:
                updateFsiForceResidual_(fsiConv.printResiduals_);
                break;
            default:
                break;
        }
    }

    if (isRobinNeumann_)
    {
        updateInterfacePCorrChange_();
    }

    if (writeForceBalance_)
    {
        updateForceBalance_();
    }

    // Write combined diagnostics row if enabled
    const auto& fsiConvCfg =
        sim_.controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceCriteria_.fsiConvergence_;
    if (fsiConvCfg.diagnosticsEnabled_)
    {
        if (!diagnosticsStream_)
            initializeDiagnosticsFile_();

        writeDiagnosticsRow_(sim_.getIterationCount(),
                             fsiInterfaceResidualNorm_,
                             fsiForceResidualNorm_,
                             interfacePCorrChangeNorm_,
                             robinAccelTermRms_,
                             robinPCorrFluxRms_);
    }
}

void fsiConvergence::updateFsiInterfaceResidual_(bool writeResiduals)
{
    solidDisplacementEquation* solidEq = nullptr;
    for (auto& equation : sim_.equationVector_)
    {
        if (equation->getID() == equationID::solidDisplacement)
        {
            solidEq = dynamic_cast<solidDisplacementEquation*>(equation.get());
            break;
        }
    }

    if (!solidEq)
    {
        return;
    }

    auto& DField = solidEq->DRef().stkFieldRef();
    // raw (un-relaxed) solid correction X^k, transferred to the fluid side
    auto& XField = solidEq->solidCorrectionRef().stkFieldRef();
    std::unordered_set<label> visited;

    for (const auto& domain : sim_.domainVector_)
    {
        for (const interface* interf : domain->interfacesRef())
        {
            if (!interf->isFluidSolidType())
            {
                continue;
            }
            if (!visited.insert(interf->index()).second)
            {
                continue;
            }

            const label masterIdx = interf->masterZoneIndex();
            const label slaveIdx = interf->slaveZoneIndex();
            const label fluidZoneIndex =
                (sim_.domainRef(masterIdx).type() == domainType::fluid)
                    ? masterIdx
                    : slaveIdx;

            const interfaceSideInfo* fluidSide =
                interf->interfaceSideInfoPtr(fluidZoneIndex);

            // locally owned nodes only: shared nodes must enter the global
            // sums once
            stk::mesh::Selector selFluidNodes =
                DField.mesh_meta_data().locally_owned_part() &
                stk::mesh::selectUnion(fluidSide->currentPartVec_);
            const stk::mesh::BucketVector& fluidNodeBuckets =
                DField.get_mesh().get_buckets(stk::topology::NODE_RANK,
                                              selFluidNodes);

            size_t nTotal = 0;
            for (auto ib = fluidNodeBuckets.begin();
                 ib != fluidNodeBuckets.end();
                 ++ib)
            {
                nTotal += (*ib)->size();
            }
            const size_t vecSize = nTotal * SPATIAL_DIM;

            std::vector<scalar> DCurrent(vecSize, 0.0);
            std::vector<scalar> XCurrent(vecSize, 0.0);
            {
                size_t offset = 0;
                for (auto ib = fluidNodeBuckets.begin();
                     ib != fluidNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& b = **ib;
                    const scalar* Db = stk::mesh::field_data(DField, b);
                    const scalar* Xb = stk::mesh::field_data(XField, b);
                    for (size_t iNode = 0; iNode < b.size(); ++iNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            DCurrent[offset] = Db[SPATIAL_DIM * iNode + i];
                            XCurrent[offset] = Xb[SPATIAL_DIM * iNode + i];
                            offset++;
                        }
                    }
                }
            }

            const label interfIdx = interf->index();

            // r = X^k (un-relaxed correction); iteration 1 counts but cannot
            // converge
            scalar normSq = 0.0;
            scalar totalDisplNormSq = 0.0;
            for (size_t i = 0; i < vecSize; ++i)
            {
                const scalar r = XCurrent[i];
                normSq += r * r;
                totalDisplNormSq += DCurrent[i] * DCurrent[i];
            }
            messager::sumReduce(normSq);
            messager::sumReduce(totalDisplNormSq);
            const scalar norm = std::sqrt(normSq);
            const scalar totalDisplNorm = std::sqrt(totalDisplNormSq);

            // residual = min(|X|/max|X|, |X|/max|D|) (s4f dual norm)
            auto& maxNorm = fsiInterfaceResidualNormMax_[interfIdx];
            maxNorm = std::max(maxNorm, norm);
            const scalar normRel1 = norm / (maxNorm + SMALL);

            auto& maxTotalDispl = fsiInterfaceMaxTotalDisplNorm_[interfIdx];
            maxTotalDispl = std::max(maxTotalDispl, totalDisplNorm);
            const scalar normRel2 = norm / (maxTotalDispl + SMALL);

            const scalar normRel = std::min(normRel1, normRel2);
            fsiInterfaceResidualNorms_[interfIdx] = normRel;
            fsiInterfaceResidualNorm_ =
                std::max(fsiInterfaceResidualNorm_, normRel);

            if (messager::master())
            {
                std::cout << "  FSI Residual [" << interf->name() << "]"
                          << "  |r|_norm=" << std::scientific
                          << std::setprecision(4) << normRel << std::endl;
            }

            if (writeResiduals)
            {
                auto& streams = residualStreams_["fsi_interface_residual"];
                if (streams.find(interfIdx) == streams.end())
                {
                    initializeResidualFile_(
                        interfIdx, interf->name(), "fsi_interface_residual");
                }
                writeResidualLine_(
                    "fsi_interface_residual", interfIdx, normRel);
            }
        }
    }
}

void fsiConvergence::updateFsiForceResidual_(bool writeResiduals)
{
    solidDisplacementEquation* solidEq = nullptr;
    for (auto& equation : sim_.equationVector_)
    {
        if (equation->getID() == equationID::solidDisplacement)
        {
            solidEq = dynamic_cast<solidDisplacementEquation*>(equation.get());
            break;
        }
    }

    if (!solidEq)
    {
        return;
    }

    // traction: D's side-flux field, written by the solid model
    auto& tractionSTKField = solidEq->DRef().sideFluxFieldRef().stkFieldRef();
    std::unordered_set<label> visited;

    for (const auto& domain : sim_.domainVector_)
    {
        for (const interface* interf : domain->interfacesRef())
        {
            if (!interf->isFluidSolidType())
            {
                continue;
            }
            if (!visited.insert(interf->index()).second)
            {
                continue;
            }

            const label masterIdx = interf->masterZoneIndex();
            const label slaveIdx = interf->slaveZoneIndex();
            const label fluidZoneIndex =
                (sim_.domainRef(masterIdx).type() == domainType::fluid)
                    ? masterIdx
                    : slaveIdx;

            const interfaceSideInfo* fluidSide =
                interf->interfaceSideInfoPtr(fluidZoneIndex);

            stk::mesh::MetaData& metaData = tractionSTKField.mesh_meta_data();
            stk::mesh::BulkData& bulkData = tractionSTKField.get_mesh();

            // side-rank field: locally owned sides of the fluid-side
            // interface parts (each side counted once in the global sums)
            stk::mesh::Selector selFluidSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(fluidSide->currentPartVec_);
            const stk::mesh::BucketVector& fluidSideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selFluidSides);

            // gather the traction in bucket/side/ip order
            std::vector<scalar> tCurrent;
            for (const stk::mesh::Bucket* bucket : fluidSideBuckets)
            {
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        bucket->topology());
                const label numScsBip = meFC->numIntPoints_;

                for (const stk::mesh::Entity side : *bucket)
                {
                    const scalar* tractionValues =
                        stk::mesh::field_data(tractionSTKField, side);
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            tCurrent.push_back(
                                tractionValues[SPATIAL_DIM * ip + j]);
                        }
                    }
                }
            }

            const label interfIdx = interf->index();
            auto& prev = fsiInterfaceTractionPrev_[interfIdx];

            // iteration 1: global decision (no previous traction); reports 1.0
            const bool firstIteration = (sim_.getIterationCount() <= 1);
            const bool hasPrev =
                !firstIteration && (prev.size() == tCurrent.size());

            scalar normSq = 0.0;
            scalar tractionNormSq = 0.0;
            for (size_t i = 0; i < tCurrent.size(); ++i)
            {
                if (hasPrev)
                {
                    const scalar r = tCurrent[i] - prev[i];
                    normSq += r * r;
                }
                tractionNormSq += tCurrent[i] * tCurrent[i];
            }
            messager::sumReduce(normSq);
            messager::sumReduce(tractionNormSq);
            const scalar norm = std::sqrt(normSq);
            const scalar tractionNorm = std::sqrt(tractionNormSq);

            auto& maxNorm = fsiInterfaceTractionResidualNormMax_[interfIdx];
            maxNorm = std::max(maxNorm, tractionNorm);

            const scalar normRel =
                firstIteration ? 1.0 : norm / (maxNorm + SMALL);
            fsiForceResidualNorms_[interfIdx] = normRel;
            fsiForceResidualNorm_ = std::max(fsiForceResidualNorm_, normRel);

            prev = tCurrent;

            if (messager::master())
            {
                std::cout << "  FSI Force Residual [" << interf->name() << "]"
                          << "  |r|_norm=" << std::scientific
                          << std::setprecision(4) << normRel << std::endl;
            }

            if (writeResiduals)
            {
                auto& streams = residualStreams_["fsi_force_residual"];
                if (streams.find(interfIdx) == streams.end())
                {
                    initializeResidualFile_(
                        interfIdx, interf->name(), "fsi_force_residual");
                }
                writeResidualLine_("fsi_force_residual", interfIdx, normRel);
            }
        }
    }
}

void fsiConvergence::updateInterfacePCorrChange_()
{
    // diagnostics: p' change, Robin a_n term and Q_vol RMS at the RN interface
    // (the pressure equation is reached through the segregated wrapper)
    const stk::mesh::Field<scalar>* pCorrSTKFieldPtr = nullptr;
    const std::unordered_map<uint64_t, scalar>* robinAccelMapPtr = nullptr;
    const stk::mesh::Field<scalar>* rnWallFluxSideSTKFieldPtr = nullptr;
    for (auto& eq : sim_.equationVector_)
    {
        if (auto* sfe = dynamic_cast<segregatedFlowEquations*>(eq.get()))
        {
            pressureCorrectionEquation* pcEq = sfe->pCorrEqPtr();
            if (pcEq)
            {
                pCorrSTKFieldPtr =
                    &pcEq->flowModelPtr()->pCorrRef().stkFieldRef();
                robinAccelMapPtr = &pcEq->assemblerPtr()->robinAccelTermMap();
                rnWallFluxSideSTKFieldPtr = &pcEq->flowModelPtr()
                                                 ->rnWallFluxRef()
                                                 .sideFieldRef()
                                                 .stkFieldRef();
            }
            break;
        }
        if (auto* sfse =
                dynamic_cast<segregatedFreeSurfaceFlowEquations*>(eq.get()))
        {
            bulkPressureCorrectionEquation* bpcEq = sfse->pCorrEqPtr();
            if (bpcEq)
            {
                pCorrSTKFieldPtr = &sfse->pCorrRef().stkFieldRef();
                robinAccelMapPtr = &bpcEq->robinAccelTermMap();
                rnWallFluxSideSTKFieldPtr =
                    &sfse->rnWallFluxRef().sideFieldRef().stkFieldRef();
            }
            break;
        }
    }
    if (!pCorrSTKFieldPtr)
    {
        interfacePCorrChangeNorm_ = 0.0;
        robinAccelTermRms_ = 0.0;
        robinPCorrFluxRms_ = 0.0;
        return;
    }

    // RMS of the Robin acceleration term over all RN BIPs
    {
        scalar accelSumSq = 0.0;
        label accelCount = 0;
        for (const auto& kv : *robinAccelMapPtr)
        {
            accelSumSq += kv.second * kv.second;
            accelCount++;
        }
        messager::sumReduce(accelSumSq);
        messager::sumReduce(accelCount);
        robinAccelTermRms_ =
            (accelCount > 0)
                ? std::sqrt(accelSumSq / static_cast<scalar>(accelCount))
                : 0.0;
    }

    const auto& pCorrSTKFieldRef = *pCorrSTKFieldPtr;

    scalar sumSq = 0.0;
    label count = 0;
    std::unordered_set<label> visited;

    // RMS of Q_vol/|A| over the fluid-side RN BIPs (current area)
    scalar pCorrFluxSumSq = 0.0;
    label pCorrFluxCount = 0;
    const stk::mesh::Field<scalar>* exposedAreaVecSTKFieldPtr =
        pCorrSTKFieldRef.mesh_meta_data().get_field<scalar>(
            pCorrSTKFieldRef.mesh_meta_data().side_rank(),
            mesh::exposed_area_vector_ID);

    for (const auto& domainPtr : sim_.domainVector_)
    {
        for (const interface* interf : domainPtr->interfacesRef())
        {
            if (!interf->isFluidSolidType() || !interf->isRobinNeumann())
            {
                continue;
            }
            if (!visited.insert(interf->index()).second)
            {
                continue;
            }

            const label masterIdx = interf->masterZoneIndex();
            const label slaveIdx = interf->slaveZoneIndex();
            const label fluidZoneIndex =
                (sim_.domainRef(masterIdx).type() == domainType::fluid)
                    ? masterIdx
                    : slaveIdx;

            const interfaceSideInfo* fluidSide =
                interf->interfaceSideInfoPtr(fluidZoneIndex);

            stk::mesh::MetaData& metaData = pCorrSTKFieldRef.mesh_meta_data();
            stk::mesh::BulkData& bulkData = pCorrSTKFieldRef.get_mesh();

            // robin_pCorr_flux_rms: locally owned fluid-side sides carrying
            // the rnWallFlux side field (per-BIP values)
            if (rnWallFluxSideSTKFieldPtr && exposedAreaVecSTKFieldPtr)
            {
                stk::mesh::Selector selOwnedFluxSides =
                    metaData.locally_owned_part() &
                    stk::mesh::selectUnion(fluidSide->currentPartVec_) &
                    stk::mesh::selectField(*rnWallFluxSideSTKFieldPtr);
                const stk::mesh::BucketVector& fluxSideBuckets =
                    bulkData.get_buckets(metaData.side_rank(),
                                         selOwnedFluxSides);

                for (const stk::mesh::Bucket* bucket : fluxSideBuckets)
                {
                    MasterElement* meFC =
                        MasterElementRepo::get_surface_master_element(
                            bucket->topology());
                    const label numScsBip = meFC->numIntPoints_;

                    for (const stk::mesh::Entity side : *bucket)
                    {
                        const scalar* qVals = stk::mesh::field_data(
                            *rnWallFluxSideSTKFieldPtr, side);
                        const scalar* areaVec = stk::mesh::field_data(
                            *exposedAreaVecSTKFieldPtr, side);

                        for (label ip = 0; ip < numScsBip; ++ip)
                        {
                            scalar aMagSq = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    areaVec[ip * SPATIAL_DIM + j];
                                aMagSq += axj * axj;
                            }
                            const scalar aMag = std::sqrt(aMagSq);
                            if (!(aMag > SMALL))
                            {
                                continue;
                            }
                            const scalar uQ = qVals[ip] / aMag;
                            pCorrFluxSumSq += uQ * uQ;
                            pCorrFluxCount++;
                        }
                    }
                }
            }

            stk::mesh::Selector selFluidNodes =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(fluidSide->currentPartVec_);
            const stk::mesh::BucketVector& nodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK, selFluidNodes);

            // gather current p' in stable bucket/node order
            std::vector<scalar> pCurrent;
            for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
            {
                stk::mesh::Bucket& bucket = **ib;
                for (stk::mesh::Bucket::size_type i = 0; i < bucket.size(); ++i)
                {
                    pCurrent.push_back(
                        *stk::mesh::field_data(pCorrSTKFieldRef, bucket[i]));
                }
            }

            const label interfIdx = interf->index();
            auto& prev = fsiInterfacePCorrPrev_[interfIdx];
            if (prev.empty() || prev.size() != pCurrent.size())
            {
                prev = pCurrent;
                continue; // no prior snapshot to diff against on this call
            }

            for (size_t i = 0; i < pCurrent.size(); ++i)
            {
                const scalar r = pCurrent[i] - prev[i];
                sumSq += r * r;
                count++;
            }

            prev = pCurrent;
        }
    }

    messager::sumReduce(sumSq);
    messager::sumReduce(count);

    messager::sumReduce(pCorrFluxSumSq);
    messager::sumReduce(pCorrFluxCount);
    robinPCorrFluxRms_ =
        (pCorrFluxCount > 0)
            ? std::sqrt(pCorrFluxSumSq / static_cast<scalar>(pCorrFluxCount))
            : 0.0;

    interfacePCorrChangeNorm_ =
        (count > 0) ? std::sqrt(sumSq / static_cast<scalar>(count)) : 0.0;

    // Set reference on k=2 (first iteration that produces a nonzero
    // change; k=1 has no previous snapshot).
    const label fsiIter = sim_.controlsRef().iter;
    if (fsiIter == 2 && interfacePCorrChangeNorm_ > SMALL)
        interfacePCorrChangeNormRef_ = interfacePCorrChangeNorm_;

    interfacePCorrChangeNorm_ /=
        (interfacePCorrChangeNormRef_ > SMALL ? interfacePCorrChangeNormRef_
                                              : 1.0);
}

void fsiConvergence::updateForceBalance_()
{
    // fluid: p*n*dA + tau*dA (current area); solid: transferred traction * dA0

    // the flow equations derive from flowModel (pRef/wallShearStressRef)
    flowModel* flowModelPtr = nullptr;
    for (auto& eq : sim_.equationVector_)
    {
        flowModelPtr = dynamic_cast<flowModel*>(eq.get());
        if (flowModelPtr)
        {
            break;
        }
    }
    if (!flowModelPtr)
    {
        return;
    }

    solidDisplacementEquation* solidEq = nullptr;
    for (auto& equation : sim_.equationVector_)
    {
        if (equation->getID() == equationID::solidDisplacement)
        {
            solidEq = dynamic_cast<solidDisplacementEquation*>(equation.get());
            break;
        }
    }
    if (!solidEq)
    {
        return;
    }

    const auto& pSTKFieldRef = flowModelPtr->pRef().stkFieldRef();
    const auto& wallShearStressSTKFieldRef =
        flowModelPtr->wallShearStressRef().stkFieldRef();
    const bool isPShifted = flowModelPtr->pRef().isShifted();

    auto& tractionSTKField = solidEq->DRef().sideFluxFieldRef().stkFieldRef();
    stk::mesh::MetaData& metaData = tractionSTKField.mesh_meta_data();
    stk::mesh::BulkData& bulkData = tractionSTKField.get_mesh();

    const STKScalarField* exposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::exposed_area_vector_ID);
    const STKScalarField* origExposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::original_exposed_area_vector_ID);

    std::vector<scalar> fluidForceSum(SPATIAL_DIM, 0.0);
    std::vector<scalar> solidForceSum(SPATIAL_DIM, 0.0);

    std::unordered_set<label> visited;

    std::vector<scalar> ws_p;
    std::vector<scalar> ws_face_shape_function;

    for (const auto& domain : sim_.domainVector_)
    {
        for (const interface* interf : domain->interfacesRef())
        {
            if (!interf->isFluidSolidType())
            {
                continue;
            }
            if (!visited.insert(interf->index()).second)
            {
                continue;
            }

            const label masterIdx = interf->masterZoneIndex();
            const label slaveIdx = interf->slaveZoneIndex();
            const label fluidZoneIndex =
                (sim_.domainRef(masterIdx).type() == domainType::fluid)
                    ? masterIdx
                    : slaveIdx;
            const label solidZoneIndex =
                (fluidZoneIndex == masterIdx) ? slaveIdx : masterIdx;

            // fluid side
            {
                const interfaceSideInfo* fluidSide =
                    interf->interfaceSideInfoPtr(fluidZoneIndex);

                stk::mesh::Selector selOwnedSides =
                    metaData.locally_owned_part() &
                    stk::mesh::selectUnion(fluidSide->currentPartVec_);
                const stk::mesh::BucketVector& sideBuckets =
                    bulkData.get_buckets(metaData.side_rank(), selOwnedSides);

                for (const stk::mesh::Bucket* bucket : sideBuckets)
                {
                    MasterElement* meFC =
                        MasterElementRepo::get_surface_master_element(
                            bucket->topology());
                    const label nodesPerFace = meFC->nodesPerElement_;
                    const label numScsBip = meFC->numIntPoints_;

                    ws_p.resize(nodesPerFace);
                    ws_face_shape_function.resize(numScsBip * nodesPerFace);
                    scalar* p_p = &ws_p[0];
                    scalar* p_face_shape_function = &ws_face_shape_function[0];

                    if (isPShifted)
                    {
                        meFC->shifted_shape_fcn(&p_face_shape_function[0]);
                    }
                    else
                    {
                        meFC->shape_fcn(&p_face_shape_function[0]);
                    }

                    for (const stk::mesh::Entity side : *bucket)
                    {
                        stk::mesh::Entity const* sideNodeRels =
                            bulkData.begin_nodes(side);

                        for (label ni = 0; ni < nodesPerFace; ++ni)
                        {
                            stk::mesh::Entity node = sideNodeRels[ni];
                            p_p[ni] =
                                *stk::mesh::field_data(pSTKFieldRef, node);
                        }

                        const scalar* areaVec = stk::mesh::field_data(
                            *exposedAreaVecSTKFieldPtr, side);
                        const scalar* wallShearStressBip =
                            stk::mesh::field_data(wallShearStressSTKFieldRef,
                                                  side);

                        for (label ip = 0; ip < numScsBip; ++ip)
                        {
                            const label offSetAreaVec = ip * SPATIAL_DIM;
                            const label offSetSF = ip * nodesPerFace;

                            scalar pBip = 0.0;
                            for (label ic = 0; ic < nodesPerFace; ++ic)
                            {
                                pBip += p_face_shape_function[offSetSF + ic] *
                                        p_p[ic];
                            }

                            scalar aMag = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj = areaVec[offSetAreaVec + j];
                                aMag += axj * axj;
                            }
                            aMag = std::sqrt(aMag);

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                fluidForceSum[j] +=
                                    pBip * areaVec[offSetAreaVec + j] +
                                    wallShearStressBip[offSetAreaVec + j] *
                                        aMag;
                            }
                        }
                    }
                }
            }

            // solid side
            {
                const interfaceSideInfo* solidSide =
                    interf->interfaceSideInfoPtr(solidZoneIndex);

                stk::mesh::Selector selOwnedSides =
                    metaData.locally_owned_part() &
                    stk::mesh::selectUnion(solidSide->currentPartVec_);
                const stk::mesh::BucketVector& sideBuckets =
                    bulkData.get_buckets(metaData.side_rank(), selOwnedSides);

                for (const stk::mesh::Bucket* bucket : sideBuckets)
                {
                    MasterElement* meFC =
                        MasterElementRepo::get_surface_master_element(
                            bucket->topology());
                    const label numScsBip = meFC->numIntPoints_;

                    for (const stk::mesh::Entity side : *bucket)
                    {
                        const scalar* tractionValues =
                            stk::mesh::field_data(tractionSTKField, side);
                        const scalar* origAreaVec = stk::mesh::field_data(
                            *origExposedAreaVecSTKFieldPtr, side);

                        for (label ip = 0; ip < numScsBip; ++ip)
                        {
                            // traction stored per reference area: force =
                            // traction * dA0_mag
                            scalar dA0_mag = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    origAreaVec[SPATIAL_DIM * ip + j];
                                dA0_mag += axj * axj;
                            }
                            dA0_mag = std::sqrt(dA0_mag);

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                solidForceSum[j] +=
                                    tractionValues[SPATIAL_DIM * ip + j] *
                                    dA0_mag;
                            }
                        }
                    }
                }
            }
        }
    }

    messager::sumReduce(fluidForceSum);
    messager::sumReduce(solidForceSum);

    if (!forceBalanceStream_)
    {
        initializeForceBalanceFile_();
    }

    writeForceBalanceRow_(fluidForceSum, solidForceSum);
}

void fsiConvergence::initializeForceBalanceFile_()
{
    if (!messager::master())
    {
        return;
    }

    const fs::path filePath =
        sim_.getPostProcessingDirectory() / "fsi_force_balance.out";

    forceBalanceStream_ =
        std::make_shared<std::ofstream>(filePath, std::ios::app);
    assert(forceBalanceStream_->is_open());

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    auto& fout = *forceBalanceStream_;
    fout << "\n# ===== New run =====\n";
    fout << "# Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << "# Version: " << accel::PROJECT_VERSION << '\n';
    fout << "# Git hash: " << accel::GIT_HASH << '\n';
    fout << "# Git describe: " << accel::GIT_DESCRIBE << '\n';
    fout << "# FSI interface force balance (traction integrated over "
            "exposed area, per outer iteration)\n";
    fout << "# \n";
    fout << "# " << "sim_time[s]" << '\t' << "fsi_iter";
    for (label d = 0; d < SPATIAL_DIM; ++d)
        fout << '\t' << "fluid_force_" << ("xyz"[d]);
    for (label d = 0; d < SPATIAL_DIM; ++d)
        fout << '\t' << "solid_force_" << ("xyz"[d]);
    for (label d = 0; d < SPATIAL_DIM; ++d)
        fout << '\t' << "diff_force_" << ("xyz"[d]);
    for (label d = 0; d < SPATIAL_DIM; ++d)
        fout << '\t' << "%_diff_" << ("xyz"[d]);
    fout << '\n';
}

void fsiConvergence::writeForceBalanceRow_(
    const std::vector<scalar>& fluidForce,
    const std::vector<scalar>& solidForce)
{
    if (!messager::master() || !forceBalanceStream_)
        return;

    std::vector<scalar> diff(SPATIAL_DIM);
    std::vector<scalar> pctDiff(SPATIAL_DIM);
    for (label j = 0; j < SPATIAL_DIM; ++j)
    {
        diff[j] = solidForce[j] - fluidForce[j];
        pctDiff[j] = (std::abs(fluidForce[j]) > SMALL)
                         ? 100.0 * diff[j] / fluidForce[j]
                         : 0.0;
    }

    auto& fout = *forceBalanceStream_;
    fout << std::setprecision(6) << std::scientific << sim_.getSimulationTime()
         << '\t' << sim_.getIterationCount();
    for (label j = 0; j < SPATIAL_DIM; ++j)
        fout << '\t' << fluidForce[j];
    for (label j = 0; j < SPATIAL_DIM; ++j)
        fout << '\t' << solidForce[j];
    for (label j = 0; j < SPATIAL_DIM; ++j)
        fout << '\t' << diff[j];
    for (label j = 0; j < SPATIAL_DIM; ++j)
        fout << '\t' << pctDiff[j];
    fout << std::endl;
}

bool fsiConvergence::isConverged() const
{
    const auto& fsiConv =
        sim_.controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceCriteria_.fsiConvergence_;
    if (!fsiConv.enabled_)
    {
        return true;
    }

    // defensive: never converged without a criterion
    if (fsiConv.criteria_.empty())
    {
        return false;
    }

    bool converged = true;
    for (const auto& criterion : fsiConv.criteria_)
    {
        switch (criterion)
        {
            case fsiConvergenceType::fsiInterfaceResidual:
                if (fsiInterfaceResidualNorms_.empty())
                {
                    converged = false;
                    break;
                }
                // iteration 1 reports the real ||X^1|| but may not converge:
                // the fluid has not yet seen the new interface displacement
                converged = converged && (sim_.getIterationCount() > 1) &&
                            (fsiInterfaceResidualNorm_ <=
                             fsiConv.fsiInterfaceResidualTarget_);
                break;
            case fsiConvergenceType::fsiForceResidual:
                if (fsiForceResidualNorms_.empty())
                {
                    converged = false;
                    break;
                }
                converged =
                    converged && (sim_.getIterationCount() > 1) &&
                    (fsiForceResidualNorm_ <= fsiConv.fsiForceResidualTarget_);
                break;
            default:
                converged = false;
                break;
        }
    }

    return converged;
}

void fsiConvergence::initializeResidualFile_(label interfIdx,
                                             const std::string& interfName,
                                             const std::string& criterionName)
{
    if (!messager::master())
    {
        return;
    }

    std::string baseName = criterionName + "_" + interfName;
    std::replace(baseName.begin(), baseName.end(), ' ', '_');

    const fs::path filePath = sim_.getResidualDirectory() / (baseName + ".out");

    auto stream = std::make_shared<std::ofstream>(filePath, std::ios::app);
    assert(stream->is_open());

    if (fs::exists(filePath) && fs::file_size(filePath) > 0)
    {
        residualStreams_[criterionName][interfIdx] = stream;
        return;
    }

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    auto& fout = *stream;
    fout << "# Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << "# Version: " << accel::PROJECT_VERSION << '\n';
    fout << "# Git hash: " << accel::GIT_HASH << '\n';
    fout << "# Git describe: " << accel::GIT_DESCRIBE << '\n';
    fout << "# FSI convergence residual history - interface: " << interfName
         << '\n';
    fout << "# Criterion: " << criterionName << '\n';
    fout << "# \n";
    fout << "# " << "global_iterations" << '\t' << "inner_iterations" << '\t'
         << "sim_time[s]" << '\t' << "residual_norm" << '\n';

    residualStreams_[criterionName][interfIdx] = stream;
}

void fsiConvergence::writeResidualLine_(const std::string& criterionName,
                                        label interfIdx,
                                        scalar residualNorm)
{
    if (!messager::master())
    {
        return;
    }

    auto it = residualStreams_.find(criterionName);
    if (it == residualStreams_.end())
    {
        return;
    }

    auto& streams = it->second;
    auto streamIt = streams.find(interfIdx);
    if (streamIt == streams.end())
    {
        return;
    }

    auto& fout = *(streamIt->second);
    fout << sim_.getGlobalIterationCount() << '\t' << sim_.getIterationCount()
         << '\t' << std::setprecision(3) << std::scientific
         << sim_.getSimulationTime() << '\t' << std::setprecision(6)
         << std::scientific << residualNorm << std::endl;
}

void fsiConvergence::initializeDiagnosticsFile_()
{
    if (!messager::master())
    {
        return;
    }

    const fs::path filePath =
        sim_.getPostProcessingDirectory() / "fsi_diagnostics.out";

    diagnosticsStream_ =
        std::make_shared<std::ofstream>(filePath, std::ios::app);
    assert(diagnosticsStream_->is_open());

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    auto& fout = *diagnosticsStream_;
    fout << "\n# ===== New run =====\n";
    fout << "# Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << "# Version: " << accel::PROJECT_VERSION << '\n';
    fout << "# Git hash: " << accel::GIT_HASH << '\n';
    fout << "# Git describe: " << accel::GIT_DESCRIBE << '\n';
    fout << "# Combined FSI outer-iteration diagnostics\n";
    fout << "# Coupling method: "
         << (isRobinNeumann_ ? (rnFallback_ ? "dirichlet_neumann_fallback"
                                            : "robin_neumann")
                             : "dirichlet_neumann")
         << '\n';
    fout << "# \n";
    fout << "# " << "sim_time[s]" << '\t' << "fsi_iter" << '\t'
         << "kinematic_res" << '\t' << "force_res" << '\t'
         << "interface_pCorr_change" << '\t' << "robin_accelTerm_rms" << '\t'
         << "robin_pCorr_flux_rms";
    fout << '\n';
}

void fsiConvergence::writeDiagnosticsRow_(label fsiIter,
                                          scalar kinematicRes,
                                          scalar forceRes,
                                          scalar interfacePCorrChange,
                                          scalar robinAccelTermRms,
                                          scalar robinPCorrFluxRms)
{
    if (!messager::master() || !diagnosticsStream_)
        return;

    auto& fout = *diagnosticsStream_;
    fout << std::setprecision(6) << std::scientific << sim_.getSimulationTime()
         << '\t' << fsiIter << '\t' << kinematicRes << '\t' << forceRes << '\t'
         << interfacePCorrChange << '\t' << robinAccelTermRms << '\t'
         << robinPCorrFluxRms;

    fout << std::endl;
}

} // namespace accel
