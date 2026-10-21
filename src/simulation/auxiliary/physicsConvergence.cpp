// File       : physicsConvergence.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Physics-based convergence checks for coupled simulations
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "physicsConvergence.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#include "messager.h"
#include "simulation.h"
#include "solidDisplacementEquation.h"
#include "version.h"

namespace accel
{

physicsConvergence::physicsConvergence(simulation& sim) : sim_(sim)
{
}

bool physicsConvergence::enabled() const
{
    const auto& physConv = sim_.controlsRef()
                               .solverRef()
                               .solverControl_.basicSettings_
                               .convergenceCriteria_.physicsConvergence_;
    return physConv.enabled_;
}

void physicsConvergence::resetForTimeStep()
{
    if (!enabled())
    {
        return;
    }

    fsiInterfaceDispPrev_.clear();
    fsiInterfaceResidualNormMax_.clear();
    fsiInterfaceMaxTotalDisplNorm_.clear();
    fsiInterfaceResidualNorms_.clear();
    fsiInterfaceResidualNorm_ = 0.0;

    fsiInterfaceTractionPrev_.clear();
    fsiInterfaceTractionResidualNormMax_.clear();
    fsiForceResidualNorms_.clear();
    fsiForceResidualNorm_ = 0.0;
}

void physicsConvergence::update()
{
    if (!enabled())
    {
        return;
    }

    const auto& physConv = sim_.controlsRef()
                               .solverRef()
                               .solverControl_.basicSettings_
                               .convergenceCriteria_.physicsConvergence_;

    fsiInterfaceResidualNorm_ = 0.0;
    fsiForceResidualNorm_ = 0.0;

    for (const auto& criterion : physConv.criteria_)
    {
        switch (criterion)
        {
            case physicsConvergenceType::fsiInterfaceResidual:
                updateFsiInterfaceResidual_(physConv.writeResiduals_);
                break;
            case physicsConvergenceType::fsiForceResidual:
                updateFsiForceResidual_(physConv.writeResiduals_);
                break;
            default:
                break;
        }
    }
}

void physicsConvergence::updateFsiInterfaceResidual_(bool writeResiduals)
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

            stk::mesh::Selector selFluidNodes =
                DField.mesh_meta_data().universal_part() &
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
            {
                size_t offset = 0;
                for (auto ib = fluidNodeBuckets.begin();
                     ib != fluidNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& b = **ib;
                    const scalar* Db = stk::mesh::field_data(DField, b);
                    for (size_t iNode = 0; iNode < b.size(); ++iNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            DCurrent[offset++] = Db[SPATIAL_DIM * iNode + i];
                        }
                    }
                }
            }

            const label interfIdx = interf->index();
            auto& prev = fsiInterfaceDispPrev_[interfIdx];
            if (prev.empty() || prev.size() != vecSize)
            {
                prev = DCurrent;

                // seed running max of ||D_total|| (avoids dividing by 0)
                scalar seedTotalDisplNormSq = 0.0;
                for (size_t i = 0; i < vecSize; ++i)
                {
                    seedTotalDisplNormSq += DCurrent[i] * DCurrent[i];
                }
                messager::sumReduce(seedTotalDisplNormSq);
                auto& maxTotalDispl = fsiInterfaceMaxTotalDisplNorm_[interfIdx];
                maxTotalDispl =
                    std::max(maxTotalDispl, std::sqrt(seedTotalDisplNormSq));

                fsiInterfaceResidualNorms_[interfIdx] = 1.0;
                fsiInterfaceResidualNorm_ =
                    std::max(fsiInterfaceResidualNorm_, 1.0);
                continue;
            }

            scalar normSq = 0.0;
            scalar totalDisplNormSq = 0.0;
            for (size_t i = 0; i < vecSize; ++i)
            {
                const scalar r = DCurrent[i] - prev[i];
                normSq += r * r;
                totalDisplNormSq += DCurrent[i] * DCurrent[i];
            }
            messager::sumReduce(normSq);
            messager::sumReduce(totalDisplNormSq);
            const scalar norm = std::sqrt(normSq);
            const scalar totalDisplNorm = std::sqrt(totalDisplNormSq);

            // norm 1: mismatch over its running max
            auto& maxNorm = fsiInterfaceResidualNormMax_[interfIdx];
            maxNorm = std::max(maxNorm, norm);
            const scalar normRel1 = norm / (maxNorm + SMALL);

            // norm 2: mismatch over max ||D_total||; take min of both (s4f)
            auto& maxTotalDispl = fsiInterfaceMaxTotalDisplNorm_[interfIdx];
            maxTotalDispl = std::max(maxTotalDispl, totalDisplNorm);
            const scalar normRel2 = norm / (maxTotalDispl + SMALL);

            const scalar normRel = std::min(normRel1, normRel2);
            fsiInterfaceResidualNorms_[interfIdx] = normRel;
            fsiInterfaceResidualNorm_ =
                std::max(fsiInterfaceResidualNorm_, normRel);

            prev = DCurrent;

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

void physicsConvergence::updateFsiForceResidual_(bool writeResiduals)
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

    // fluid-side traction, side-rank IP data on D's side-flux field
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

            // locally-owned fluid-side interface sides
            stk::mesh::Selector selFluidSides =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(fluidSide->currentPartVec_);
            const stk::mesh::BucketVector& fluidSideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selFluidSides);

            // flatten traction in stable bucket/side/IP order
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
            if (prev.empty() || prev.size() != tCurrent.size())
            {
                prev = tCurrent;
                fsiForceResidualNorms_[interfIdx] = 1.0;
                fsiForceResidualNorm_ = std::max(fsiForceResidualNorm_, 1.0);
                continue;
            }

            scalar normSq = 0.0;
            scalar tractionNormSq = 0.0;
            for (size_t i = 0; i < tCurrent.size(); ++i)
            {
                const scalar r = tCurrent[i] - prev[i];
                normSq += r * r;
                tractionNormSq += tCurrent[i] * tCurrent[i];
            }
            messager::sumReduce(normSq);
            messager::sumReduce(tractionNormSq);
            const scalar norm = std::sqrt(normSq);
            const scalar tractionNorm = std::sqrt(tractionNormSq);

            auto& maxNorm = fsiInterfaceTractionResidualNormMax_[interfIdx];
            maxNorm = std::max(maxNorm, tractionNorm);
            const scalar normRel = norm / (maxNorm + SMALL);
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

bool physicsConvergence::isConverged() const
{
    const auto& physConv = sim_.controlsRef()
                               .solverRef()
                               .solverControl_.basicSettings_
                               .convergenceCriteria_.physicsConvergence_;
    if (!physConv.enabled_)
    {
        return true;
    }

    if (physConv.criteria_.empty())
    {
        return false;
    }

    bool converged = true;
    for (const auto& criterion : physConv.criteria_)
    {
        switch (criterion)
        {
            case physicsConvergenceType::fsiInterfaceResidual:
                if (fsiInterfaceResidualNorms_.empty())
                {
                    converged = false;
                    break;
                }
                converged = converged && (fsiInterfaceResidualNorm_ <=
                                          physConv.fsiInterfaceResidualTarget_);
                break;
            case physicsConvergenceType::fsiForceResidual:
                if (fsiForceResidualNorms_.empty())
                {
                    converged = false;
                    break;
                }
                converged = converged && (fsiForceResidualNorm_ <=
                                          physConv.fsiForceResidualTarget_);
                break;
            default:
                converged = false;
                break;
        }
    }

    return converged;
}

void physicsConvergence::initializeResidualFile_(
    label interfIdx,
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
    fout << COMMENT << "Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << COMMENT << "Version: " << accel::PROJECT_VERSION << '\n';
    fout << COMMENT << "Git hash: " << accel::GIT_HASH << '\n';
    fout << COMMENT << "Git describe: " << accel::GIT_DESCRIBE << '\n';
    fout << COMMENT
         << "Physics convergence residual history — interface: " << interfName
         << '\n';
    fout << COMMENT << "Criterion: " << criterionName << '\n';
    fout << COMMENT << '\n';
    fout << COMMENT << "global_iterations" << '\t' << "inner_iterations" << '\t'
         << "sim_time[s]" << '\t' << "residual_norm" << '\n';

    residualStreams_[criterionName][interfIdx] = stream;
}

void physicsConvergence::writeResidualLine_(const std::string& criterionName,
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

} // namespace accel
