// File       : physicsConvergence.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Physics-based convergence checks for coupled simulations
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "physicsConvergence.h"
#include "git_revision.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#include "messager.h"
#include "simulation.h"
#include "solidDisplacementEquation.h"

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
    fsiInterfaceResidualNorms_.clear();
    fsiInterfaceResidualNorm_ = 0.0;
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
                break;
            default:
                break;
        }
    }
}

void physicsConvergence::updateFsiInterfaceResidual_(bool writeResiduals)
{
#ifdef HAS_INTERFACE
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
                fsiInterfaceResidualNorms_[interfIdx] = 1.0;
                fsiInterfaceResidualNorm_ =
                    std::max(fsiInterfaceResidualNorm_, 1.0);
                continue;
            }

            scalar normSq = 0.0;
            for (size_t i = 0; i < vecSize; ++i)
            {
                const scalar r = DCurrent[i] - prev[i];
                normSq += r * r;
            }
            messager::sumReduce(normSq);
            const scalar norm = std::sqrt(normSq);

            auto& maxNorm = fsiInterfaceResidualNormMax_[interfIdx];
            maxNorm = std::max(maxNorm, norm);
            const scalar normRel = norm / (maxNorm + SMALL);
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
#endif /* HAS_INTERFACE */
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
                converged = false;
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
    fout << "# Accel solver timestamp: "
         << std::put_time(std::localtime(&in_time_t), "%c\n");
    fout << "# Git revision: " << accel::git_revision << '\n';
    fout << "# Physics convergence residual history — interface: " << interfName
         << '\n';
    fout << "# Criterion: " << criterionName << '\n';
    fout << "# \n";
    fout << "# " << "global_iterations" << '\t' << "inner_iterations" << '\t'
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
