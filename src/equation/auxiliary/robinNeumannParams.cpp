// File       : robinNeumannParams.cpp
// Created    : Thu Oct 01 2026 20:14:34 (+0300)
// Author     : Mohammad Nasreddine
// Description: Robin-Neumann FSI pressure BC helpers shared by the
//              single-phase and VOF pressure-correction equations
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "robinNeumannParams.h"
#include "domain.h"
#include "interface.h"
#include "messager.h"
#include "robinNeumannState.h"
#include "simulation.h"
#include "solidDisplacementEquation.h"
#include "version.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>

namespace accel
{

bool anyRobinNeumannFsi(
    const std::vector<std::shared_ptr<domain>>& domainVector)
{
    for (const auto& domainPtr : domainVector)
    {
        for (const interface* interf : domainPtr->interfacesRef())
        {
            if (interf->isFluidSolidType() && interf->isRobinNeumann())
                return true;
        }
    }
    return false;
}

robinNeumannParams computeRobinNeumannParams(
    simulation& sim,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    scalar dt)
{
    robinNeumannParams result;

    if (!anyRobinNeumannFsi(domainVector))
    {
        return result;
    }

    auto* solidEq = dynamic_cast<solidDisplacementEquation*>(
        sim.equationPtr(equationID::solidDisplacement));

    scalar E_val = 0.0;
    scalar nu_val = 0.0;
    scalar rhoS_val = 0.0;
    bool hasLocalValues = false;

    if (solidEq)
    {
        const stk::mesh::MetaData& metaData = solidEq->meshRef().metaDataRef();
        const stk::mesh::BulkData& bulkData = solidEq->meshRef().bulkDataRef();

        // scope selector to the solid domains' zone parts
        stk::mesh::PartVector solidParts;
        for (const auto& domainPtr : sim.domainVector())
        {
            if (domainPtr->type() == domainType::solid)
            {
                const auto& parts = domainPtr->zonePtr()->interiorParts();
                solidParts.insert(solidParts.end(), parts.begin(), parts.end());
            }
        }

        stk::mesh::Selector selSolidOwned =
            metaData.locally_owned_part() & stk::mesh::selectUnion(solidParts);

        const stk::mesh::BucketVector& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selSolidOwned);

        for (const stk::mesh::Bucket* bucket : nodeBuckets)
        {
            if (bucket->size() > 0)
            {
                stk::mesh::Entity node = (*bucket)[0];
                const scalar* Eptr =
                    stk::mesh::field_data(solidEq->ERef().stkFieldRef(), node);
                const scalar* nuPtr =
                    stk::mesh::field_data(solidEq->nuRef().stkFieldRef(), node);
                const scalar* rhoPtr = stk::mesh::field_data(
                    solidEq->rhoRef().stkFieldRef(), node);
                if (Eptr && nuPtr && rhoPtr)
                {
                    E_val = *Eptr;
                    nu_val = *nuPtr;
                    rhoS_val = *rhoPtr;
                    hasLocalValues = true;
                }
                break;
            }
        }
    }

    // take the (E, nu, rho_s) triplet of one node, from the lowest rank that
    // owns a solid node, so the three values always belong together
    int ownerRank = hasLocalValues ? messager::myProcNo() : messager::nProcs();
    messager::minReduce(ownerRank);
    if (ownerRank < messager::nProcs())
    {
        messager::broadcast(E_val, ownerRank);
        messager::broadcast(nu_val, ownerRank);
        messager::broadcast(rhoS_val, ownerRank);
    }

    if (E_val > 0.0 && rhoS_val > 0.0)
    {
        const scalar mu_s = E_val / (2.0 * (1.0 + nu_val));
        const scalar lambda_s =
            E_val * nu_val / ((1.0 + nu_val) * (1.0 - 2.0 * nu_val));
        const scalar a_p = std::sqrt((lambda_s + 2.0 * mu_s) / rhoS_val);
        scalar h_s = a_p * dt;

        // physical shell thickness (Tukovic Sec. 2.4)
        bool thicknessFound = false;
        for (const auto& domainPtr : domainVector)
        {
            for (const interface* interf : domainPtr->interfacesRef())
            {
                if (interf->isFluidSolidType() && interf->isRobinNeumann() &&
                    interf->usePhysicalSolidThickness())
                {
                    h_s = interf->solidThickness();
                    thicknessFound = true;
                    break;
                }
            }
            if (thicknessFound)
                break;
        }

        result.hs = h_s;
        result.rhoSolid = rhoS_val;
        result.valid = true;
    }

    return result;
}

void writeRobinC1Diagnostics(
    const robinNeumannState& assemblerState,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    std::shared_ptr<std::ofstream>& stream,
    scalar simTime)
{
    // input-only decision: identical on all ranks (reductions stay matched)
    bool writeRequested = false;
    for (const auto& domainPtr : domainVector)
    {
        for (const interface* interf : domainPtr->interfacesRef())
        {
            if (interf->isFluidSolidType() && interf->isRobinNeumann() &&
                interf->writeRobinAlpha())
            {
                writeRequested = true;
                break;
            }
        }
        if (writeRequested)
            break;
    }

    if (!writeRequested || domainVector.empty())
    {
        return;
    }

    // global reduction (collective, all ranks)
    scalar c1Min = assemblerState.robinC1MinLocal();
    scalar c1Max = assemblerState.robinC1MaxLocal();
    scalar c1WtSum = assemblerState.robinC1WtSumLocal();
    scalar c1AreaSum = assemblerState.robinC1AreaSumLocal();

    messager::minReduce(c1Min);
    messager::maxReduce(c1Max);
    messager::sumReduce(c1WtSum);
    messager::sumReduce(c1AreaSum);

    // no RN BIP assembled on any rank: min is still the sentinel
    if (c1Min == std::numeric_limits<scalar>::max())
    {
        c1Min = 0.0;
    }

    const scalar c1Mean = (c1AreaSum > 0.0) ? c1WtSum / c1AreaSum : 0.0;

    if (!messager::master())
    {
        return;
    }

    if (!stream)
    {
        simulation& sim = domainVector[0]->simulationRef();
        const fs::path filePath =
            sim.getPostProcessingDirectory() / "robin_alpha.out";

        stream = std::make_shared<std::ofstream>(filePath, std::ios::app);
        assert(stream->is_open());

        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        auto& fout = *stream;
        fout << "\n# ===== New run =====\n";
        fout << "# Accel solver timestamp: "
             << std::put_time(std::localtime(&in_time_t), "%c\n");
        fout << "# Version: " << accel::PROJECT_VERSION << '\n';
        fout << "# Git hash: " << accel::GIT_HASH << '\n';
        fout << "# Git describe: " << accel::GIT_DESCRIBE << '\n';
        fout << "# Robin alpha scale (c1 = rho_s*h_s/rho_f) statistics over "
                "all RN BIPs of the timestep\n";
        fout << "# (min/max over all assembly passes; mean area-weighted and "
                "pass-count weighted)\n";
        fout << "# Coupling method: robin_neumann\n";
        fout << "# \n";
        fout << "# " << "sim_time[s]" << '\t' << "c1_min" << '\t' << "c1_max"
             << '\t' << "c1_mean" << '\n';
    }

    auto& fout = *stream;
    fout << std::setprecision(6) << std::scientific << simTime << '\t' << c1Min
         << '\t' << c1Max << '\t' << c1Mean << std::endl;
}

} // namespace accel
