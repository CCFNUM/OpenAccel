// File       : robinNeumannState.h
// Created    : Thu Oct 01 2026 20:14:34 (+0300)
// Author     : Mohammad Nasreddine
// Description: Robin-Neumann FSI pressure BC state shared by the single-phase
//              and VOF pressure-correction assemblers
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef ROBINNEUMANNSTATE_H
#define ROBINNEUMANNSTATE_H

// code
#include "types.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <unordered_map>

namespace accel
{

class domain;
class interfaceSideInfo;

// phase of the wall kernel; the defaults give the single-phase kernel
struct wallPhase
{
    const STKScalarField* rho = nullptr;    // phase density
    const STKScalarField* alpha = nullptr;  // phase volume fraction (VOF)
    const STKScalarField* rhoMix = nullptr; // mixture density in c1 (VOF)
    scalar densityScale = 1.0;
};

// wall kernels of the pressure-correction assemblers (robinNeumannWall.h)
template <class Assembler>
void assembleWallInterfaceSide(Assembler& self,
                               const domain* domain,
                               const interfaceSideInfo* interfaceSideInfoPtr,
                               typename Assembler::Context* ctx,
                               const wallPhase& phase,
                               const bool robin);

template <class Assembler>
void accumulateRnWallFluxSide(Assembler& self,
                              const domain* domain,
                              const interfaceSideInfo* interfaceSideInfoPtr);

template <class Assembler>
void accumulateRnWallFluxSides(Assembler& self, const domain* domain);

class robinNeumannState
{
public:
    // clear the parameters and warnings of the last step (preTimeStep)
    void invalidateRobinParams()
    {
        robinHs_ = 0.0;
        robinRhoSolid_ = 0.0;
        robinParamsValid_ = false;
        rnWarned_.clear();
        robinAccelTerm_.clear();
    }

    // h_s = a_p*dt and rho_s, set once per timestep (Tukovic et al. 2014)
    void setRobinParams(scalar hs, scalar rhoSolid)
    {
        robinHs_ = hs;
        robinRhoSolid_ = rhoSolid;
        robinParamsValid_ = true;
    }

    bool robinParamsValid() const
    {
        return robinParamsValid_;
    }

    scalar robinHs() const
    {
        return robinHs_;
    }

    scalar robinRhoSolid() const
    {
        return robinRhoSolid_;
    }

    // E is stored at the first corrector; a flag, since E = 0 is a valid value
    void requestFreezeE()
    {
        freezeE_ = true;
    }

    void clearFreezeE()
    {
        freezeE_ = false;
    }

    bool freezeE() const
    {
        return freezeE_;
    }

    // true at the first call per interface and timestep
    bool shouldWarnRn(label interfaceIndex)
    {
        return rnWarned_.insert(interfaceIndex).second;
    }

    // per-BIP rho_s*h_s*a_n, key (side id << 8) | ip; read by fsiConvergence
    const std::unordered_map<uint64_t, scalar>& robinAccelTermMap() const
    {
        return robinAccelTerm_;
    }

    void setRobinAccelTerm(uint64_t key, scalar accelTerm)
    {
        robinAccelTerm_[key] = accelTerm;
    }

    // c1 statistics of a timestep, area weighted over all assembly passes
    void resetC1Stats()
    {
        robinC1Min_ = std::numeric_limits<scalar>::max();
        robinC1Max_ = 0.0;
        robinC1WtSum_ = 0.0;
        robinAreaSum_ = 0.0;
    }

    void recordC1(scalar c1val, scalar area)
    {
        robinC1Min_ = std::min(robinC1Min_, c1val);
        robinC1Max_ = std::max(robinC1Max_, c1val);
        robinC1WtSum_ += c1val * area;
        robinAreaSum_ += area;
    }

    // Rank-local values. Caller must perform the global MPI reduction.
    scalar robinC1MinLocal() const
    {
        return robinC1Min_;
    }

    scalar robinC1MaxLocal() const
    {
        return robinC1Max_;
    }

    scalar robinC1WtSumLocal() const
    {
        return robinC1WtSum_;
    }

    scalar robinC1AreaSumLocal() const
    {
        return robinAreaSum_;
    }

private:
    // Robin-Neumann FSI pressure BC state (see setRobinParams())
    bool robinParamsValid_ = false;
    scalar robinHs_ = 0.0;
    scalar robinRhoSolid_ = 0.0;

    // compute and store E in the next assembly pass (see requestFreezeE())
    bool freezeE_ = true;

    // interfaces already warned about invalid RN params this timestep
    std::set<label> rnWarned_;

    // c1 statistics (rank-local; see resetC1Stats())
    scalar robinC1Min_ = std::numeric_limits<scalar>::max();
    scalar robinC1Max_ = 0.0;
    scalar robinC1WtSum_ = 0.0;
    scalar robinAreaSum_ = 0.0;

    std::unordered_map<uint64_t, scalar> robinAccelTerm_;
};

// rho_f in c1: local BIP density unless reference_fluid_density is set
inline scalar robinFluidDensity(const bool useReferenceFluidDensity,
                                const scalar referenceFluidDensity,
                                const scalar rhoLocal)
{
    return useReferenceFluidDensity
               ? referenceFluidDensity
               : std::max(rhoLocal, static_cast<scalar>(1e-10));
}

// BIP coefficients shared by the RN kernels and accumulateRnWallFlux()
struct robinBipCoefficients
{
    scalar aMag = 0.0;    // |A| of the exposed area vector at the BIP
    scalar invAMag = 0.0; // 1/(|A| + SMALL)
    scalar du_n = 0.0;    // sum_j du_j n_j^2
    scalar c1 = 0.0;      // rho_s*h_s/rho_f
    scalar deltaBn = 0.0; // |(x_bip - x_nearest).n|
};

inline robinBipCoefficients
computeRobinBipCoefficients(const scalar* areaVecIp,    // [SPATIAL_DIM]
                            const scalar* duBipVec,     // [SPATIAL_DIM]
                            const scalar* coordBip,     // [SPATIAL_DIM]
                            const scalar* coordNearest, // [SPATIAL_DIM]
                            const scalar rhoFluid,
                            const scalar robinRhoSolid,
                            const scalar robinHs)
{
    robinBipCoefficients rc;

    scalar aMagSq = 0.0;
    for (label j = 0; j < SPATIAL_DIM; ++j)
    {
        const scalar axj = areaVecIp[j];
        aMagSq += axj * axj;
    }
    rc.aMag = std::sqrt(aMagSq);
    rc.invAMag = 1.0 / (rc.aMag + SMALL);

    // normal component of the diagonal momentum coefficient
    for (label j = 0; j < SPATIAL_DIM; ++j)
    {
        const scalar nj = areaVecIp[j] * rc.invAMag;
        rc.du_n += duBipVec[j] * nj * nj;
    }

    rc.c1 = robinRhoSolid * robinHs / rhoFluid;

    // normal distance nearestNode-BIP (Tukovic Fig. 2), ~0 on flat faces
    for (label j = 0; j < SPATIAL_DIM; ++j)
    {
        const scalar dx = coordBip[j] - coordNearest[j];
        rc.deltaBn += dx * areaVecIp[j] * rc.invAMag;
    }
    rc.deltaBn = std::abs(rc.deltaBn);

    return rc;
}

} // namespace accel

#endif // ROBINNEUMANNSTATE_H
