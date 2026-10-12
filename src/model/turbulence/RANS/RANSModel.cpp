// File       : RANSModel.cpp
// Created    : Mon Mar 25 2024 16:48:19 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "RANSModel.h"
#include "messager.h"
#include "realm.h"
#include "zoneTransformation.h"

namespace accel
{

RANSModel::RANSModel(realm* realm) : turbulenceModel(realm)
{
}

void RANSModel::setupCurvatureCorrection(const std::shared_ptr<domain> domain)
{
    // f_r is a plain node field: allocate it now so that it exists when the
    // output file is set up (fields are matched by name in output_fields)
    frRef();

    if (!SijPtr_)
    {
        SijPtr_ =
            std::make_unique<nodeTensorField>(&this->realmRef(),
                                              "strain_rate_tensor",
                                              controlsRef().getNumberOfStates(),
                                              false);
        SijPtr_->setupGradientField();

        const auto& basic =
            controlsRef()
                .solverRef()
                .solverControl_.basicSettings_.interpolationSchemeType_;
        SijPtr_->setInterpolationScheme(basic.pressureInterpolationType_);
        SijPtr_->setGradientInterpolationScheme(
            basic.pressureGradientInterpolationType_);

        // S_ij is rebuilt every iteration: full gradient, no relaxation lag
        SijPtr_->setGradURF(1.0);
    }

    SijPtr_->setZone(domain->index());
}

void RANSModel::updateCurvatureCorrection(const std::shared_ptr<domain> domain)
{
    // Smirnov-Menter: full Lagrangian strain derivative in the
    // rotating reference frame and the Spalart-Shur production multiplier.
    // The velocity solver stores absolute velocity; convect with relative
    // velocity and retain the explicit rotating-frame strain contribution.
    auto& bulkData = meshRef().bulkDataRef();
    auto& metaData = meshRef().metaDataRef();
    const label iZone = domain->index();
    const scalar Cscale = domain->turbulence_.curvatureCorrectionCoeff_;
    const scalar beta = betaStar();
    const scalar invCmu = 1.0 / Cmu();

    // Embed a planar 2D tensor in 3D to retain out-of-plane rotation.
    static constexpr scalar eps[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}},
                                            {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
                                            {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

    scalar OmRot[3] = {0.0, 0.0, 0.0};
    scalar origin[3] = {0.0, 0.0, 0.0};
    const bool frameRotating = domain->zonePtr()->frameRotating();
    if (frameRotating)
    {
        const auto& rot = domain->zonePtr()->transformationRef().rotation();
#if SPATIAL_DIM == 2
        // Planar frames rotate about z (zoneTransformation convention).
        OmRot[2] = rot.omega_;
        origin[0] = rot.origin_[0];
        origin[1] = rot.origin_[1];
#else
        for (label m = 0; m < 3; ++m)
        {
            OmRot[m] = rot.omega_ * rot.axis_[m];
            origin[m] = rot.origin_[m];
        }
#endif
    }

    const auto& controls = controlsRef();
    const bool transient = controls.isTransient();
    const label timeStepCount = controls.getTimeStepCount();
    const bool hasOld = transient && timeStepCount > 1;
    const bool useBDF2 =
        transient && timeStepCount > 2 &&
        controls.solverRef().solverControl_.basicSettings_.transientScheme_ ==
            transientSchemeType::secondOrderBackwardEuler;

    scalar dtInv = 0.0;
    scalar timeCoeff[3] = {0.0, 0.0, 0.0};
    if (hasOld)
    {
        const scalar dt = controls.getTimestep();
        dtInv = 1.0 / dt;
        if (useBDF2)
        {
            const auto c = BDF2::coeff(dt, controls.getTimestep(-1));
            for (label n = 0; n < 3; ++n)
            {
                timeCoeff[n] = c[n];
            }
        }
        else
        {
            timeCoeff[0] = 1.0;
            timeCoeff[1] = -1.0;
        }
    }

    const bool kEpsilon =
        domain->turbulence_.option_ == turbulenceOption::kEpsilon;
    const auto* coordsField = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, getCoordinatesID_(domain));
    const auto* gradUField = URef().gradRef().stkFieldPtr();

    const stk::mesh::Selector selector =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selector);

    // S_ij and its gradient must use the current converged/iterated velocity
    // field; the STK time states hold the previous physical timestep(s).
    auto* strainField = SijPtr_->stkFieldPtr();
    const scalar* gradU = nullptr;
    const scalar* g = nullptr;
    scalar* strain = nullptr;
    scalar* Sflat = nullptr;
    for (const auto* bucket : buckets)
    {
        gradU = stk::mesh::field_data(*gradUField, *bucket);
        strain = stk::mesh::field_data(*strainField, *bucket);
        for (stk::mesh::Bucket::size_type node = 0; node < bucket->size();
             ++node)
        {
            g = gradU + SPATIAL_DIM * SPATIAL_DIM * node;
            Sflat = strain + SPATIAL_DIM * SPATIAL_DIM * node;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    Sflat[SPATIAL_DIM * i + j] =
                        0.5 * (g[SPATIAL_DIM * i + j] + g[SPATIAL_DIM * j + i]);
                }
            }
        }
    }
    SijPtr_->synchronizeGhostedEntities(iZone);
    SijPtr_->updateGradientField(iZone);

    const auto* UField = URef().stkFieldPtr();
    const auto* scaleField =
        kEpsilon ? epsilonRef().stkFieldPtr() : omegaRef().stkFieldPtr();
    const auto* kField = kEpsilon ? kRef().stkFieldPtr() : nullptr;
    const auto* gradStrainField = SijPtr_->gradRef().stkFieldPtr();
    const auto* oldStrainField =
        hasOld ? SijPtr_->prevTimeRef().stkFieldPtr() : nullptr;
    const auto* olderStrainField =
        useBDF2 ? SijPtr_->prevTimeRef().prevTimeRef().stkFieldPtr() : nullptr;
    auto* correctionField = frRef().stkFieldPtr();

    // Keep the field-data pointers and fixed-size work arrays outside the
    // bucket/node loops. Rebind field pointers for each STK bucket.
    const scalar* U = nullptr;
    const scalar* x = nullptr;
    const scalar* scales = nullptr;
    const scalar* kVals = nullptr;
    const scalar* gradStrain = nullptr;
    const scalar* strainOld = nullptr;
    const scalar* strainOlder = nullptr;
    const scalar* gradSflat = nullptr;
    scalar* fr = nullptr;
    scalar vel[3] = {};
    scalar dx[3] = {};
    scalar S[3][3] = {};
    scalar W[3][3] = {};
    scalar DS[3][3] = {};

    for (const auto* bucket : buckets)
    {
        U = stk::mesh::field_data(*UField, *bucket);
        x = stk::mesh::field_data(*coordsField, *bucket);
        scales = stk::mesh::field_data(*scaleField, *bucket);
        kVals = kEpsilon ? stk::mesh::field_data(*kField, *bucket) : nullptr;
        gradU = stk::mesh::field_data(*gradUField, *bucket);
        strain = stk::mesh::field_data(*strainField, *bucket);
        gradStrain = stk::mesh::field_data(*gradStrainField, *bucket);
        strainOld =
            hasOld ? stk::mesh::field_data(*oldStrainField, *bucket) : nullptr;
        strainOlder = useBDF2
                          ? stk::mesh::field_data(*olderStrainField, *bucket)
                          : nullptr;
        fr = stk::mesh::field_data(*correctionField, *bucket);

        for (stk::mesh::Bucket::size_type node = 0; node < bucket->size();
             ++node)
        {
            g = gradU + SPATIAL_DIM * SPATIAL_DIM * node;
            Sflat = strain + SPATIAL_DIM * SPATIAL_DIM * node;
            gradSflat =
                gradStrain + SPATIAL_DIM * SPATIAL_DIM * SPATIAL_DIM * node;

            for (label m = 0; m < SPATIAL_DIM; ++m)
            {
                vel[m] = U[SPATIAL_DIM * node + m];
                dx[m] = x[SPATIAL_DIM * node + m] - origin[m];
            }
            if (frameRotating)
            {
                for (label i = 0; i < 3; ++i)
                {
                    for (label j = 0; j < 3; ++j)
                    {
                        for (label k = 0; k < 3; ++k)
                        {
                            vel[i] -= eps[i][j][k] * OmRot[j] * dx[k];
                        }
                    }
                }
            }

            scalar S2 = 0.0;
            scalar W2 = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    S[i][j] = Sflat[SPATIAL_DIM * i + j];
                    W[i][j] =
                        0.5 * (g[SPATIAL_DIM * i + j] - g[SPATIAL_DIM * j + i]);
                    S2 += S[i][j] * S[i][j];
                    W2 += W[i][j] * W[i][j];
                }
            }
            S2 *= 2.0;
            W2 *= 2.0;
            const scalar SMag = std::sqrt(S2);
            const scalar WMag = std::sqrt(W2);
            const scalar rStar = SMag / (WMag + SMALL);

            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    scalar value = 0.0;
                    if (i < SPATIAL_DIM && j < SPATIAL_DIM)
                    {
                        for (label k = 0; k < SPATIAL_DIM; ++k)
                        {
                            value +=
                                vel[k] *
                                gradSflat[SPATIAL_DIM * (SPATIAL_DIM * i + j) +
                                          k];
                        }

                        // Full material derivative. Previous strain values
                        // are physical-time STK states, not outer iterations.
                        if (hasOld)
                        {
                            const label idx = SPATIAL_DIM * i + j;
                            const scalar current = Sflat[idx];
                            const scalar old =
                                strainOld[SPATIAL_DIM * SPATIAL_DIM * node +
                                          idx];
                            const scalar older =
                                useBDF2 ? strainOlder[SPATIAL_DIM *
                                                          SPATIAL_DIM * node +
                                                      idx]
                                        : 0.0;
                            value += dtInv * (timeCoeff[0] * current +
                                              timeCoeff[1] * old +
                                              timeCoeff[2] * older);
                        }
                    }
                    if (frameRotating)
                    {
                        for (label m = 0; m < 3; ++m)
                        {
                            for (label n = 0; n < 3; ++n)
                            {
                                value += OmRot[m] * (eps[i][m][n] * S[j][n] +
                                                     eps[j][m][n] * S[i][n]);
                            }
                        }
                    }
                    DS[i][j] = value;
                }
            }

            scalar numerator = 0.0;
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    scalar WS = 0.0;
                    for (label k = 0; k < 3; ++k)
                    {
                        WS += W[i][k] * S[j][k];
                    }
                    numerator += WS * DS[i][j];
                }
            }
            numerator *= 2.0;

            // For k-epsilon, obtain an equivalent omega = epsilon/(Cmu*k).
            const scalar omega =
                kEpsilon ? scales[node] * invCmu / std::max(kVals[node], SMALL)
                         : scales[node];
            const scalar D2 = std::max(S2, beta * omega * omega);
            const scalar D = std::sqrt(D2);
            const scalar rTilde = numerator / (WMag * D2 * D + SMALL);

            const scalar fRot = (1.0 + cR1_) * 2.0 * rStar / (1.0 + rStar) *
                                    (1.0 - cR3_ * std::atan(cR2_ * rTilde)) -
                                cR1_;
            const scalar limited = std::max(0.0, std::min(fRot, 1.25));
            fr[node] = std::max(0.0, 1.0 + Cscale * (limited - 1.0));
        }
    }

    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {correctionField});
    }
}

} /* namespace accel */
