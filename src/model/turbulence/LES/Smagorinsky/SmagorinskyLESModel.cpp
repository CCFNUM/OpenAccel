// File       : SmagorinskyLESModel.cpp
// Created    : Fri Jul 31 2026
// Author     : GitHub Copilot
// Description: LES Smagorinsky turbulence model
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "SmagorinskyLESModel.h"
#include "mesh.h"
#include "messager.h"

#include <algorithm>
#include <cmath>

namespace accel
{

SmagorinskyLESModel::SmagorinskyLESModel(realm* realm)
    : turbulenceModel(realm)
{
}

void SmagorinskyLESModel::updateTurbulentDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    const scalar urf = this->mutRef().urf();

    // model coefficients are owned per domain (see `struct turbulence`)
    const scalar Cs = domain->turbulence_.smagorinskyConstant_;
    const scalar APlus = domain->turbulence_.vanDriestAPlus_;

    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    STKScalarField* mutSTKFieldPtr = this->mutRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = this->muRef().stkFieldPtr();
    const STKScalarField* dualNodalVolumeSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                   mesh::dual_nodal_volume_ID);
    const STKScalarField& gradUSTKFieldRef =
        this->URef().gradRef().stkFieldRef();

    if (!dualNodalVolumeSTKFieldPtr)
    {
        errorMsg("SmagorinskyLESModel: dual_nodal_volume field is not available");
        return;
    }

    // Van Driest near-wall damping. It is only meaningful when the domain
    // actually bounds a no-slip wall: without one, `minimum_distance_to_wall`
    // carries no wall information (it is zero-filled) and applying the damping
    // would erroneously annihilate the SGS viscosity everywhere. In that case
    // the model falls back to undamped Smagorinsky.
    const STKScalarField* yMinSTKFieldPtr =
        this->yMinRef().stkFieldPtr();

    const bool applyWallDamping =
        !this->collectNoSlipWallParts_(domain).empty() &&
        yMinSTKFieldPtr != nullptr && muSTKFieldPtr != nullptr;

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        const scalar* rho = stk::mesh::field_data(*rhoSTKFieldPtr, b);
        const scalar* dualVol =
            stk::mesh::field_data(*dualNodalVolumeSTKFieldPtr, b);
        scalar* mut = stk::mesh::field_data(*mutSTKFieldPtr, b);

        const scalar* mu =
            applyWallDamping ? stk::mesh::field_data(*muSTKFieldPtr, b)
                             : nullptr;
        const scalar* yMin =
            applyWallDamping ? stk::mesh::field_data(*yMinSTKFieldPtr, b)
                             : nullptr;

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            const scalar* dudx = stk::mesh::field_data(gradUSTKFieldRef, b[k]);

            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 * (dudx[SPATIAL_DIM * i + j] +
                               dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            const scalar delta = std::cbrt(dualVol[k]);

            // Van Driest damping of the SGS mixing length:
            //   nu          = mu / rho
            //   uTauLocal   = sqrt(nu * |S|)
            //   yPlus       = y * uTauLocal / nu = y * sqrt(|S| / nu)
            //   fD          = 1 - exp(-yPlus / A+)
            //   l           = Cs * delta * fD
            //
            // uTauLocal is exact at the wall once the damping is effective,
            // since the near-wall stress is then carried entirely by the
            // molecular viscosity (tau_w = mu * |S|). The residual coupling is
            // self-correcting: over-predicting yPlus raises mut, which lowers
            // the resolved |S|, which lowers yPlus again.
            //
            // Degenerate inputs (non-positive rho, mu or |S|) carry no wall
            // information, so damping is skipped and the undamped Smagorinsky
            // value is retained rather than silently zeroing the model.
            scalar fD = 1.0;
            if (applyWallDamping && rho[k] > SMALL && mu[k] > SMALL &&
                sijMag > SMALL)
            {
                const scalar nu = mu[k] / rho[k];
                const scalar y = std::max(yMin[k], scalar(0));
                const scalar yPlus = y * std::sqrt(sijMag / nu);

                fD = 1.0 - std::exp(-yPlus / APlus);
            }

            const scalar lMix = Cs * delta * fD;
            const scalar mut_val = rho[k] * lMix * lMix * sijMag;

            mut[k] = urf * mut_val + (1.0 - urf) * mut[k];
        }
    }
}

} /* namespace accel */
