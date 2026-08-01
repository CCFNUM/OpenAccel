// File       : SmagorinskyLESModel.cpp
// Created    : Fri Jul 31 2026
// Author     : GitHub Copilot
// Description: LES Smagorinsky turbulence model
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "SmagorinskyLESModel.h"
#include "mesh.h"
#include "messager.h"

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
    const scalar Cs = this->Cs();

    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    STKScalarField* mutSTKFieldPtr = this->mutRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
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
            const scalar mut_val = rho[k] * (Cs * delta) * (Cs * delta) * sijMag;

            mut[k] = urf * mut_val + (1.0 - urf) * mut[k];
        }
    }
}

} /* namespace accel */
