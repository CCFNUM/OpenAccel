// File       : solidDisplacementAssemblerElemBoundaryConditions.cpp
// Description: Boundary terms for the solid displacement equation
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"

namespace accel
{

#ifndef USE_CVFEM_SOLID_MECHANICS
void solidDisplacementAssembler::assembleElemTermsBoundaryWallSpecifiedFlux_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto& sideTractionField = phi_->sideFluxFieldRef().stkFieldRef();
    const auto& exposedAreaVectorField = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    const stk::mesh::Selector selectedSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(boundary->parts());
    const auto& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selectedSides);

    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;
    std::vector<scalar> shapeFunctions;
    std::array<scalar, SPATIAL_DIM> integratedForce{};
    label assembledFaces = 0;

    for (const stk::mesh::Bucket* bucket : sideBuckets)
    {
        MasterElement* faceMasterElement =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label nodesPerSide = faceMasterElement->nodesPerElement_;
        const label integrationPoints = faceMasterElement->numIntPoints_;
        const label dofsPerSide = nodesPerSide * SPATIAL_DIM;

        lhs.assign(dofsPerSide * dofsPerSide, 0.0);
        rhs.resize(dofsPerSide);
        scratchIds.resize(dofsPerSide);
        scratchVals.resize(dofsPerSide);
        connectedNodes.resize(nodesPerSide);
        shapeFunctions.resize(integrationPoints * nodesPerSide);

        if (phi_->isShifted())
            faceMasterElement->shifted_shape_fcn(shapeFunctions.data());
        else
            faceMasterElement->shape_fcn(shapeFunctions.data());

        for (stk::mesh::Entity side : *bucket)
        {
            ++assembledFaces;
            std::fill(rhs.begin(), rhs.end(), 0.0);

            const stk::mesh::Entity* sideNodes = bulkData.begin_nodes(side);
            STK_ThrowAssert(bulkData.num_nodes(side) == nodesPerSide);
            for (label node = 0; node < nodesPerSide; ++node)
                connectedNodes[node] = sideNodes[node];

            const scalar* traction =
                stk::mesh::field_data(sideTractionField, side);
            const scalar* areaVector =
                stk::mesh::field_data(exposedAreaVectorField, side);

            for (label ip = 0; ip < integrationPoints; ++ip)
            {
                scalar areaMagnitudeSquared = 0.0;
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    const scalar area =
                        areaVector[ip * SPATIAL_DIM + dim];
                    areaMagnitudeSquared += area * area;
                }
                const scalar areaMagnitude =
                    std::sqrt(areaMagnitudeSquared);

                for (label node = 0; node < nodesPerSide; ++node)
                {
                    const scalar weight =
                        shapeFunctions[ip * nodesPerSide + node] *
                        areaMagnitude;
                    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                    {
                        const scalar force =
                            weight * traction[ip * SPATIAL_DIM + dim];
                        rhs[node * SPATIAL_DIM + dim] += force;
                        integratedForce[dim] += force;
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }

    const char* femDebug = std::getenv("OPENACCEL_FEM_DEBUG");
    if (femDebug && std::atoi(femDebug) > 0)
    {
        std::array<scalar, SPATIAL_DIM> globalIntegratedForce{};
        label globalAssembledFaces = 0;
        stk::all_reduce_sum(bulkData.parallel(),
                            integratedForce.data(),
                            globalIntegratedForce.data(),
                            SPATIAL_DIM);
        stk::all_reduce_sum(bulkData.parallel(),
                            &assembledFaces,
                            &globalAssembledFaces,
                            1);

        if (messager::myProcNo() == 0)
        {
            std::cout << "[FEM DEBUG] Boundary `" << boundary->name()
                      << "`: faces=" << globalAssembledFaces << ", force=";
            for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                std::cout << (dim == 0 ? " [" : ", ")
                          << globalIntegratedForce[dim];
            std::cout << "]\n";
        }
    }
}
#endif

} /* namespace accel */
