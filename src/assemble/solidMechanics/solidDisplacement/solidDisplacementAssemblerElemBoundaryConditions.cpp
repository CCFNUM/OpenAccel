// File       : solidDisplacementAssemblerElemBoundaryConditions.cpp
// Description: Boundary terms for the solid displacement equation
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"
#include "types.h"

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

    // Follower-pressure load-stiffness tangent (Bonet & Wood, Nonlinear
    // Continuum Mechanics for Finite Element Analysis, 2nd ed. 2008, §8.5,
    // eq. 8.65): a physical (constant Cauchy) pressure BC makes the external
    // force a function of the current displacement, so a fully consistent
    // Newton tangent needs an extra "pressure stiffness" contribution here,
    // in addition to the interior material/geometric stiffness. See the
    // guarded TODO below -- not yet implemented (Picard for now).
    auto& bc = model_->DRef().boundaryConditionRef(domain->index(),
                                                    boundary->index());
    auto& followerPressureData = bc.template data<1>("follower_pressure");
    const bool followerPressure =
        (followerPressureData.type() == inputDataType::constant) &&
        (*followerPressureData.value() > 0.5);

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

            if (followerPressure)
            {
                // TODO: pressure load-stiffness tangent (Bonet & Wood 2008
                // §8.5) -- Picard for now.
                //
                // `lhs` is intentionally left at zero (as in the dead-load
                // case above): the RHS just computed already carries the
                // correct follower-load force, because
                // solidMechanicsModel::updateDisplacementBoundarySideFieldTr
                // action_ recomputes the deformed exposed-area vector from
                // `coordinates + displacement` on every call, i.e. every
                // nonlinear iteration -- so the *converged* solution is
                // physically correct. What is missing here is only the
                // consistent Newton tangent contribution: since the applied
                // load itself depends on displacement (through the current
                // normal and current area), the exact Jacobian of this
                // boundary's residual acquires an extra "pressure/load
                // stiffness" term, on top of the interior material and
                // geometric stiffness already assembled elsewhere. Omitting
                // it degrades this boundary's contribution to Newton's
                // convergence rate to linear/Picard instead of quadratic;
                // it does not change the converged answer.
                //
                // Filling this in requires, per face/integration point: the
                // current unit normal n and current unit tangent t (both
                // derivable from the same deformed coordinates gathered in
                // updateDisplacementBoundarySideFieldTraction_, e.g. via a
                // second call to mesh::computeDeformedExposedAreaVector()
                // from here, or by having that data passed through), plus
                // the edge shape-function parametric derivative dN/ds; then
                // K_ab^ij = p * (N^a dN^b/ds t_i n_j - N^a dN^b/ds n_i t_j)
                // written into lhs[dofsPerSide*dofsPerSide] before the call
                // below, guarded by this same `followerPressure` flag so the
                // dead-load path (lhs left at zero) is unaffected.
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
