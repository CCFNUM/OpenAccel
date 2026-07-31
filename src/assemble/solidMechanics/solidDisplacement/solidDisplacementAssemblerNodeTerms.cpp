// File       : solidDisplacementAssemblerNodeTerms.cpp
// Created    : Fri Jan 31 2026 10:00:00 (+0100)
// Author     : Mhamad Mahdi Alloush
// Modified   : Added mass-proportional Rayleigh damping (alpha * rho * d(D)/dt)
//              matching solids4foam's "dampingCoeff" convention.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"
#include "solidMechanicsModel.h"

namespace accel
{

void solidDisplacementAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    // No mass terms for steady state
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

void solidDisplacementAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    // First order scheme not typically used for structural dynamics
    // (requires 3 time levels for second derivative)
    // Fall through to second order implementation
    assembleNodeTermsFusedSecondOrderUnsteady_(domain, ctx);
#else
    // FEM: First order scheme not typically used for structural dynamics
    // (requires 3 time levels for second derivative)
    // Fall through to second order implementation for now
    assembleNodeTermsFusedSecondOrderUnsteady_(domain, ctx);
#endif
}

void solidDisplacementAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    // ========================================================================
    // Second-order time discretization for structural dynamics:
    //
    //   rho * d^2 D / dt^2 ~ rho / dt^2 * (D^{k+1} - 2*D^n + D^{n-1})
    //
    // where D^{k+1} is the new outer-iteration displacement, D^n is the
    // converged displacement at the end of the previous time step, and
    // D^{n-1} is the converged displacement two time steps back.
    //
    // ------------------------------------------------------------------------
    // Mass-proportional Rayleigh damping (alpha)
    // ------------------------------------------------------------------------
    // When damping_coeff = alpha > 0, the additional term
    //
    //   alpha * rho * d(D)/dt ~ alpha * rho / dt * (D^{k+1} - D^n)
    //
    // is added to the equation. The discretization on the first time
    // derivative is BDF1 (Euler implicit), which matches the default
    // behaviour of solids4foam's fvm::ddt(D) in its solid models. The term
    // contributes
    //   LHS diagonal: + alpha * rho * V_0 / dt
    //   RHS         : - alpha * rho * V_0 / dt * (D^{k+1} - D^n)
    // in the residual form used elsewhere in this assembly; at convergence
    // (D^{k+1} -> D^n+1) it correctly enforces the damped momentum balance.
    // Reduces the under-damped oscillation seen in undamped FSI runs and
    // matches the behaviour of s4f's "dampingCoeff" entry in solidProperties.
    // ========================================================================

    auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // Space for LHS/RHS (single node)
    const label lhsSize = SPATIAL_DIM * SPATIAL_DIM;
    const label rhsSize = SPATIAL_DIM;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // Pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get displacement fields (current iter, prev time, prev-prev time)
    const STKScalarField* DSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* DSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* DSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();

    // Rayleigh mass-proportional damping coefficient (1/s). When > 0, the
    // term alpha * rho * d(D)/dt is added to the equation (BDF1 stencil).
    const scalar dampingCoeff = domain->solidMechanics_.dampingCoeff_;
    const bool useDamping = (dampingCoeff > 0.0);

    // Get density field (constant for solid mechanics)
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();

    // Geometric fields - use original volume for total Lagrangian
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Time step
    const scalar dt = mesh.controlsRef().getTimestep();
    const scalar rDeltaT2 = 1.0 / (dt * dt);
    const scalar rDeltaT = 1.0 / dt; // needed for damping (1st time deriv)

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        // Field chunks in bucket
        scalar* Db = stk::mesh::field_data(*DSTKFieldPtr, nodeBucket);
        scalar* DbOld = stk::mesh::field_data(*DSTKFieldPtrOld, nodeBucket);
        scalar* DbOldOld =
            stk::mesh::field_data(*DSTKFieldPtrOldOld, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // Get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            // Zero lhs/rhs
            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // Get nodal values
            const scalar rho = rhob[iNode];
            const scalar vol = volb[iNode];

            // LHS coefficient: rho * V / dt^2
            const scalar lhsfac = rho * vol * rDeltaT2;

            // Damping LHS factor (alpha * rho * V_0 / dt). Computed once per
            // node and zeroed when damping is disabled, so the inner-loop
            // branch is predictable and cheap.
            const scalar lhsfacDamp =
                useDamping ? (dampingCoeff * rho * vol * rDeltaT) : 0.0;

            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const scalar Di = Db[SPATIAL_DIM * iNode + i];
                const scalar DiOld = DbOld[SPATIAL_DIM * iNode + i];
                const scalar DiOldOld = DbOldOld[SPATIAL_DIM * iNode + i];

                // LHS: add to diagonal
                p_lhs[i * SPATIAL_DIM + i] += lhsfac;

                // RHS: rho*V/dt^2 * (D - 2*D_old + D_old_old)
                p_rhs[i] -= lhsfac * (Di - 2.0 * DiOld + DiOldOld);

                // Mass-proportional Rayleigh damping
                //   LHS diag : + alpha * rho * V / dt
                //   RHS      : - alpha * rho * V / dt * (D - D_old)
                // (residual form; at convergence D -> D^{n+1} this enforces
                //  alpha*rho*(D^{n+1} - D^n)/dt, matching s4f's
                //  dampingCoeff*rho*fvm::ddt(D))
                if (useDamping)
                {
                    p_lhs[i * SPATIAL_DIM + i] += lhsfacDamp;
                    p_rhs[i] -= lhsfacDamp * (Di - DiOld);
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} // namespace accel
