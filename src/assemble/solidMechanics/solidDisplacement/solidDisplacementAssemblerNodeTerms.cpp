// File       : solidDisplacementAssemblerNodeTerms.cpp
// Created    : Fri Jan 31 2026 10:00:00 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"
#include "solidMechanicsModel.h"

namespace accel
{

void solidDisplacementAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
    // No mass terms for steady state
}

void solidDisplacementAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    // solid_transient_scheme selects the d2dt2 operator (shared dispatcher)
    assembleNodeTermsFusedSecondOrderUnsteady_(domain, ctx);
}

void solidDisplacementAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
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

    // Get displacement field
    const STKScalarField* DSTKFieldPtr = phi_->stkFieldPtr();

    // Rayleigh mass-proportional damping coefficient [1/s]
    const scalar dampingCoeff = domain->solidMechanics_.dampingCoeff_;
    const bool useDamping = (dampingCoeff > 0.0);

    // Get density field (constant for solid mechanics)
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();

    // Geometric fields - use original volume for total Lagrangian
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Time step
    const scalar dt = mesh.controlsRef().getTimestep();

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

    // solid_transient_scheme selects the d2dt2 operator
    if (field_broker_->controlsRef().getSolidTransientScheme() ==
        solidTransientSchemeType::backwardD2dt2)
    {
        // backward d2dt2: (2.25D - 6D1 + 5.5D2 - 2D3 + 0.25D4)/dt^2, 5 states
        const STKScalarField* DSTKFieldPtr1 = phi_->prevTimeRef().stkFieldPtr();
        const STKScalarField* DSTKFieldPtr2 =
            phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
        const STKScalarField* DSTKFieldPtr3 =
            phi_->prevTimeRef().prevTimeRef().prevTimeRef().stkFieldPtr();
        const STKScalarField* DSTKFieldPtr4 = phi_->prevTimeRef()
                                                  .prevTimeRef()
                                                  .prevTimeRef()
                                                  .prevTimeRef()
                                                  .stkFieldPtr();

        const scalar dt0 = mesh.controlsRef().getTimestep(-1);

        // BDF2 outer coefficients (S4F backwardD2dt2Scheme)
        const scalar safedt0 = (dt0 > 0.0) ? dt0 : dt;
        const scalar coefft = 1.0 + dt / (dt + safedt0);
        const scalar coefft00 = dt * dt / (safedt0 * (dt + safedt0));
        const scalar coefft0 = coefft + coefft00;

        // inner BDF2 (constant dt): (1.5 f^n - 2 f^{n-1} + 0.5 f^{n-2})/dt
        const scalar inner_c0 = 1.5;
        const scalar inner_c1 = -2.0;
        const scalar inner_c2 = 0.5;

        const scalar rDeltaT = 1.0 / dt;
        const scalar rDeltaT2 = rDeltaT * rDeltaT;

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            const scalar* Db = stk::mesh::field_data(*DSTKFieldPtr, nodeBucket);
            const scalar* Db1 =
                stk::mesh::field_data(*DSTKFieldPtr1, nodeBucket);
            const scalar* Db2 =
                stk::mesh::field_data(*DSTKFieldPtr2, nodeBucket);
            const scalar* Db3 =
                stk::mesh::field_data(*DSTKFieldPtr3, nodeBucket);
            const scalar* Db4 =
                stk::mesh::field_data(*DSTKFieldPtr4, nodeBucket);
            const scalar* rhob =
                stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            const scalar* volb =
                stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                stk::mesh::Entity node = nodeBucket[iNode];
                connectedNodes[0] = node;

                const scalar rho = rhob[iNode];
                const scalar vol = volb[iNode];
                const scalar fac = rho * vol * rDeltaT2;

                // BDF2-consistent Rayleigh damping alpha*rho*dD/dt
                const scalar lhsfacDamp =
                    useDamping ? (dampingCoeff * rho * vol * inner_c0 * rDeltaT)
                               : 0.0;

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    lhs[i * SPATIAL_DIM + i] = fac * coefft * inner_c0;

                    // RHS: nested BDF2 expansion of the d2dt2 stencil
                    const scalar Di = Db[SPATIAL_DIM * iNode + i];
                    const scalar Di1 = Db1[SPATIAL_DIM * iNode + i];
                    const scalar Di2 = Db2[SPATIAL_DIM * iNode + i];
                    const scalar Di3 = Db3[SPATIAL_DIM * iNode + i];
                    const scalar Di4 = Db4[SPATIAL_DIM * iNode + i];

                    // BDF2 at each level:
                    const scalar bdf2_n =
                        inner_c0 * Di + inner_c1 * Di1 + inner_c2 * Di2;
                    const scalar bdf2_n1 =
                        inner_c0 * Di1 + inner_c1 * Di2 + inner_c2 * Di3;
                    const scalar bdf2_n2 =
                        inner_c0 * Di2 + inner_c1 * Di3 + inner_c2 * Di4;

                    rhs[i] = -fac * (coefft * bdf2_n - coefft0 * bdf2_n1 +
                                     coefft00 * bdf2_n2);

                    if (useDamping)
                    {
                        lhs[i * SPATIAL_DIM + i] += lhsfacDamp;
                        rhs[i] -=
                            dampingCoeff * rho * vol * rDeltaT *
                            (inner_c0 * Di + inner_c1 * Di1 + inner_c2 * Di2);
                    }
                }

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
        // default: Euler d2dt2 rho*V0/dt^2*(D - 2D_old + D_old2), as S4F
        // Pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];

        const STKScalarField* DSTKFieldPtrOld =
            phi_->prevTimeRef().stkFieldPtr();
        const STKScalarField* DSTKFieldPtrOldOld =
            phi_->prevTimeRef().prevTimeRef().stkFieldPtr();

        const scalar rDeltaT2 = 1.0 / (dt * dt);
        const scalar rDeltaT = 1.0 / dt; // first time derivative of the damping

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // Field chunks in bucket
            scalar* Db = stk::mesh::field_data(*DSTKFieldPtr, nodeBucket);
            scalar* DbOld = stk::mesh::field_data(*DSTKFieldPtrOld, nodeBucket);
            scalar* DbOldOld =
                stk::mesh::field_data(*DSTKFieldPtrOldOld, nodeBucket);
            scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
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

                // LHS coefficient: rho*V/dt^2
                const scalar lhsfac = rho * vol * rDeltaT2;

                // damping LHS coefficient: alpha*rho*V/dt
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

                    // Rayleigh damping alpha*rho*V/dt (BDF1)
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
    }
}

} // namespace accel
