// File : navierStokesAssemblerElemBoundaryConditions.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "flowModel.h"
#include "navierStokesAssembler.h"
#include "thermoModel.h"

namespace accel
{

void navierStokesAssembler::assembleElemTermsBoundary_(const domain* domain,
                                                       Context* ctx)
{
    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        const boundaryPhysicalType type = boundary->type();

        const boundaryConditionType UBCType =
            model_->URef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();
        const boundaryConditionType pBCType =
            model_->pRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                {
                    assembleElemTermsBoundarySymmetry_(domain, boundary, ctx);
                }
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::noSlip:
                            {
                                assembleElemTermsBoundaryWallNoSlip_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::slip:
                            {
                                assembleElemTermsBoundarySlipWall_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg("invalid velocity boundary "
                                     "condition at wall");
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::staticPressure:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    case boundaryConditionType::zeroGradient:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedVelocity_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        case boundaryConditionType::normalSpeed:
                        case boundaryConditionType::massFlowRate:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::zeroGradient:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedVelocity_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        case boundaryConditionType::specifiedDirection:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::staticPressure:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    case boundaryConditionType::totalPressure:
                                        {
                                            assembleElemTermsBoundaryInletSpecifiedTotalPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at inlet");
                                }
                            }
                            break;

                        default:
                            errorMsg("invalid velocity boundary "
                                     "condition at inlet");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::staticPressure:
                                    case boundaryConditionType::
                                        averageStaticPressure:
                                        {
                                            assembleElemTermsBoundaryOutletSpecifiedPressure_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    case boundaryConditionType::zeroGradient:
                                        {
                                            assembleElemTermsBoundaryOutletOutflow_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at outlet");
                                }
                                break;
                            }
                            break;

                        case boundaryConditionType::massFlowRate:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::massFlowRate:
                                        {
                                            assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at outlet");
                                }
                                break;
                            }
                            break;

                        default:
                            errorMsg("invalid velocity boundary "
                                     "condition at outlet");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::specifiedDirection:
                            {
                                switch (pBCType)
                                {
                                    case boundaryConditionType::staticPressure:
                                    case boundaryConditionType::totalPressure:
                                        {
                                            assembleElemTermsBoundaryOpening_(
                                                domain, boundary, ctx);
                                        }
                                        break;

                                    default:
                                        errorMsg("invalid pressure boundary "
                                                 "condition at opening");
                                }
                                break;
                            }

                        default:
                            errorMsg("invalid velocity boundary "
                                     "condition at opening");
                    }
                    break;
                }
        }
    }
}

// FIXME: [2024-10-28] Equation terms assembly kernels.
// Macros are being used here to avoid error prone code duplication and to make
// the code where these macros are used more readable. This is not an ideal
// approach, assembly kernels should be generated with a source code generator
// approach using some high-level DSL (also helps with porting to different
// architectures).
#define IP_EXPLICIT_ADVECTIVE_FLUX__(flux)                                     \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * explicit advective flux                                             \
         * */                                                                  \
        for (label i = 0; i < SPATIAL_DIM; ++i)                                \
        {                                                                      \
            const label indexR = nearestNode * SPATIAL_DIM + i;                \
            p_rhs[indexR] -= (flux);                                           \
        }                                                                      \
    } while (0)

#define IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(flux)                              \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * implicit advective flux                                             \
         * */                                                                  \
        for (label i = 0; i < SPATIAL_DIM; ++i)                                \
        {                                                                      \
            const label indexR = nearestNode * SPATIAL_DIM + i;                \
            const label rowR = indexR * nodesPerElement * SPATIAL_DIM;         \
            p_rhs[indexR] -= (flux);                                           \
            p_lhs[rowR + nearestNode * SPATIAL_DIM + i] +=                     \
                tmDot; /* upwind lhs */                                        \
        }                                                                      \
    } while (0)

#define IP_IMPLICIT_ADVECTIVE_FLUX__(flux)                                     \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * implicit advective flux                                             \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerSide; ++ic)                            \
        {                                                                      \
            const label inn = faceNodeOrdinals[ic];                            \
            const scalar r = p_velocity_face_shape_function[offSetSF + ic];    \
            for (label i = 0; i < SPATIAL_DIM; ++i)                            \
            {                                                                  \
                const label indexR = nearestNode * SPATIAL_DIM + i;            \
                const label rowR = indexR * nodesPerElement * SPATIAL_DIM;     \
                p_rhs[indexR] -= (flux);                                       \
                p_lhs[rowR + inn * SPATIAL_DIM + i] += tmDot * r;              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define IP_FULL_STRESS__()                                                     \
    do                                                                         \
    {                                                                          \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label off =                                                  \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
                                                                               \
            /* divU = sum_k dN/dx_k * u_k (at this ip, for node ic) */         \
            scalar divU_ic = 0.0;                                              \
            for (label k = 0; k < SPATIAL_DIM; ++k)                            \
            {                                                                  \
                const scalar dndxk = p_dndx[off + k];                          \
                const scalar uk = p_U[ic * SPATIAL_DIM + k];                   \
                divU_ic += dndxk * uk;                                         \
            }                                                                  \
                                                                               \
            for (label i = 0; i < SPATIAL_DIM; ++i)                            \
            {                                                                  \
                const scalar Ai = areaVec[faceOffSet + i];                     \
                const label indexR = nearestNode * SPATIAL_DIM + i;            \
                const label rowR = indexR * nodesPerElement * SPATIAL_DIM;     \
                const scalar uxi = p_U[ic * SPATIAL_DIM + i];                  \
                const scalar dndxi = p_dndx[off + i];                          \
                                                                               \
                /* shear parts: -mu*dui/dxj*Aj and -mu*duj/dxi*Aj */           \
                for (label j = 0; j < SPATIAL_DIM; ++j)                        \
                {                                                              \
                    const scalar Aj = areaVec[faceOffSet + j];                 \
                    const scalar dndxj = p_dndx[off + j];                      \
                    const scalar uxj = p_U[ic * SPATIAL_DIM + j];              \
                                                                               \
                    /* -mu * dui/dxj * Aj */                                   \
                    scalar lhsfac = -muEffBip * dndxj * Aj;                    \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi;                             \
                                                                               \
                    /* -mu * duj/dxi * Aj (transpose) */                       \
                    lhsfac = -muEffBip * dndxi * Aj;                           \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                }                                                              \
                                                                               \
                /* bulk (Stokes) term: lambda*divU*Ai*mu */                    \
                const scalar divUStress =                                      \
                    (2.0 / 3.0) * muEffBip * divU_ic * Ai * comp;              \
                p_rhs[indexR] -= divUStress;                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

#define IP_FULL_STRESS_FIXED_VEL__()                                           \
    do                                                                         \
    {                                                                          \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label off =                                                  \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
                                                                               \
            /* divU = sum_k dN/dx_k * u_k (at this ip, for node ic) */         \
            scalar divU_ic = 0.0;                                              \
            for (label k = 0; k < SPATIAL_DIM; ++k)                            \
            {                                                                  \
                const scalar dndxk = p_dndx[off + k];                          \
                const scalar uk = p_U[ic * SPATIAL_DIM + k];                   \
                divU_ic += dndxk * uk;                                         \
            }                                                                  \
                                                                               \
            for (label i = 0; i < SPATIAL_DIM; ++i)                            \
            {                                                                  \
                const scalar Ai = areaVec[faceOffSet + i];                     \
                const label indexR = nearestNode * SPATIAL_DIM + i;            \
                const label rowR = indexR * nodesPerElement * SPATIAL_DIM;     \
                const scalar uxi = p_U[ic * SPATIAL_DIM + i];                  \
                const scalar dndxi = p_dndx[off + i];                          \
                                                                               \
                /* shear parts: -mu*dui/dxj*Aj and -mu*duj/dxi*Aj */           \
                for (label j = 0; j < SPATIAL_DIM; ++j)                        \
                {                                                              \
                    const scalar Aj = areaVec[faceOffSet + j];                 \
                    const scalar dndxj = p_dndx[off + j];                      \
                    const scalar uxj = p_U[ic * SPATIAL_DIM + j];              \
                                                                               \
                    /* -mu * dui/dxj * Aj */                                   \
                    scalar lhsfac = -muEffBip * dndxj * Aj;                    \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] +=                      \
                        lhsfac * p_bcMultiplier[ic];                           \
                    p_rhs[indexR] -= lhsfac * uxi;                             \
                                                                               \
                    /* -mu * duj/dxi * Aj (transpose) */                       \
                    lhsfac = -muEffBip * dndxi * Aj;                           \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] +=                      \
                        lhsfac * p_bcMultiplier[ic];                           \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                }                                                              \
                                                                               \
                /* bulk (Stokes) term: lambda*divU*Ai*mu */                    \
                const scalar divUStress =                                      \
                    (2.0 / 3.0) * muEffBip * divU_ic * Ai * comp;              \
                p_rhs[indexR] -= divUStress;                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

#define IP_ZERO_NORMAL_STRESS__()                                              \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * zero normal stress                                                  \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar dndxi = p_dndx[offSetDnDx + i];               \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar om_nxinxi = 1.0 - nxi * nxi;                  \
                                                                               \
                    /* -mu*dui/dxj*Aj(1.0-nini); sneak in divU (explicit) */   \
                    scalar lhsfac = -muEffBip * dndxj * axj * om_nxinxi;       \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * om_nxinxi;    \
                                                                               \
                    /* -mu*duj/dxi*Aj(1.0-nini) */                             \
                    lhsfac = -muEffBip * dndxi * axj * om_nxinxi;              \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                                                                               \
                    /* now we need the -nx*ny*Fy - nx*nz*Fz part */            \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                            const scalar dndxl = p_dndx[offSetDnDx + l];       \
                                                                               \
                            /* +ni*nl*mu*dul/dxj*Aj; sneak in divU (explicit): \
                             */                                                \
                            lhsfac = muEffBip * dndxj * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                                                                               \
                            /* +ni*nl*mu*duj/dxl*Aj */                         \
                            lhsfac = muEffBip * dndxl * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;      \
                            p_rhs[indexR] -= lhsfac * uxj;                     \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

// FIXME: [2024-11-24] The divUstress term below should use axi
// instead of axj because Kronecker delta will project only the diagonal part.
// In the second part where i != l the divU stress should be zero since
// Kronecker delta is zero for indices not equal to each other (NEEDS CHECKING).
// This issue exists in various places of the code, not just in the kernels
// below.

// WARNING: [2024-11-25] This implementation does not account for the transpose
// term of the velocity gradient. An alternative approach includes the transpose
// term in the assembly.
#define IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__()                           \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * zero normal stress (no grad(U).transpose() term)                    \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar om_nxinxi = 1.0 - nxi * nxi;                  \
                                                                               \
                    /* -mu*dui/dxj*Aj(1.0-nini); sneak in divU (explicit) */   \
                    scalar lhsfac = -muEffBip * dndxj * axj * om_nxinxi;       \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * om_nxinxi;    \
                                                                               \
                    /* -nx*ny*Fy - nx*nz*Fz part */                            \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                                                                               \
                            /* +ni*nl*mu*dul/dxj*Aj; sneak in divU (explicit)  \
                             */                                                \
                            lhsfac = muEffBip * dndxj * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

// Directional stress term for specified direction inlet boundary conditions.
// Applies a transformation to the full stress tensor (including transpose
// gradient terms) that couples velocity only through the normal component.
// This ensures flow enters in the prescribed direction without imposing
// artificial tangential stress constraints. Requires 'dir' array in scope.
#define IP_DIRECTIONAL_STRESS__()                                              \
    do                                                                         \
    {                                                                          \
        /* Compute direction dot normal */                                     \
        scalar dir_dot_n = 0.0;                                                \
        for (label m = 0; m < SPATIAL_DIM; ++m)                                \
        {                                                                      \
            dir_dot_n += dir[SPATIAL_DIM * ip + m] * p_nx[m];                  \
        }                                                                      \
        const scalar inv_dir_dot_n = 1.0 / (dir_dot_n + SMALL);                \
                                                                               \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
                                                                               \
            /* Compute trace-like term: Σⱼ Aⱼ ∂N/∂xⱼ */             \
            scalar trace_term = 0.0;                                           \
            for (label k = 0; k < SPATIAL_DIM; ++k)                            \
            {                                                                  \
                const scalar dndxk = p_dndx[offSetDnDx + k];                   \
                const scalar axk = areaVec[faceOffSet + k];                    \
                trace_term += axk * dndxk;                                     \
            }                                                                  \
                                                                               \
            /* Compute area-direction dot product: Σₖ Aₖ dₖ */          \
            scalar area_dot_dir = 0.0;                                         \
            for (label k = 0; k < SPATIAL_DIM; ++k)                            \
            {                                                                  \
                const scalar axk = areaVec[faceOffSet + k];                    \
                const scalar dk = dir[SPATIAL_DIM * ip + k];                   \
                area_dot_dir += axk * dk;                                      \
            }                                                                  \
                                                                               \
            /* For each momentum equation */                                   \
            for (label i = 0; i < SPATIAL_DIM; ++i)                            \
            {                                                                  \
                const scalar dndxi = p_dndx[offSetDnDx + i];                   \
                const scalar di = dir[SPATIAL_DIM * ip + i];                   \
                                                                               \
                /* Project stress coefficients onto direction:                 \
                 * includes both gradient and transpose gradient terms */      \
                const scalar stress_proj_dir =                                 \
                    muEffBip * (di * trace_term + dndxi * area_dot_dir);       \
                                                                               \
                /* Apply transformation to couple through normal velocity */   \
                const label indexR = nearestNode * SPATIAL_DIM + i;            \
                const label rowR = indexR * nodesPerElement * SPATIAL_DIM;     \
                                                                               \
                for (label j = 0; j < SPATIAL_DIM; ++j)                        \
                {                                                              \
                    const scalar nxj = p_nx[j];                                \
                    const scalar uxj = p_U[ic * SPATIAL_DIM + j];              \
                                                                               \
                    /* Transformed coefficient: only n·u contributes */       \
                    const scalar lhsfac =                                      \
                        stress_proj_dir * nxj * inv_dir_dot_n;                 \
                                                                               \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define IP_ZERO_NORMAL_STRESS_FIXED_VEL__()                                    \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * zero normal stress                                                  \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar dndxi = p_dndx[offSetDnDx + i];               \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar om_nxinxi = 1.0 - nxi * nxi;                  \
                                                                               \
                    /* -mu*dui/dxj*Aj(1.0-nini) */                             \
                    scalar lhsfac = -muEffBip * dndxj * axj * om_nxinxi;       \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] +=                      \
                        lhsfac * p_bcMultiplier[ic];                           \
                    p_rhs[indexR] -= lhsfac * uxi;                             \
                                                                               \
                    /* -mu*duj/dxi*Aj(1.0-nini) */                             \
                    lhsfac = -muEffBip * dndxi * axj * om_nxinxi;              \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] +=                      \
                        lhsfac * p_bcMultiplier[ic];                           \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                                                                               \
                    /* now we need the -nx*ny*Fy - nx*nz*Fz part */            \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                            const scalar dndxl = p_dndx[offSetDnDx + l];       \
                                                                               \
                            /* +ni*nl*mu*dul/dxj*Aj                            \
                             */                                                \
                            lhsfac = muEffBip * dndxj * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] +=              \
                                lhsfac * p_bcMultiplier[ic];                   \
                            p_rhs[indexR] -= lhsfac * uxl;                     \
                                                                               \
                            /* +ni*nl*mu*duj/dxl*Aj */                         \
                            lhsfac = muEffBip * dndxl * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + j] +=              \
                                lhsfac * p_bcMultiplier[ic];                   \
                            p_rhs[indexR] -= lhsfac * uxj;                     \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define IP_ZERO_TANGENTIAL_STRESS__()                                          \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * zero tangential stress                                              \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    /* matrix entries */                                       \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar dndxi = p_dndx[offSetDnDx + i];               \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar nxinxi = nxi * nxi;                           \
                                                                               \
                    /* -mu*dui/dxj*Aj*ni*ni; sneak in divU (explicit): */      \
                    scalar lhsfac = -muEffBip * dndxj * axj * nxinxi;          \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * nxinxi;       \
                                                                               \
                    /* -mu*duj/dxi*Aj*ni*ni */                                 \
                    lhsfac = -muEffBip * dndxi * axj * nxinxi;                 \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                                                                               \
                    /* now we need the +nx*ny*Fy + nx*nz*Fz part */            \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                            const scalar dndxl = p_dndx[offSetDnDx + l];       \
                                                                               \
                            /* -ni*nl*mu*dul/dxj*Aj; sneak in divU (explicit): \
                             */                                                \
                            lhsfac = -muEffBip * dndxj * axj * nxinxl;         \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                                                                               \
                            /* -ni*nl*mu*duj/dxl*Aj */                         \
                            lhsfac = -muEffBip * dndxl * axj * nxinxl;         \
                            p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;      \
                            p_rhs[indexR] -= lhsfac * uxj;                     \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

// WARNING: [2024-11-25] This implementation does account for _positive_
// contributions of the transpose term of the velocity gradient (clipped to
// > 0.0).
// FIXME: [2024-11-25] Current implementation differs from the
// legacy implementation
#define IP_ZERO_TANGENTIAL_STRESS_CLIPPED_GRADU_TRANSPOSE__()                  \
    do                                                                         \
    {                                                                          \
        /*                                                                     \
         * zero tangential stress (velocity gradient transpose terms clipped)  \
         * */                                                                  \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    /* matrix entries */                                       \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar dndxi = p_dndx[offSetDnDx + i];               \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar nxinxi = nxi * nxi;                           \
                                                                               \
                    /* -mu*dui/dxj*Aj*ni*ni; sneak in divU (explicit): */      \
                    scalar lhsfac = -muEffBip * dndxj * axj * nxinxi;          \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * nxinxi;       \
                                                                               \
                    /* -mu*duj/dxi*Aj*ni*ni (clipped to > 0.0) */              \
                    lhsfac = -muEffBip * dndxi * axj * nxinxi;                 \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] +=                      \
                        (lhsfac > 0.0) ? lhsfac : 0.0;                         \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                                                                               \
                    /* +nx*ny*Fy + nx*nz*Fz part */                            \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                            const scalar dndxl = p_dndx[offSetDnDx + l];       \
                                                                               \
                            /* -ni*nl*mu*dul/dxj*Aj; sneak in divU (explicit): \
                             */                                                \
                            lhsfac = -muEffBip * dndxj * axj * nxinxl;         \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                                                                               \
                            /* -ni*nl*mu*duj/dxl*Aj (clipped to > 0.0) */      \
                            lhsfac = -muEffBip * dndxl * axj * nxinxl;         \
                            p_lhs[rowR + ic * SPATIAL_DIM + j] +=              \
                                (lhsfac > 0.0) ? lhsfac : 0.0;                 \
                            p_rhs[indexR] -= lhsfac * uxj;                     \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

void navierStokesAssembler::assembleElemTermsBoundaryInletSpecifiedVelocity_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& nodalSideUSTKFieldRef =
        model_->URef().nodeSideFieldRef().stkFieldRef();
    const auto& sideUSTKFieldRef = model_->URef().sideFieldRef().stkFieldRef();

    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_bcMultiplier[ni] = 1.0;

                // gather vectors
                const label offSet = ni * SPATIAL_DIM;

                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coords[offSet + j] = coords[j];
                    p_U[offSet + j] = U[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather vectors
                scalar* U = stk::mesh::field_data(nodalSideUSTKFieldRef, node);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    // overwrite velocity value at boundary with the dirichlet
                    // value
                    p_U[ic * SPATIAL_DIM + i] = U[i];
                }
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // calc vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                //================================
                // advection: entering domain
                //================================
                const scalar tmDot = mDot[ip];
                IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot *
                                             UbcVec[ip * SPATIAL_DIM + i]);

                //================================
                // stress: full stress
                //================================
                IP_FULL_STRESS_FIXED_VEL__();
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::
    assembleElemTermsBoundaryInletSpecifiedVelocityAndPressure_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_muEff;
    std::vector<scalar> ws_bcMultiplier;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& nodalSideUSTKFieldRef =
        model_->URef().nodeSideFieldRef().stkFieldRef();
    const auto& sideUSTKFieldRef = model_->URef().sideFieldRef().stkFieldRef();
    const auto& sidePSTKFieldRef = model_->pRef().sideFieldRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_bcMultiplier.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_bcMultiplier = &ws_bcMultiplier[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // set 1 on all nodes initially
                p_bcMultiplier[ni] = 1.0;

                // gather vectors
                const label offSet = ni * SPATIAL_DIM;

                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coords[offSet + j] = coords[j];
                    p_U[offSet + j] = U[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);

                // set 0 the boundary nodes
                p_bcMultiplier[ic] = 0.0;

                // gather vectors
                scalar* U = stk::mesh::field_data(nodalSideUSTKFieldRef, node);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    // overwrite velocity value at boundary with the dirichlet
                    // value
                    p_U[ic * SPATIAL_DIM + i] = U[i];
                }
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // calc vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                // FIXME: [2024-10-30] Flow reversal check is not applicable
                // if velocity is specified at boundary.
                if (rfflag[ip] == 1)
                {
                    errorMsg("invalid outflow at supersonic inlet");
                }
                else
                {
                    //================================
                    // advection: entering domain
                    //================================
                    const scalar tmDot = mDot[ip];
                    IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot *
                                                 UbcVec[ip * SPATIAL_DIM + i]);

                    //================================
                    // stress: full stress
                    //================================
                    IP_FULL_STRESS_FIXED_VEL__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundaryOutletSpecifiedPressure_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& sidePSTKFieldRef = model_->pRef().sideFieldRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // zero out vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                if (rfflag[ip] == 1)
                {
                    // slip-wall
                }
                else
                {
                    const scalar tmDot = mDot[ip];
                    IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(
                        tmDot * p_U[SPATIAL_DIM * nearestNode + i]);
                    IP_ZERO_NORMAL_STRESS__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundaryOutletOutflow_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // zero out vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                if (rfflag[ip] == 1)
                {
                    // slip-wall
                }
                else
                {
                    const scalar tmDot = mDot[ip];
                    IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(
                        tmDot * p_U[SPATIAL_DIM * nearestNode + i]);
                    IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundarySymmetry_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // vectors
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nBoundaryFaces = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundaryFaces;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* faceElemOrds =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = faceElemOrds[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elem_node_rels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elem_node_rels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // form unit normal
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                // Symmetry: assemble complete viscous stress without
                // projection. The full stress contribution is assembled at the
                // boundary, including both gradient and transpose terms:
                // τ = μ(∇U + ∇U^T) - (2/3)μ(∇·U)I
                // Note: The normal component of the momentum residual will be
                // zeroed later in applySymmetryConditions_() to enforce the
                // symmetry constraint (U·n = 0)
                IP_FULL_STRESS__();
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundarySlipWall_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    // Slip walls use the same implementation as symmetry planes
    // but with different stress formulation (zero tangential stress)
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // vectors
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nBoundaryFaces = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundaryFaces;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* faceElemOrds =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = faceElemOrds[0];

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elem_node_rels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elem_node_rels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // form unit normal
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];
                }

                // Slip wall: zero tangential stress with clipped transpose
                // for enhanced diagonal dominance
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundaryWallNoSlip_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    if (domain->turbulence_.option_ == turbulenceOption::laminar)
    {
        // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM
        // and nodesPerElem*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values; both boundary and opposing surface
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_nx = &nx[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_coords;
        std::vector<scalar> ws_p;
        std::vector<scalar> ws_muEff;
        std::vector<scalar> ws_bcMultiplier;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
        const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
        const auto& nodalSideUSTKFieldRef =
            model_->URef().nodeSideFieldRef().stkFieldRef();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // define vector of parent topos; should always be UNITY in size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();

        // shifted ip's for gradients?
        const bool isUGradientShifted = model_->URef().isGradientShifted();

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // extract connected element topology
            sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
            STK_ThrowAssert(parentTopo.size() == 1);
            stk::topology theElemTopo = parentTopo[0];

            // volume master element
            MasterElement* meSCS =
                accel::MasterElementRepo::get_surface_master_element(
                    theElemTopo);
            const label nodesPerElement = meSCS->nodesPerElement_;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // resize some things; matrix related
            const label lhsSize =
                nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
            const label rhsSize = nodesPerElement * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerElement);

            // algorithm related; element/face
            ws_U.resize(nodesPerElement * SPATIAL_DIM);
            ws_coords.resize(nodesPerElement * SPATIAL_DIM);
            ws_muEff.resize(nodesPerSide);
            ws_p.resize(nodesPerElement);
            ws_bcMultiplier.resize(nodesPerElement);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
            ws_det_j.resize(numScsBip);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_p = &ws_p[0];
            scalar* p_bcMultiplier = &ws_bcMultiplier[0];
            scalar* p_coords = &ws_coords[0];
            scalar* p_muEff = &ws_muEff[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];
            scalar* p_dndx = &ws_dndx[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    p_lhs[p] = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    p_rhs[p] = 0.0;
                }

                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                // extract the connected element to this exposed face; should be
                // single in size!
                stk::mesh::Entity const* faceElemRels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);

                // get element; its face ordinal number and populate
                // face_node_ordinals
                stk::mesh::Entity element = faceElemRels[0];
                const label faceOrdinal =
                    bulkData.begin_element_ordinals(side)[0];

                // mapping from ip to nodes for this ordinal;
                const label* ipNodeMap =
                    meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

                // populate faceNodeOrdinals
                const label* faceNodeOrdinals =
                    meSCS->side_node_ordinals(faceOrdinal);

                //==========================================
                // gather nodal data off of element
                //==========================================
                stk::mesh::Entity const* elemNodeRels =
                    bulkData.begin_nodes(element);
                label numNodes = bulkData.num_nodes(element);

                // sanity check on num nodes
                STK_ThrowAssert(numNodes == nodesPerElement);
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = elemNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                    p_bcMultiplier[ni] = 1.0;

                    // gather vectors
                    const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    const scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                        p_coords[offSet + j] = coords[j];
                    }
                }

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    const label ic = faceNodeOrdinals[ni];

                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);

                    // set 0 the boundary nodes
                    p_bcMultiplier[ic] = 0.0;

                    // gather vectors
                    scalar* U =
                        stk::mesh::field_data(nodalSideUSTKFieldRef, node);

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        // overwrite velocity value at boundary with the
                        // dirichlet value
                        p_U[ic * SPATIAL_DIM + i] = U[i];
                    }
                }

                // compute dndx
                scalar scs_error = 0.0;
                if (isUGradientShifted)
                {
                    meSCS->shifted_face_grad_op(1,
                                                faceOrdinal,
                                                &p_coords[0],
                                                &p_dndx[0],
                                                &ws_det_j[0],
                                                &scs_error);
                }
                else
                {
                    meSCS->face_grad_op(1,
                                        faceOrdinal,
                                        &p_coords[0],
                                        &p_dndx[0],
                                        &ws_det_j[0],
                                        &scs_error);
                }

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // offset for bip area vector and types of shape function
                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF = ip * nodesPerSide;

                    // zero out vector quantities
                    scalar asq = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[faceOffSet + j];
                        asq += axj * axj;
                    }
                    const scalar amag = std::sqrt(asq);

                    // interpolate to bip
                    scalar muEffBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF + ic];

                        muEffBip += r_vel * p_muEff[ic];
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_nx[i] = areaVec[faceOffSet + i] / amag;
                    }

                    //================================
                    // stress: zero normal stress
                    //================================
                    IP_ZERO_NORMAL_STRESS_FIXED_VEL__();
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
        // space for LHS/RHS; nodesPerSide*SPATIAL_DIM*nodesPerSide*SPATIAL_DIM
        // and nodesPerSide*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // bip values
        std::vector<scalar> uBip(SPATIAL_DIM);
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_uBip = &uBip[0];
        scalar* p_nx = &nx[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_p;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& sideUSTKFieldRef =
            model_->URef().sideFieldRef().stkFieldRef();
        const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
        const auto& uWallCoeffsSTKFieldRef =
            model_->uWallCoeffsRef().stkFieldRef();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);

        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // mapping from ip to nodes for this ordinal; face perspective (use
            // with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            // resize some things; matrix related
            const label lhsSize =
                nodesPerSide * SPATIAL_DIM * nodesPerSide * SPATIAL_DIM;
            const label rhsSize = nodesPerSide * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerSide);

            // algorithm related; element
            ws_U.resize(nodesPerSide * SPATIAL_DIM);
            ws_p.resize(nodesPerSide);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_p = &ws_p[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nBoundaryFaces =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundaryFaces;
                 ++iSide)
            {
                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    p_lhs[p] = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    p_rhs[p] = 0.0;
                }

                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* face_node_rels =
                    bulkData.begin_nodes(side);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = face_node_rels[ni];
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                    // gather vectors
                    scalar* U = stk::mesh::field_data(USTKFieldRef, node);

                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                    }
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const scalar* uWallCoeffsBip =
                    stk::mesh::field_data(uWallCoeffsSTKFieldRef, side);
                const scalar* UbcVec =
                    stk::mesh::field_data(sideUSTKFieldRef, side);

                // loop over face nodes
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label offSetAreaVec = ip * SPATIAL_DIM;
                    const label offSetSF = ip * nodesPerSide;

                    const label nearestNode = faceIpNodeMap[ip];

                    // zero out vector quantities; squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] = 0.0;

                        const scalar axj = areaVec[offSetAreaVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    // form unit normal
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_nx[j] = areaVec[offSetAreaVec + j] / aMag;
                    }

                    // interpolate to bip
                    // TODO: [2024-10-31] Redundant interpolation; cleanup
                    // candidate
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF + ic];

                        const label offSet = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_uBip[j] += r_vel * p_U[offSet + j];
                        }
                    }

                    //================================
                    // stress: zero normal stress
                    //================================
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        label indexR = nearestNode * SPATIAL_DIM + i;
                        label rowR = indexR * nodesPerSide * SPATIAL_DIM;

                        scalar uiTan = 0.0;
                        scalar uiBcTan = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar ninj = p_nx[i] * p_nx[j];
                            if (i == j)
                            {
                                const scalar om_nini = 1.0 - ninj;
                                uiTan += om_nini * p_uBip[j];
                                uiBcTan +=
                                    om_nini * UbcVec[ip * SPATIAL_DIM + j];
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    p_lhs[rowR + ic * SPATIAL_DIM + i] +=
                                        uWallCoeffsBip[ip] * om_nini *
                                        p_velocity_face_shape_function
                                            [offSetSF + ic];
                                }
                            }
                            else
                            {
                                // WARNING: [2024-11-25] Transpose term is not
                                // assembled in the legacy implementation
                                uiTan -= ninj * p_uBip[j];
                                uiBcTan -= ninj * UbcVec[ip * SPATIAL_DIM + j];
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    p_lhs[rowR + ic * SPATIAL_DIM + j] -=
                                        uWallCoeffsBip[ip] * ninj *
                                        p_velocity_face_shape_function
                                            [offSetSF + ic];
                                }
                            }
                        }
                        p_rhs[indexR] -= uWallCoeffsBip[ip] * (uiTan - uiBcTan);
                    }
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundaryOpening_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& sidePSTKFieldRef = model_->pRef().sideFieldRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();
    const auto& sideFlowDirectionSTKFieldRef =
        model_->URef().sideFlowDirectionFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* dir =
                stk::mesh::field_data(sideFlowDirectionSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // zero out vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to scs point; operate on saved off ws_field
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                }

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];

                    const label icNdim = inn * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r_vel * p_U[icNdim + j];
                    }
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                const scalar tmDot = mDot[ip];

                if (tmDot > 0.0) // outflow
                {
                    IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(
                        tmDot * p_U[SPATIAL_DIM * nearestNode + i]);
                    IP_ZERO_NORMAL_STRESS__();
                }
                else // inflow
                {
                    scalar num = 0.0;
                    scalar den = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const scalar d = dir[SPATIAL_DIM * ip + i];
                        num += p_uBip[i] * p_nx[i];
                        den += d * p_nx[i];
                    }
                    const scalar Iu_ipI = num / (den + SMALL);
                    IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot * Iu_ipI *
                                                 dir[SPATIAL_DIM * ip + i]);
                    IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleElemTermsBoundaryInletSpecifiedPressure_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];
    scalar* p_uBip = &uBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& sidePSTKFieldRef = model_->pRef().sideFieldRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();
    const auto& sideFlowDirectionSTKFieldRef =
        model_->URef().sideFlowDirectionFieldRef().stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        model_->URef().reversalFlagRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_muEff.resize(nodesPerSide);
        ws_p.resize(nodesPerElement);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* dir =
                stk::mesh::field_data(sideFlowDirectionSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // zero out vector quantities
                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                // interpolate to scs point; operate on saved off ws_field
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                }

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];

                    const label icNdim = inn * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r_vel * p_U[icNdim + j];
                    }
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                if (rfflag[ip] == 1)
                {
                    // slip-wall
                }
                else
                {
                    scalar num = 0.0;
                    scalar den = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const scalar d = dir[SPATIAL_DIM * ip + i];
                        num += p_uBip[i] * p_nx[i];
                        den += d * p_nx[i];
                    }

                    const scalar tmDot = mDot[ip];
                    const scalar Iu_ipI = num / (den + SMALL);

                    // implicit advective part
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF + ic];

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            label indexR = nearestNode * SPATIAL_DIM + i;
                            label rowR = indexR * nodesPerElement * SPATIAL_DIM;

                            const scalar di = dir[SPATIAL_DIM * ip + i];
                            const scalar ni = p_nx[i];

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar nj = p_nx[j];
                                if (i != j)
                                {
                                    p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                        tmDot / (den + SMALL) * r_vel * di * nj;
                                }
                                else
                                {
                                    p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                        std::max(tmDot / (den + SMALL) * r_vel *
                                                     di * ni,
                                                 0.0);
                                }
                            }
                        }
                    }

                    IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot * Iu_ipI *
                                                 dir[SPATIAL_DIM * ip + i]);
                    IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::
    assembleElemTermsBoundaryInletSpecifiedTotalPressure_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    if (!domain->isMaterialCompressible())
    {
        auto& mesh = model_->meshRef();

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        // hard-coded
        const scalar comp = 0.0;

        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM
        // and nodesPerElem*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values; both boundary and opposing surface
        std::vector<scalar> uBip(SPATIAL_DIM);
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_nx = &nx[0];
        scalar* p_uBip = &uBip[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_coords;
        std::vector<scalar> ws_p;
        std::vector<scalar> ws_muEff;
        std::vector<scalar> ws_rho;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
        const auto& sidePSTKFieldRef =
            model_->pRef().sideFieldRef().stkFieldRef();
        const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
        const auto& mDotSideSTKFieldRef =
            model_->mDotRef().sideFieldRef().stkFieldRef();
        const auto& sideFlowDirectionSTKFieldRef =
            model_->URef().sideFlowDirectionFieldRef().stkFieldRef();
        const auto& reversalFlowFlagSTKFieldRef =
            model_->URef().reversalFlagRef().stkFieldRef();

        // Get thermal fields for linearization
        const auto& rhoSTKFieldRef = model_->rhoRef().stkFieldRef();

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // define vector of parent topos; should always be UNITY in size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();
        const bool isUGradientShifted = model_->URef().isGradientShifted();

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // extract connected element topology
            sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
            STK_ThrowAssert(parentTopo.size() == 1);
            stk::topology theElemTopo = parentTopo[0];

            // volume master element
            MasterElement* meSCS =
                accel::MasterElementRepo::get_surface_master_element(
                    theElemTopo);
            const label nodesPerElement = meSCS->nodesPerElement_;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // resize some things; matrix related
            const label lhsSize =
                nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
            const label rhsSize = nodesPerElement * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerElement);

            // algorithm related; element/face
            ws_U.resize(nodesPerElement * SPATIAL_DIM);
            ws_coords.resize(nodesPerElement * SPATIAL_DIM);
            ws_p.resize(nodesPerElement);
            ws_muEff.resize(nodesPerSide);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
            ws_det_j.resize(numScsBip);
            ws_rho.resize(nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_p = &ws_p[0];
            scalar* p_coords = &ws_coords[0];
            scalar* p_muEff = &ws_muEff[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];
            scalar* p_dndx = &ws_dndx[0];
            scalar* p_rho = &ws_rho[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    p_lhs[p] = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    p_rhs[p] = 0.0;
                }

                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* mDot =
                    stk::mesh::field_data(mDotSideSTKFieldRef, side);
                const scalar* dir =
                    stk::mesh::field_data(sideFlowDirectionSTKFieldRef, side);
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const label* rfflag =
                    stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

                // extract the connected element to this exposed face; should be
                // single in size!
                stk::mesh::Entity const* faceElemRels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);

                // get element; its face ordinal number and populate
                // face_node_ordinals
                stk::mesh::Entity element = faceElemRels[0];
                const label faceOrdinal =
                    bulkData.begin_element_ordinals(side)[0];

                // mapping from ip to nodes for this ordinal;
                const label* ipNodeMap =
                    meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

                // populate faceNodeOrdinals
                const label* faceNodeOrdinals =
                    meSCS->side_node_ordinals(faceOrdinal);

                //==========================================
                // gather nodal data off of element
                //==========================================
                stk::mesh::Entity const* elemNodeRels =
                    bulkData.begin_nodes(element);
                label numNodes = bulkData.num_nodes(element);

                // sanity check on num nodes
                STK_ThrowAssert(numNodes == nodesPerElement);
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = elemNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                    // gather vectors
                    scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                        p_coords[offSet + j] = coords[j];
                    }
                }

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);
                    p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                }

                // compute dndx for velocity gradient
                scalar scs_error = 0.0;
                if (isUGradientShifted)
                {
                    meSCS->shifted_face_grad_op(1,
                                                faceOrdinal,
                                                &p_coords[0],
                                                &p_dndx[0],
                                                &ws_det_j[0],
                                                &scs_error);
                }
                else
                {
                    meSCS->face_grad_op(1,
                                        faceOrdinal,
                                        &p_coords[0],
                                        &p_dndx[0],
                                        &ws_det_j[0],
                                        &scs_error);
                }

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // offset for bip area vector and types of shape function
                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF = ip * nodesPerSide;

                    // zero out vector quantities
                    scalar asq = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[faceOffSet + j];
                        asq += axj * axj;
                    }
                    const scalar amag = std::sqrt(asq);

                    // interpolate to scs point; operate on saved off ws_field
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar muEffBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF + ic];

                        muEffBip += r_vel * p_muEff[ic];

                        const label icNdim = inn * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_uBip[j] += r_vel * p_U[icNdim + j];
                        }
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_nx[i] = areaVec[faceOffSet + i] / amag;
                    }

                    if (rfflag[ip] == 1)
                    {
                        // slip-wall
                    }
                    else
                    {
                        scalar num = 0.0;
                        scalar den = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar d = dir[SPATIAL_DIM * ip + i];
                            num += p_uBip[i] * p_nx[i];
                            den += d * p_nx[i];
                        }

                        const scalar tmDot = mDot[ip];
                        const scalar Iu_ipI = num / (den + SMALL);

                        // implicit advective part
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const label inn = faceNodeOrdinals[ic];

                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                label indexR = nearestNode * SPATIAL_DIM + i;
                                label rowR =
                                    indexR * nodesPerElement * SPATIAL_DIM;

                                const scalar di = dir[SPATIAL_DIM * ip + i];
                                const scalar ni = p_nx[i];

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar nj = p_nx[j];
                                    if (i != j)
                                    {
                                        p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                            tmDot / (den + SMALL) * r_vel * di *
                                            nj;
                                    }
                                    else
                                    {
                                        p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                            std::max(tmDot / (den + SMALL) *
                                                         r_vel * di * ni,
                                                     0.0);
                                    }
                                }
                            }
                        }

                        IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot * Iu_ipI *
                                                     dir[SPATIAL_DIM * ip + i]);
                        IP_DIRECTIONAL_STRESS__();

                        //==========================================
                        // Total Pressure Linearization
                        //==========================================

                        scalar rhoBip = 0.0;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];
                            rhoBip += r_vel * p_rho[ic];
                        }

                        // linearization coefficient: coef = -rho
                        const scalar coef = -rhoBip;

                        // add linearization: ∂P_static/∂U coupling to lhs
                        // matrix for total pressure inlet boundary condition,
                        // the static pressure depends on velocity magnitude
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const label inn = faceNodeOrdinals[ic];
                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                const label indexR =
                                    nearestNode * SPATIAL_DIM + i;
                                const label rowR =
                                    indexR * nodesPerElement * SPATIAL_DIM;
                                const scalar axi = areaVec[faceOffSet + i];

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                        coef * p_uBip[j] * axi * r_vel;
                                }
                            }
                        }
                    }
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
        auto& mesh = model_->meshRef();

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        // hard-coded
        const scalar comp = 1.0;

        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM
        // and nodesPerElem*SPATIAL_DIM
        std::vector<scalar> lhs;
        std::vector<scalar> rhs;
        std::vector<label> scratchIds;
        std::vector<scalar> scratchVals;
        std::vector<stk::mesh::Entity> connectedNodes;

        // ip values; both boundary and opposing surface
        std::vector<scalar> uBip(SPATIAL_DIM);
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_nx = &nx[0];
        scalar* p_uBip = &uBip[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_coords;
        std::vector<scalar> ws_p;
        std::vector<scalar> ws_muEff;
        std::vector<scalar> ws_T;
        std::vector<scalar> ws_cp;
        std::vector<scalar> ws_rho;

        // master element
        std::vector<scalar> ws_velocity_face_shape_function;
        std::vector<scalar> ws_pressure_face_shape_function;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_det_j;

        // Get fields
        const auto& USTKFieldRef = model_->URef().stkFieldRef();
        const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
        const auto& sidePSTKFieldRef =
            model_->pRef().sideFieldRef().stkFieldRef();
        const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
        const auto& mDotSideSTKFieldRef =
            model_->mDotRef().sideFieldRef().stkFieldRef();
        const auto& sideFlowDirectionSTKFieldRef =
            model_->URef().sideFlowDirectionFieldRef().stkFieldRef();
        const auto& reversalFlowFlagSTKFieldRef =
            model_->URef().reversalFlagRef().stkFieldRef();

        // Get thermal fields for linearization
        const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
        const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
        const auto& rhoSTKFieldRef = model_->rhoRef().stkFieldRef();

        // Get material properties for linearization
        const scalar M =
            domain->materialRef()
                .thermodynamicProperties_.equationOfState_.molarMass_;
        const scalar R = thermoModel::universalGasConstant_;
        const scalar Rs = R / M;

        // Get geometric fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // define vector of parent topos; should always be UNITY in size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        // shifted ip's for fields?
        const bool isUShifted = model_->URef().isShifted();
        const bool isPShifted = model_->pRef().isShifted();
        const bool isUGradientShifted = model_->URef().isGradientShifted();

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);
        for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
             ib != sideBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideBucket = **ib;

            // extract connected element topology
            sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
            STK_ThrowAssert(parentTopo.size() == 1);
            stk::topology theElemTopo = parentTopo[0];

            // volume master element
            MasterElement* meSCS =
                accel::MasterElementRepo::get_surface_master_element(
                    theElemTopo);
            const label nodesPerElement = meSCS->nodesPerElement_;

            // face master element
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
                    sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // resize some things; matrix related
            const label lhsSize =
                nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
            const label rhsSize = nodesPerElement * SPATIAL_DIM;
            lhs.resize(lhsSize);
            rhs.resize(rhsSize);
            scratchIds.resize(rhsSize);
            scratchVals.resize(rhsSize);
            connectedNodes.resize(nodesPerElement);

            // algorithm related; element/face
            ws_U.resize(nodesPerElement * SPATIAL_DIM);
            ws_coords.resize(nodesPerElement * SPATIAL_DIM);
            ws_p.resize(nodesPerElement);
            ws_muEff.resize(nodesPerSide);
            ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
            ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
            ws_det_j.resize(numScsBip);
            ws_T.resize(nodesPerSide);
            ws_cp.resize(nodesPerSide);
            ws_rho.resize(nodesPerSide);

            // pointers
            scalar* p_lhs = &lhs[0];
            scalar* p_rhs = &rhs[0];
            scalar* p_U = &ws_U[0];
            scalar* p_p = &ws_p[0];
            scalar* p_coords = &ws_coords[0];
            scalar* p_muEff = &ws_muEff[0];
            scalar* p_velocity_face_shape_function =
                &ws_velocity_face_shape_function[0];
            scalar* p_pressure_face_shape_function =
                &ws_pressure_face_shape_function[0];
            scalar* p_dndx = &ws_dndx[0];
            scalar* p_T = &ws_T[0];
            scalar* p_cp = &ws_cp[0];
            scalar* p_rho = &ws_rho[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_velocity_face_shape_function[0]);
            }

            if (isPShifted)
            {
                meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_pressure_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    p_lhs[p] = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    p_rhs[p] = 0.0;
                }

                // get face
                stk::mesh::Entity side = sideBucket[iSide];

                // pointer to face data
                const scalar* mDot =
                    stk::mesh::field_data(mDotSideSTKFieldRef, side);
                const scalar* dir =
                    stk::mesh::field_data(sideFlowDirectionSTKFieldRef, side);
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const scalar* pbc =
                    stk::mesh::field_data(sidePSTKFieldRef, side);
                const label* rfflag =
                    stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

                // extract the connected element to this exposed face; should be
                // single in size!
                stk::mesh::Entity const* faceElemRels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);

                // get element; its face ordinal number and populate
                // face_node_ordinals
                stk::mesh::Entity element = faceElemRels[0];
                const label faceOrdinal =
                    bulkData.begin_element_ordinals(side)[0];

                // mapping from ip to nodes for this ordinal;
                const label* ipNodeMap =
                    meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

                // populate faceNodeOrdinals
                const label* faceNodeOrdinals =
                    meSCS->side_node_ordinals(faceOrdinal);

                //==========================================
                // gather nodal data off of element
                //==========================================
                stk::mesh::Entity const* elemNodeRels =
                    bulkData.begin_nodes(element);
                label numNodes = bulkData.num_nodes(element);

                // sanity check on num nodes
                STK_ThrowAssert(numNodes == nodesPerElement);
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = elemNodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    // gather scalars
                    p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                    // gather vectors
                    scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                        p_coords[offSet + j] = coords[j];
                    }
                }

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < numSideNodes; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_muEff[ni] =
                        *stk::mesh::field_data(muEffSTKFieldRef, node);
                    p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                    p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                    p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                }

                // compute dndx for velocity gradient
                scalar scs_error = 0.0;
                if (isUGradientShifted)
                {
                    meSCS->shifted_face_grad_op(1,
                                                faceOrdinal,
                                                &p_coords[0],
                                                &p_dndx[0],
                                                &ws_det_j[0],
                                                &scs_error);
                }
                else
                {
                    meSCS->face_grad_op(1,
                                        faceOrdinal,
                                        &p_coords[0],
                                        &p_dndx[0],
                                        &ws_det_j[0],
                                        &scs_error);
                }

                // loop over boundary ips
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    const label nearestNode = ipNodeMap[ip];

                    // offset for bip area vector and types of shape function
                    const label faceOffSet = ip * SPATIAL_DIM;
                    const label offSetSF = ip * nodesPerSide;

                    // zero out vector quantities
                    scalar asq = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[faceOffSet + j];
                        asq += axj * axj;
                    }
                    const scalar amag = std::sqrt(asq);

                    // interpolate to scs point; operate on saved off ws_field
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] = 0.0;
                    }

                    // interpolate to bip
                    scalar muEffBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const label inn = faceNodeOrdinals[ic];

                        const scalar r_vel =
                            p_velocity_face_shape_function[offSetSF + ic];

                        muEffBip += r_vel * p_muEff[ic];

                        const label icNdim = inn * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_uBip[j] += r_vel * p_U[icNdim + j];
                        }
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_nx[i] = areaVec[faceOffSet + i] / amag;
                    }

                    if (rfflag[ip] == 1)
                    {
                        // slip-wall
                    }
                    else
                    {
                        scalar num = 0.0;
                        scalar den = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar d = dir[SPATIAL_DIM * ip + i];
                            num += p_uBip[i] * p_nx[i];
                            den += d * p_nx[i];
                        }

                        const scalar tmDot = mDot[ip];
                        const scalar Iu_ipI = num / (den + SMALL);

                        // implicit advective part
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const label inn = faceNodeOrdinals[ic];

                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                label indexR = nearestNode * SPATIAL_DIM + i;
                                label rowR =
                                    indexR * nodesPerElement * SPATIAL_DIM;

                                const scalar di = dir[SPATIAL_DIM * ip + i];
                                const scalar ni = p_nx[i];

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar nj = p_nx[j];
                                    if (i != j)
                                    {
                                        p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                            tmDot / (den + SMALL) * r_vel * di *
                                            nj;
                                    }
                                    else
                                    {
                                        p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                            std::max(tmDot / (den + SMALL) *
                                                         r_vel * di * ni,
                                                     0.0);
                                    }
                                }
                            }
                        }

                        IP_EXPLICIT_ADVECTIVE_FLUX__(tmDot * Iu_ipI *
                                                     dir[SPATIAL_DIM * ip + i]);
                        IP_DIRECTIONAL_STRESS__();

                        //==========================================
                        // Total Pressure Linearization
                        //==========================================
                        // Interpolate temperature and Cp to boundary IP
                        scalar TBip = 0.0;
                        scalar cpBip = 0.0;
                        scalar rhoBip = 0.0;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];
                            TBip += r_vel * p_T[ic];
                            cpBip += r_vel * p_cp[ic];
                            rhoBip += r_vel * p_rho[ic];
                        }

                        // Compute gamma = cp / cv = cp / (cp - R)
                        const scalar gammaBip = cpBip / (cpBip - Rs);

                        // Calculate total temperature locally from static
                        // temperature and velocity T_total = T_static +
                        // U²/(2*Cp)
                        scalar U2 = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            U2 += p_uBip[j] * p_uBip[j];
                        }
                        const scalar T0Bip = TBip + U2 / (2.0 * cpBip);

                        // calculate local mach number
                        const scalar Ma =
                            std::sqrt(U2) / std::sqrt(gammaBip * Rs * TBip);

                        // add linearization: ∂P_static/∂U coupling to lhs
                        // matrix for total pressure inlet boundary condition,
                        // the static pressure depends on velocity magnitude
                        // through the isentropic relation. The linearization is
                        // with respect to velocity components at the boundary
                        // IP.
                        scalar coef = 0;
                        if (Ma < 0.3)
                        {
                            // linearization coefficient same as incompressible:
                            // coef = -rho
                            coef = -rhoBip;
                        }
                        else
                        {
                            // linearization coefficient: coef =
                            // -(gamma/(gamma-1)) / (T0 * cp) * p_abs
                            coef = -gammaBip / (gammaBip - 1.0) /
                                   (T0Bip * cpBip) * pbc[ip];
                        }

                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const label inn = faceNodeOrdinals[ic];
                            const scalar r_vel =
                                p_velocity_face_shape_function[offSetSF + ic];

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                const label indexR =
                                    nearestNode * SPATIAL_DIM + i;
                                const label rowR =
                                    indexR * nodesPerElement * SPATIAL_DIM;
                                const scalar axi = areaVec[faceOffSet + i];

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_lhs[rowR + inn * SPATIAL_DIM + j] +=
                                        coef * p_uBip[j] * axi * r_vel;
                                }
                            }
                        }
                    }
                }

                this->applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void navierStokesAssembler::
    assembleElemTermsBoundaryOutletSpecifiedMassFlowRate_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS; nodesPerElem*SPATIAL_DIM*nodesPerElem*SPATIAL_DIM and
    // nodesPerElem*SPATIAL_DIM
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // ip values; both boundary and opposing surface
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_muEff;

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_pressure_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& USTKFieldRef = model_->URef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& muEffSTKFieldRef = model_->muEffRef().stkFieldRef();
    const auto& mDotSideSTKFieldRef =
        model_->mDotRef().sideFieldRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();
    const bool isPShifted = model_->pRef().isShifted();

    // shifted ip's for gradients?
    const bool isUGradientShifted = model_->URef().isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element/face
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_pressure_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_p = &ws_p[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_pressure_face_shape_function =
            &ws_pressure_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_pressure_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_pressure_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            // mapping from ip to nodes for this ordinal;
            const label* ipNodeMap =
                meSCS->ipNodeMap(faceOrdinal); // use with elem_node_rels

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                // offset for bip area vector and types of shape function
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                // interpolate to scs point; operate on saved off ws_field
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                }

                // interpolate to bip
                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];

                    const label icNdim = inn * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r_vel * p_U[icNdim + j];
                    }
                }

                const scalar tmDot = mDot[ip];
                IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(
                    tmDot * p_U[SPATIAL_DIM * nearestNode + i]);
                IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__();
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

#undef IP_EXPLICIT_ADVECTIVE_FLUX__
#undef IP_IMPLICIT_ADVECTIVE_FLUX__
#undef IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__
#undef IP_FULL_STRESS__
#undef IP_FULL_STRESS_FIXED_VEL__
#undef IP_ZERO_NORMAL_STRESS__
#undef IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__
#undef IP_ZERO_TANGENTIAL_STRESS_CLIPPED_GRADU_TRANSPOSE__
#undef IP_ZERO_NORMAL_STRESS_FIXED_VEL__
#undef IP_DIRECTIONAL_STRESS__

} // namespace accel
