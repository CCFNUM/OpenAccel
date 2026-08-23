// File       : solidDisplacementAssemblerElemTerms.cpp
// Created    : Thu Dec 04 2025 10:15:23 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Element interior terms for linear elastic stress
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"

#include "linear_elasticity.hpp"
#include "sfem_GeneratedNeoHookeanOgden_element_api.hpp"
#include "sfem_GeneratedModifiedMooneyRivlin_c_abi.hpp"

namespace accel
{

#if SPATIAL_DIM == 2
namespace
{

void assembleQuad4LinearElasticity(
    const std::vector<std::vector<geom_t>>& coordinates,
    const scalar mu,
    const scalar lambda,
    std::vector<scalar>& lhs)
{
    constexpr scalar invSqrt3 = 0.57735026918962576451;
    const scalar gaussPoints[2] = {-invSqrt3, invSqrt3};
    const scalar c11 = lambda + 2.0 * mu;
    const scalar c12 = lambda;
    const scalar c33 = mu;

    for (const scalar xi : gaussPoints)
    {
        for (const scalar eta : gaussPoints)
        {
            scalar dN_dxi[4] = {
                -0.25 * (1.0 - eta),
                0.25 * (1.0 - eta),
                0.25 * (1.0 + eta),
                -0.25 * (1.0 + eta)};
            scalar dN_deta[4] = {
                -0.25 * (1.0 - xi),
                -0.25 * (1.0 + xi),
                0.25 * (1.0 + xi),
                0.25 * (1.0 - xi)};

            scalar j00 = 0.0;
            scalar j01 = 0.0;
            scalar j10 = 0.0;
            scalar j11 = 0.0;
            for (label node = 0; node < 4; ++node)
            {
                const scalar x = coordinates[0][node];
                const scalar y = coordinates[1][node];
                j00 += dN_dxi[node] * x;
                j01 += dN_deta[node] * x;
                j10 += dN_dxi[node] * y;
                j11 += dN_deta[node] * y;
            }

            const scalar detJ = j00 * j11 - j01 * j10;
            STK_ThrowRequireMsg(
                detJ > 0.0,
                "Invalid QUAD4 FEM solid element: non-positive Jacobian");

            const scalar invJ00 = j11 / detJ;
            const scalar invJ01 = -j01 / detJ;
            const scalar invJ10 = -j10 / detJ;
            const scalar invJ11 = j00 / detJ;

            scalar dN_dx[4];
            scalar dN_dy[4];
            for (label node = 0; node < 4; ++node)
            {
                dN_dx[node] = invJ00 * dN_dxi[node] +
                              invJ10 * dN_deta[node];
                dN_dy[node] = invJ01 * dN_dxi[node] +
                              invJ11 * dN_deta[node];
            }

            for (label a = 0; a < 4; ++a)
            {
                const scalar Bax[3] = {dN_dx[a], 0.0, dN_dy[a]};
                const scalar Bay[3] = {0.0, dN_dy[a], dN_dx[a]};

                for (label b = 0; b < 4; ++b)
                {
                    const scalar Bbx[3] = {dN_dx[b], 0.0, dN_dy[b]};
                    const scalar Bby[3] = {0.0, dN_dy[b], dN_dx[b]};

                    const scalar kxx =
                        Bax[0] * (c11 * Bbx[0] + c12 * Bbx[1]) +
                        Bax[1] * (c12 * Bbx[0] + c11 * Bbx[1]) +
                        Bax[2] * c33 * Bbx[2];
                    const scalar kxy =
                        Bax[0] * (c11 * Bby[0] + c12 * Bby[1]) +
                        Bax[1] * (c12 * Bby[0] + c11 * Bby[1]) +
                        Bax[2] * c33 * Bby[2];
                    const scalar kyx =
                        Bay[0] * (c11 * Bbx[0] + c12 * Bbx[1]) +
                        Bay[1] * (c12 * Bbx[0] + c11 * Bbx[1]) +
                        Bay[2] * c33 * Bbx[2];
                    const scalar kyy =
                        Bay[0] * (c11 * Bby[0] + c12 * Bby[1]) +
                        Bay[1] * (c12 * Bby[0] + c11 * Bby[1]) +
                        Bay[2] * c33 * Bby[2];

                    const label rowX = a * SPATIAL_DIM;
                    const label rowY = rowX + 1;
                    const label colX = b * SPATIAL_DIM;
                    const label colY = colX + 1;
                    const scalar weight = detJ;

                    lhs[rowX * 8 + colX] += weight * kxx;
                    lhs[rowX * 8 + colY] += weight * kxy;
                    lhs[rowY * 8 + colX] += weight * kyx;
                    lhs[rowY * 8 + colY] += weight * kyy;
                }
            }
        }
    }
}

} // namespace

#endif

namespace
{

void assembleSfemNeoHookeanElement(
    const smesh::ElemType sfemElementType,
    const label nodesPerElement,
    const scalar lambda,
    const scalar mu,
    const std::vector<std::vector<scalar>>& coordinates,
    const std::vector<scalar>& displacement,
    std::vector<scalar>& rhs,
    std::vector<scalar>& lhs)
{
    const label dofsPerElement = nodesPerElement * SPATIAL_DIM;
    const label lhsSize = dofsPerElement * dofsPerElement;

    std::vector<scalar> gradient(dofsPerElement, 0.0);
    std::vector<scalar> hessian(lhsSize, 0.0);
    std::vector<const scalar*> coordinateStreams(dofsPerElement);
    std::vector<const scalar*> displacementStreams(dofsPerElement);
    std::vector<scalar*> gradientStreams(dofsPerElement);
    std::vector<scalar*> hessianStreams(lhsSize);

    // The generated SFEM kernels expect a structure-of-arrays interface.
    // We pass a single element batch here so OpenAccel keeps the same
    // element-wise assemble -> insert-to-global workflow as the linear FEM
    // path.
    for (label node = 0; node < nodesPerElement; ++node)
    {
        for (label dim = 0; dim < SPATIAL_DIM; ++dim)
        {
            const label dof = node * SPATIAL_DIM + dim;
            coordinateStreams[dof] = &coordinates[dim][node];
            displacementStreams[dof] = &displacement[dof];
            gradientStreams[dof] = &gradient[dof];
        }
    }

    for (label entry = 0; entry < lhsSize; ++entry)
        hessianStreams[entry] = &hessian[entry];

#if SPATIAL_DIM == 2
    const int gradientStatus =
        sfem::codegen::neohookean_ogden_gradient_2d_element_soa<scalar, 1>(
            sfemElementType,
            1,
            coordinateStreams.data(),
            lambda,
            mu,
            displacementStreams.data(),
            gradientStreams.data());
    const int hessianStatus =
        sfem::codegen::neohookean_ogden_hessian_2d_element_soa<scalar, 1>(
            sfemElementType,
            1,
            coordinateStreams.data(),
            lambda,
            mu,
            displacementStreams.data(),
            hessianStreams.data());
#else
    const int gradientStatus =
        sfem::codegen::neohookean_ogden_gradient_3d_element_soa<scalar, 1>(
            sfemElementType,
            1,
            coordinateStreams.data(),
            lambda,
            mu,
            displacementStreams.data(),
            gradientStreams.data());
    const int hessianStatus =
        sfem::codegen::neohookean_ogden_hessian_3d_element_soa<scalar, 1>(
            sfemElementType,
            1,
            coordinateStreams.data(),
            lambda,
            mu,
            displacementStreams.data(),
            hessianStreams.data());
#endif

    STK_ThrowRequireMsg(
        gradientStatus == SFEM_SUCCESS && hessianStatus == SFEM_SUCCESS,
        "sfem Neo-Hookean elemental kernel failed");

    lhs = std::move(hessian);
    for (label row = 0; row < dofsPerElement; ++row)
        rhs[row] = -gradient[row];
}

#if SPATIAL_DIM == 2
// Unlike the Neo-Hookean kernels, the generated modified Mooney-Rivlin
// kernels only expose the mesh-level C ABI (*_isoparametric_mesh_soa,
// declared in sfem_GeneratedModifiedMooneyRivlin_c_abi.hpp) -- there is no
// single-element "element_soa" convenience wrapper to mirror
// assembleSfemNeoHookeanElement's call above. We therefore present this one
// element as a trivial local "mesh" of `nodesPerElement` nodes, using the
// same self-referential local connectivity already built for the
// linear_elasticity_apply_aos fallback below (elementConnectivity[node] =
// &localNodeIds[node], localNodeIds[node] = node).
//
// The hessian is BSR (block sparse row) over that local node graph: it needs
// an explicit CSR (rowptr, colidx) pattern, and writes DIM x DIM,
// row-major blocks into `values` at each (row, col) graph entry via
// atomic += (so `values` must be pre-zeroed). For our local 1-element
// "mesh" every node is coupled to every other node, so the graph is fully
// dense: rowptr[i] = i * nodesPerElement, colidx lists 0..nodesPerElement-1
// for every row. The block for local nodes (i, j) then lands at
// values[(i * nodesPerElement + j) * SPATIAL_DIM * SPATIAL_DIM + bi *
// SPATIAL_DIM + bj], which we scatter into the dense, node-major
// `lhs[dofsPerElement * dofsPerElement]` OpenAccel expects.
void assembleSfemModifiedMooneyRivlinElement(
    const smesh::ElemType sfemElementType,
    const label nodesPerElement,
    const scalar c1,
    const scalar c2,
    const scalar kappa,
    const std::vector<std::vector<geom_t>>& coordinates,
    std::vector<idx_t*>& elementConnectivity,
    const std::vector<scalar>& displacement,
    std::vector<scalar>& rhs,
    std::vector<scalar>& lhs)
{
    const label dofsPerElement = nodesPerElement * SPATIAL_DIM;

    // SoA displacement/gradient-output views into the existing node-major
    // `displacement`/`gradient` buffers: component d of node `n` lives at
    // offset n * SPATIAL_DIM + d, i.e. exactly a stride-SPATIAL_DIM SoA
    // array starting at offset d -- no separate copy needed.
    std::vector<scalar> gradient(dofsPerElement, 0.0);
    const scalar* ux = displacement.data() + 0;
    const scalar* uy = displacement.data() + 1;
    scalar* outx = gradient.data() + 0;
    scalar* outy = gradient.data() + 1;

    std::vector<const geom_t*> pointStreams(SPATIAL_DIM);
    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
        pointStreams[dim] = coordinates[dim].data();

    const int gradientStatus = modified_mooney_rivlin_gradient_2d_isoparametric_mesh_soa(
        sfemElementType,
        1,
        nodesPerElement,
        elementConnectivity.data(),
        pointStreams.data(),
        c1,
        c2,
        kappa,
        SPATIAL_DIM,
        ux,
        uy,
        SPATIAL_DIM,
        outx,
        outy);

    // Dense CSR graph over the local pseudo-mesh's `nodesPerElement` nodes
    // (see comment above): every node coupled to every other node.
    std::vector<count_t> rowptr(nodesPerElement + 1);
    std::vector<idx_t> colidx(nodesPerElement * nodesPerElement);
    for (label i = 0; i < nodesPerElement; ++i)
    {
        rowptr[i] = i * nodesPerElement;
        for (label j = 0; j < nodesPerElement; ++j)
            colidx[i * nodesPerElement + j] = j;
    }
    rowptr[nodesPerElement] = nodesPerElement * nodesPerElement;

    std::vector<scalar> values(
        nodesPerElement * nodesPerElement * SPATIAL_DIM * SPATIAL_DIM, 0.0);

    const int hessianStatus = modified_mooney_rivlin_hessian_bsr_2d_isoparametric_mesh_soa(
        sfemElementType,
        1,
        nodesPerElement,
        elementConnectivity.data(),
        pointStreams.data(),
        c1,
        c2,
        kappa,
        SPATIAL_DIM,
        ux,
        uy,
        rowptr.data(),
        colidx.data(),
        values.data());

    STK_ThrowRequireMsg(
        gradientStatus == SFEM_SUCCESS && hessianStatus == SFEM_SUCCESS,
        "sfem modified Mooney-Rivlin elemental kernel failed");

    // Scatter the dense BSR blocks into OpenAccel's dense, node-major lhs.
    lhs.assign(dofsPerElement * dofsPerElement, 0.0);
    for (label i = 0; i < nodesPerElement; ++i)
    {
        for (label j = 0; j < nodesPerElement; ++j)
        {
            const scalar* block =
                &values[(i * nodesPerElement + j) * SPATIAL_DIM * SPATIAL_DIM];
            for (label bi = 0; bi < SPATIAL_DIM; ++bi)
            {
                const label row = i * SPATIAL_DIM + bi;
                for (label bj = 0; bj < SPATIAL_DIM; ++bj)
                {
                    const label col = j * SPATIAL_DIM + bj;
                    lhs[row * dofsPerElement + col] =
                        block[bi * SPATIAL_DIM + bj];
                }
            }
        }
    }

    for (label row = 0; row < dofsPerElement; ++row)
        rhs[row] = -gradient[row];
}
#endif // SPATIAL_DIM == 2

} // namespace

void solidDisplacementAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
    if (field_broker_->controlsRef().isCvfemSolidMechanics())
    {
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Workspace for LHS/RHS
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // Nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_D;
    std::vector<scalar> ws_E;
    std::vector<scalar> ws_nu;

    // Geometry related
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;

    // Workspace for neo-hookean model (only used if neoHookean)
    std::vector<scalar> ws_gradTens;

    // const-size containers
    std::vector<scalar> kd(SPATIAL_DIM * SPATIAL_DIM);

    // fill Kronecker delta
    for (label i = 0; i < SPATIAL_DIM; ++i)
    {
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            kd[i * SPATIAL_DIM + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // pointers to const-size arrays
    scalar* p_kd = &kd[0];

    // Get fields
    const auto& DSTKFieldRef = phi_->stkFieldRef();
    const auto& ESTKFieldRef = model_->ERef().stkFieldRef();
    const auto& nuSTKFieldRef = model_->nuRef().stkFieldRef();
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // Check if plane stress or plane strain
    const bool planeStress = domain->solidMechanics_.planeStress_;

    // Get the solid mechanics option (linear elastic vs simplified neo-hookean)
    const solidMechanicsOption solidMechOption =
        domain->solidMechanics_.option_;

    // Shifted integration points?
    const bool isShifted = phi_->isShifted();
    const bool isGradientShifted = phi_->isGradientShifted();

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // Extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        // Resize arrays
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_D.resize(nodesPerElement * SPATIAL_DIM);
        ws_E.resize(nodesPerElement);
        ws_nu.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);

        // Resize workspace for neo-hookean model if needed
        if (solidMechOption == solidMechanicsOption::neoHookean)
        {
            ws_gradTens.resize(SPATIAL_DIM * SPATIAL_DIM);
        }

        // Pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_D = &ws_D[0];
        scalar* p_E = &ws_E[0];
        scalar* p_nu = &ws_nu[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_gradTens =
            solidMechOption == solidMechanicsOption::neoHookean
                ? &ws_gradTens[0]
                : nullptr;

        // Extract shape functions for stress terms
        if (isShifted)
        {
            meSCS->shifted_shape_fcn(p_shape_function);
        }
        else
        {
            meSCS->shape_fcn(p_shape_function);
        }

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            stk::mesh::Entity elem = elementBucket[iElement];

            // Zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
                p_lhs[p] = 0.0;
            for (label p = 0; p < rhsSize; ++p)
                p_rhs[p] = 0.0;

            // Gather nodal data
            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(elem);
            label numNodes = bulkData.num_nodes(elem);
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];
                connectedNodes[ni] = node;

                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);
                const scalar* D = stk::mesh::field_data(DSTKFieldRef, node);

                p_E[ni] = *stk::mesh::field_data(ESTKFieldRef, node);
                p_nu[ni] = *stk::mesh::field_data(nuSTKFieldRef, node);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    p_D[ni * SPATIAL_DIM + i] = D[i];
                }
            }

            // Compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(1, p_coordinates, p_scs_areav, &scs_error);

            // compute dndx
            if (isGradientShifted)
            {
                meSCS->shifted_grad_op(1,
                                       &ws_coordinates[0],
                                       &ws_dndx[0],
                                       &ws_deriv[0],
                                       &ws_det_j[0],
                                       &scs_error);
            }
            else
            {
                meSCS->grad_op(1,
                               &ws_coordinates[0],
                               &ws_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scs_error);
            }

            // Loop over integration points
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // Interpolate Lame parameters to integration point
                scalar muIp = 0.0;
                scalar lambdaIp = 0.0;
                const label offSetSF = ip * nodesPerElement;

                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF + ic];
                    const scalar E = p_E[ic];
                    const scalar nu = p_nu[ic];

                    // Compute Lame parameters from E and nu
                    const scalar mu = E / (2.0 * (1.0 + nu));
                    scalar lambda;
                    if (planeStress)
                    {
                        lambda = nu * E / ((1.0 + nu) * (1.0 - nu));
                    }
                    else
                    {
                        lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    }

                    muIp += r * mu;
                    lambdaIp += r * lambda;
                }

                // ================================================================
                // Choose model based on solidMechanicsOption
                // ================================================================
                if (solidMechOption == solidMechanicsOption::linearElastic)
                {
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label icNdim = ic * SPATIAL_DIM;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const label indexL = il * SPATIAL_DIM + i;
                            const label indexR = ir * SPATIAL_DIM + i;

                            const label rowL =
                                indexL * nodesPerElement * SPATIAL_DIM;
                            const label rowR =
                                indexR * nodesPerElement * SPATIAL_DIM;

                            scalar lhs_riC_i = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    p_scs_areav[ip * SPATIAL_DIM + j];
                                const scalar Dxj = p_D[ic * SPATIAL_DIM + j];

                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;

                                // First mu term: -mu * dN/dx_j * A_j (diagonal
                                // contribution)
                                const scalar lhsfacDiff_i =
                                    -muIp * p_dndx[offSetDnDx + j] * axj;
                                lhs_riC_i += lhsfacDiff_i;

                                // Second mu term: -mu * dN/dx_i * A_j
                                // (off-diagonal contribution)
                                const scalar lhsfacDiff_j =
                                    -muIp * p_dndx[offSetDnDx + i] * axj;

                                p_lhs[rowL + icNdim + j] += lhsfacDiff_j;
                                p_lhs[rowR + icNdim + j] -= lhsfacDiff_j;

                                p_rhs[indexL] -= lhsfacDiff_j * Dxj;
                                p_rhs[indexR] += lhsfacDiff_j * Dxj;

                                // Lambda divergence term: -lambda * delta_ij *
                                // dN/dx_l * A_j
                                for (label l = 0; l < SPATIAL_DIM; ++l)
                                {
                                    const scalar lhsfacDiv =
                                        -lambdaIp * p_kd[i * SPATIAL_DIM + j] *
                                        p_dndx[offSetDnDx + l] * axj;

                                    p_lhs[rowL + icNdim + l] += lhsfacDiv;
                                    p_lhs[rowR + icNdim + l] -= lhsfacDiv;

                                    const scalar dxl =
                                        p_D[ic * SPATIAL_DIM + l];
                                    p_rhs[indexL] -= lhsfacDiv * dxl;
                                    p_rhs[indexR] += lhsfacDiv * dxl;
                                }
                            }

                            // Accumulated diagonal term
                            p_lhs[rowL + icNdim + i] += lhs_riC_i;
                            p_lhs[rowR + icNdim + i] -= lhs_riC_i;

                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];
                            p_rhs[indexL] -= lhs_riC_i * dxi;
                            p_rhs[indexR] += lhs_riC_i * dxi;
                        }
                    }
                }
                else if (solidMechOption ==
                         solidMechanicsOption::neoHookean)
                {
                    // Zero out gradient tensor
                    for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                    {
                        p_gradTens[i] = 0.0;
                    }

                    // Compute full lagged gradient tensor: F = I + grad(u)
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;
                                p_gradTens[i * SPATIAL_DIM + j] +=
                                    p_dndx[offSetDnDx + j] * dxi;
                            }
                        }
                    }

                    // Add identity
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_gradTens[i * SPATIAL_DIM + j] +=
                                p_kd[i * SPATIAL_DIM + j];
                        }
                    }

                    // Compute I = tr(F^T * F)
                    scalar Iconst = 0.0;
                    if (SPATIAL_DIM == 3)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            Iconst += p_gradTens[i] * p_gradTens[i] +
                                      p_gradTens[i + 3] * p_gradTens[i + 3] +
                                      p_gradTens[i + 6] * p_gradTens[i + 6];
                        }
                    }
                    else
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            Iconst += p_gradTens[i] * p_gradTens[i] +
                                      p_gradTens[i + 2] * p_gradTens[i + 2];
                        }
                    }

                    // Beta parameter
                    const scalar beta = 1.0 - 1.0 / (1.0 + Iconst);

                    // Jacobian J = det(F)
                    scalar J;
                    if (SPATIAL_DIM == 3)
                    {
                        J = p_gradTens[0] * (p_gradTens[4] * p_gradTens[8] -
                                             p_gradTens[7] * p_gradTens[5]) -
                            p_gradTens[1] * (p_gradTens[3] * p_gradTens[8] -
                                             p_gradTens[6] * p_gradTens[5]) +
                            p_gradTens[2] * (p_gradTens[3] * p_gradTens[7] -
                                             p_gradTens[6] * p_gradTens[4]);
                    }
                    else
                    {
                        J = p_gradTens[0] * p_gradTens[3] -
                            p_gradTens[2] * p_gradTens[1];
                    }

                    // Assemble stress contributions
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label icNdim = ic * SPATIAL_DIM;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const label indexL = il * SPATIAL_DIM + i;
                            const label indexR = ir * SPATIAL_DIM + i;

                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];

                            scalar lhs_riC_i = 0.0;

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    p_scs_areav[ip * SPATIAL_DIM + j];
                                const scalar dxj = p_D[ic * SPATIAL_DIM + j];

                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;

                                // First set of terms: mu*beta/J * F_ik * dN/dxk
                                // * Aj
                                for (label k = 0; k < SPATIAL_DIM; ++k)
                                {
                                    const scalar factor =
                                        p_gradTens[SPATIAL_DIM * i + k] / J;
                                    const scalar lhsfacDiff_k =
                                        -(muIp * beta) * factor *
                                        p_dndx[offSetDnDx + k] * axj;

                                    p_lhs[indexL * rhsSize + icNdim + j] +=
                                        lhsfacDiff_k;
                                    p_lhs[indexR * rhsSize + icNdim + j] -=
                                        lhsfacDiff_k;

                                    p_rhs[indexL] -= lhsfacDiff_k * dxj;
                                    p_rhs[indexR] += lhsfacDiff_k * dxj;
                                }

                                // Second set of terms: mu*beta/J * dN/dxj * Aj
                                const scalar lhsfacDiff_j =
                                    -((muIp * beta) / J) *
                                    p_dndx[offSetDnDx + j] * axj;
                                lhs_riC_i += lhsfacDiff_j;

                                // Third set of terms: mu*beta/J * dN/dxi * Aj
                                const scalar lhsfacDiff_i =
                                    -((muIp * beta) / J) *
                                    p_dndx[offSetDnDx + i] * axj;
                                p_lhs[indexL * rhsSize + icNdim + j] +=
                                    lhsfacDiff_i;
                                p_lhs[indexR * rhsSize + icNdim + j] -=
                                    lhsfacDiff_i;
                                p_rhs[indexL] -= lhsfacDiff_i * dxj;
                                p_rhs[indexR] += lhsfacDiff_i * dxj;
                            }

                            // Accumulated diagonal term
                            p_lhs[indexL * rhsSize + icNdim + i] += lhs_riC_i;
                            p_lhs[indexR * rhsSize + icNdim + i] -= lhs_riC_i;
                            p_rhs[indexL] -= lhs_riC_i * dxi;
                            p_rhs[indexR] += lhs_riC_i * dxi;
                        }
                    }

                    // Additional RHS terms (volumetric contribution)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label indexL = il * SPATIAL_DIM + i;
                        const label indexR = ir * SPATIAL_DIM + i;

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj =
                                p_scs_areav[ip * SPATIAL_DIM + j];
                            const scalar rhsTerm =
                                lambdaIp * p_kd[i * SPATIAL_DIM + j] * axj *
                                    (J - (1.0 + 0.75 * muIp / lambdaIp)) +
                                (muIp * beta / J) * p_kd[i * SPATIAL_DIM + j] *
                                    axj;

                            p_rhs[indexL] += rhsTerm;
                            p_rhs[indexR] -= rhsTerm;
                        }
                    }
                }
            }

            // Apply coefficients to system (stress terms only, mass in
            // NodeTerms)
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
    }
    else
    {
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const solidMechanicsOption solidMechOption =
        domain->solidMechanics_.option_;
    STK_ThrowRequireMsg(
        solidMechOption == solidMechanicsOption::linearElastic ||
            solidMechOption == solidMechanicsOption::neoHookean ||
            solidMechOption == solidMechanicsOption::modifiedMooneyRivlin,
        "Unsupported FEM solid mechanics option");
#if SPATIAL_DIM == 3
    STK_ThrowRequireMsg(
        solidMechOption != solidMechanicsOption::modifiedMooneyRivlin,
        "modified_mooney_rivlin is only wired for the 2D FEM path");
#endif

    const auto& DSTKFieldRef = phi_->stkFieldRef();
    const auto& ESTKFieldRef = model_->ERef().stkFieldRef();
    const auto& nuSTKFieldRef = model_->nuRef().stkFieldRef();
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    const stk::mesh::Selector selAllElements =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // Keep the FEM constitutive parameters aligned with the OpenAccel input
    // switch so the lambda term follows the same plane_stress selection as the
    // CVFEM solid path.
    const bool usePlaneStressForFEM = domain->solidMechanics_.planeStress_;

    for (const stk::mesh::Bucket* bucket : elementBuckets)
    {
        const stk::topology topology = bucket->topology();
        smesh::ElemType sfemElementType = smesh::INVALID;
        bool useOpenAccelQuad4 = false;

#if SPATIAL_DIM == 3
        if (topology == stk::topology::HEX_8)
            sfemElementType = smesh::HEX8;
        else if (topology == stk::topology::TET_4)
            sfemElementType = smesh::TET4;
        else if (topology == stk::topology::TET_10)
            sfemElementType = smesh::TET10;
#else
        if (topology == stk::topology::TRI_3_2D)
            sfemElementType = smesh::TRI3;
        else if (topology == stk::topology::QUAD_4_2D)
        {
            sfemElementType = smesh::QUAD4;
            useOpenAccelQuad4 = true;
        }
#endif

        STK_ThrowRequireMsg(
            sfemElementType != smesh::INVALID || useOpenAccelQuad4,
            "Unsupported FEM solid element topology: " << topology.name());

        const label nodesPerElement = topology.num_nodes();
        const label dofsPerElement = nodesPerElement * SPATIAL_DIM;
        const label lhsSize = dofsPerElement * dofsPerElement;

        lhs.resize(lhsSize);
        rhs.resize(dofsPerElement);
        scratchIds.resize(dofsPerElement);
        scratchVals.resize(dofsPerElement);
        connectedNodes.resize(nodesPerElement);

        std::vector<idx_t> localNodeIds(nodesPerElement);
        std::vector<idx_t*> elementConnectivity(nodesPerElement);
        std::vector<std::vector<geom_t>> coordinates(
            SPATIAL_DIM, std::vector<geom_t>(nodesPerElement));
        std::vector<std::vector<scalar>> hyperCoordinates(
            SPATIAL_DIM, std::vector<scalar>(nodesPerElement));
        std::vector<geom_t*> coordinatePointers(SPATIAL_DIM);
        std::vector<real_t> basis(dofsPerElement, 0.0);
        std::vector<real_t> column(dofsPerElement, 0.0);
        std::vector<scalar> displacement(dofsPerElement, 0.0);

        for (label node = 0; node < nodesPerElement; ++node)
        {
            localNodeIds[node] = node;
            elementConnectivity[node] = &localNodeIds[node];
        }
        for (label dim = 0; dim < SPATIAL_DIM; ++dim)
            coordinatePointers[dim] = coordinates[dim].data();

        for (stk::mesh::Entity elem : *bucket)
        {
            std::fill(lhs.begin(), lhs.end(), 0.0);
            std::fill(rhs.begin(), rhs.end(), 0.0);
            std::fill(displacement.begin(), displacement.end(), 0.0);

            scalar mu = 0.0;
            scalar lambda = 0.0;
            const stk::mesh::Entity* nodeRels = bulkData.begin_nodes(elem);
            STK_ThrowAssert(bulkData.num_nodes(elem) == nodesPerElement);

            for (label node = 0; node < nodesPerElement; ++node)
            {
                connectedNodes[node] = nodeRels[node];
                const scalar* xyz =
                    stk::mesh::field_data(coordinatesRef, nodeRels[node]);
                const scalar* D =
                    stk::mesh::field_data(DSTKFieldRef, nodeRels[node]);
                const scalar E =
                    *stk::mesh::field_data(ESTKFieldRef, nodeRels[node]);
                const scalar nu =
                    *stk::mesh::field_data(nuSTKFieldRef, nodeRels[node]);

                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                {
                    coordinates[dim][node] = xyz[dim];
                    hyperCoordinates[dim][node] = xyz[dim];
                    displacement[node * SPATIAL_DIM + dim] = D[dim];
                }

                mu += E / (2.0 * (1.0 + nu));
                lambda += usePlaneStressForFEM
                              ? nu * E / ((1.0 + nu) * (1.0 - nu))
                              : nu * E /
                                    ((1.0 + nu) * (1.0 - 2.0 * nu));
            }

            mu /= nodesPerElement;
            lambda /= nodesPerElement;

            if (solidMechOption == solidMechanicsOption::neoHookean)
            {
                assembleSfemNeoHookeanElement(sfemElementType,
                                              nodesPerElement,
                                              lambda,
                                              mu,
                                              hyperCoordinates,
                                              displacement,
                                              rhs,
                                              lhs);
            }
#if SPATIAL_DIM == 2
            else if (solidMechOption ==
                     solidMechanicsOption::modifiedMooneyRivlin)
            {
                // c1/c2/kappa are simple per-domain material constants (see
                // domain::material::mechanicalProperties_, Task 2), unlike
                // mu/lambda above which come from the E/nu node fields --
                // not needed for this branch.
                const auto& mechProps =
                    domain->materialRef().mechanicalProperties_;
                assembleSfemModifiedMooneyRivlinElement(sfemElementType,
                                                        nodesPerElement,
                                                        mechProps.c1_,
                                                        mechProps.c2_,
                                                        mechProps.kappa_,
                                                        coordinates,
                                                        elementConnectivity,
                                                        displacement,
                                                        rhs,
                                                        lhs);
            }
            else if (useOpenAccelQuad4)
            {
                assembleQuad4LinearElasticity(coordinates, mu, lambda, lhs);
            }
#endif
            else
            {
#if SPATIAL_DIM == 2
                // sfem provides the isoparametric matrix-free operator. Applying
                // it to each local basis vector recovers the element matrix while
                // OpenAccel retains ownership of the global block matrix.
#endif
                for (label col = 0; col < dofsPerElement; ++col)
                {
                    std::fill(basis.begin(), basis.end(), 0.0);
                    std::fill(column.begin(), column.end(), 0.0);
                    basis[col] = 1.0;

                    const int status = linear_elasticity_apply_aos(
                        sfemElementType,
                        1,
                        nodesPerElement,
                        elementConnectivity.data(),
                        coordinatePointers.data(),
                        mu,
                        lambda,
                        basis.data(),
                        column.data());
                    STK_ThrowRequireMsg(status == SFEM_SUCCESS,
                                        "sfem linear elasticity kernel failed");

                    for (label row = 0; row < dofsPerElement; ++row)
                        lhs[row * dofsPerElement + col] = column[row];
                }
            }

            static bool femDebugPrinted = false;
            const char* femDebug = std::getenv("OPENACCEL_FEM_DEBUG");
            if (!femDebugPrinted && femDebug && std::atoi(femDebug) > 0 &&
                messager::myProcNo() == 0)
            {
                femDebugPrinted = true;
                scalar maxEntry = 0.0;
                scalar maxAsymmetry = 0.0;
                scalar maxTranslationResidual = 0.0;
                scalar minDiagonal = std::numeric_limits<scalar>::max();
                scalar maxDiagonal = -std::numeric_limits<scalar>::max();
                scalar probeEnergy = 0.0;

                std::vector<scalar> probe(dofsPerElement);
                for (label row = 0; row < dofsPerElement; ++row)
                    probe[row] = scalar((row * 17 + 3) % 11) - 5.0;

                for (label row = 0; row < dofsPerElement; ++row)
                {
                    const scalar diagonal =
                        lhs[row * dofsPerElement + row];
                    minDiagonal = std::min(minDiagonal, diagonal);
                    maxDiagonal = std::max(maxDiagonal, diagonal);

                    scalar probeRow = 0.0;
                    for (label col = 0; col < dofsPerElement; ++col)
                    {
                        const scalar value =
                            lhs[row * dofsPerElement + col];
                        maxEntry = std::max(maxEntry, std::abs(value));
                        maxAsymmetry = std::max(
                            maxAsymmetry,
                            std::abs(value -
                                     lhs[col * dofsPerElement + row]));
                        probeRow += value * probe[col];
                    }
                    probeEnergy += probe[row] * probeRow;

                    for (label component = 0; component < SPATIAL_DIM;
                         ++component)
                    {
                        scalar translationResidual = 0.0;
                        for (label node = 0; node < nodesPerElement; ++node)
                        {
                            translationResidual +=
                                lhs[row * dofsPerElement +
                                    node * SPATIAL_DIM + component];
                        }
                        maxTranslationResidual =
                            std::max(maxTranslationResidual,
                                     std::abs(translationResidual));
                    }
                }

                std::cout << "\n[FEM DEBUG] First FEM element\n"
                          << "  topology: " << topology.name() << '\n'
                          << "  nodes: " << nodesPerElement << '\n'
                          << "  mu: " << mu << '\n'
                          << "  lambda: " << lambda << '\n';
                for (label node = 0; node < nodesPerElement; ++node)
                {
                    std::cout << "  x[" << node << "]:";
                    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                        std::cout << ' ' << coordinates[dim][node];
                    std::cout << '\n';
                }
                std::cout << "  max_abs_K: " << maxEntry << '\n'
                          << "  max_abs_K_minus_KT: " << maxAsymmetry << '\n'
                          << "  max_translation_residual: "
                          << maxTranslationResidual << '\n'
                          << "  diagonal_min/max: " << minDiagonal << " / "
                          << maxDiagonal << '\n'
                          << "  probe_energy_vKv: " << probeEnergy << '\n'
                          << std::endl;
            }

            // OpenAccel solves for an increment, not for the total field:
            // K * deltaD = f_ext - K * D. Boundary assembly adds f_ext.
            if (solidMechOption == solidMechanicsOption::linearElastic)
            {
                for (label row = 0; row < dofsPerElement; ++row)
                {
                    for (label col = 0; col < dofsPerElement; ++col)
                    {
                        rhs[row] -= lhs[row * dofsPerElement + col] *
                                    displacement[col];
                    }
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
    }
}

} // namespace accel
