// File       : navierStokesEquation.cpp
// Created    : Wed Mar 13 2024 13:41:00 (+0100)
// Author     : Fabian Wermelinger
// Description: Coupled momentum equation implementation details
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "navierStokesEquation.h"
#include "realm.h"
#include "scaling.h"
#include "interface.h"
#include "ipInfo.h"
#include "simulation.h"

// std
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>

namespace accel
{

navierStokesEquation::navierStokesEquation(realm* realm, flowModel* model)
    : equation("Coupled Navier-Stokes", true), model_(model),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(model))
{
    // Collect equation names
    std::vector<std::string> eqNames;

#if SPATIAL_DIM == 2
    eqNames.push_back("U-Mom");
    eqNames.push_back("V-Mom");
#elif SPATIAL_DIM == 3
    eqNames.push_back("U-Mom");
    eqNames.push_back("V-Mom");
    eqNames.push_back("W-Mom");
#endif

    // set
    this->setEquationName(eqNames);

    // set relaxation for velocity
    model_->URef().setURF(model_->controlsRef()
                              .solverRef()
                              .solverControl_.basicSettings_.convergenceControl_
                              .relaxationParameters_.velocityRelaxationFactor_);
}

void navierStokesEquation::checkDomain(const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool navierStokesEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void navierStokesEquation::setup()
{
    // setup DU coefficients and body forces
    FOREACH_DOMAIN_PTR(assembler_->setupDUCoefficients);
    FOREACH_DOMAIN(model_->setupBodyForces);

    // force creation of beta field for U if not yet
    if (model_->URef().blendingFactorPtr() == nullptr)
    {
        model_->URef().setupBlendingFactorField(
            /*enable only a dummy field*/ true);
    }

    if (model_->rhoRef().blendingFactorPtr() == nullptr)
    {
        model_->rhoRef().setupBlendingFactorField(
            /*enable only a dummy field*/ true);
    }

    // setup required fields
    FOREACH_DOMAIN(model_->setupDensity);
    FOREACH_DOMAIN(model_->setupVelocity);
    FOREACH_DOMAIN(model_->setupMassFlowRate);
    FOREACH_DOMAIN(model_->setupMomentumFlowRate);
    FOREACH_DOMAIN(model_->setupDynamicViscosity);

    // setup assembler
    assembler_->setup(&model_->URef(), null, domainVector_, nullptr);

    // setup linear solver
    linearSystem::setupSolver(this->name(),
                              model_->meshRef(),
                              this->domainZones_(),
                              this->fallbackName());

    equation::isCreated_ = true;
}

void navierStokesEquation::initialize()
{
    // raw initialization of velocity
    FOREACH_DOMAIN(model_->initializeVelocity);

    // update velocity gradient
    FOREACH_DOMAIN(model_->updateVelocityGradientField);

    // update high-res fields
    FOREACH_DOMAIN(model_->updateVelocityBlendingFactorField);

    // update velocity scale
    model_->URef().updateScale();
}

void navierStokesEquation::postInitialize()
{
    // initialization of properties
    FOREACH_DOMAIN(model_->initializeDensity);
    FOREACH_DOMAIN_IF(model_->updateDensityGradientField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN_IF(model_->updateDensityBlendingFactorField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN(model_->initializeDynamicViscosity);

    // post initialization
    model_->rhoRef().updateScale();
    model_->muRef().updateScale();

    // raw initialization of mass flux: can be safely done now after
    // initialization of velocity and density at all active domains
    FOREACH_DOMAIN(model_->initializeMassFlowRate);

    // transform mDot to relative: availability of a moving frame
    // will be checked inside. Then, update divergence fields
    FOREACH_DOMAIN(model_->transformMassFlowRateToRelative);

    // in case of transient, old density must be updated before div update
    // (not after a restart that restored the old density state: that would
    // overwrite it, see nodeField::initialize)
    if (model_->controlsRef().isTransient() &&
        !model_->rhoRef().restoredWithPrevTimeState())
    {
        FOREACH_DOMAIN(model_->updateDensityPrevTimeField);
    }

    // now update div field
    FOREACH_DOMAIN(model_->updateMassDivergenceField);

    equation::isInitialized_ = true;
}

void navierStokesEquation::preSolve()
{
    // velocity update
    FOREACH_DOMAIN(model_->updateVelocityPrevIterField);
    FOREACH_DOMAIN(model_->updateVelocity);

    // properties update
    FOREACH_DOMAIN_IF(model_->updateDensity, domain->isMaterialCompressible());
    FOREACH_DOMAIN(model_->updateDynamicViscosity);
    FOREACH_DOMAIN(model_->updateEffectiveDynamicViscosity);

    // update density-related fields
    FOREACH_DOMAIN_IF(model_->updateDensityGradientField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN_IF(model_->updateDensityBlendingFactorField,
                      domain->isMaterialCompressible());

    // update body forces
    FOREACH_DOMAIN(model_->computeBodyForces);
    FOREACH_DOMAIN(model_->redistributeBodyForces);
}

void navierStokesEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assemble
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get(), {});

    if (!model_->controlsRef()
             .solverRef()
             .solverControl_.expertParameters_.disableMomentumPredictor_)
    {
        // solve linear system
        if (ctx->getGraph()->isGraphMember())
        {
            linearSystem::solve();
        }
        messager::barrier();

        // correction
        auto& U = model_->URef();
        for (const auto& domain : domainVector_)
        {
            scalar effectiveRelaxationFactor = 1.0;

            // find the ramp value (explicit relaxation) for compressible
            // domains
            if (domain->isMaterialCompressible())
            {
                scalar iter = model_->controlsRef().iter;
                label rampIter = 20;
                if (iter <= rampIter)
                {
                    effectiveRelaxationFactor *=
                        std::max(scalar(iter) / scalar(rampIter), 0.1);
                }
            }

            this->template correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
                domain.get(),
                ctx->coeffs().get(),
                stk::topology::NODE_RANK,
                U.stkFieldRef(),
                effectiveRelaxationFactor);

            // synchronize
            U.synchronizeGhostedEntities(domain->index());
        }

        // post correction

        // velocity gradient must be only updated after corrector step ...

        // update scale
        U.updateScale();
    }
}

// ---------------------------------------------------------------------------
// DIAGNOSTIC (read-only): discrete x-momentum and kinetic-energy budgets of
// the converged step that has just been completed.
//
// Active only when the environment variable OPENACCEL_BUDGET_FILE names an
// output file; otherwise it is never called and the solver is unchanged.
// It reads the solver's own fields and recomputes, ip by ip, exactly the terms
// the momentum assembler builds (navierStokesAssemblerElemTerms.cpp,
// ...NodeTerms.cpp second-order unsteady, ...ElemBoundaryConditions.cpp no-slip
// wall), using the stored mass flow rates mDot. Nothing is written to any
// field. Called at the start of preTimeStep, before the previous-time states
// are shifted: U = U^{n}, U_old = U^{n-1}, U_oldold = U^{n-2}, and the momentum
// source still holds the force applied during step n.
//
// Sign convention: every term is reported as its contribution to the nodal
// right-hand side (net force on the control volume); at convergence the terms
// of each node, and hence of each wall-parallel level, sum to zero.
// ---------------------------------------------------------------------------
namespace
{

void diagnoseDiscreteBudget(flowModel* model,
                            const std::shared_ptr<domain> domain,
                            const std::string& path)
{
    const label nd = SPATIAL_DIM;
    stk::mesh::MetaData& metaData = model->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = model->meshRef().bulkDataRef();

    const auto& Uf = model->URef().stkFieldRef();
    const auto& UOld = model->URef().prevTimeRef().stkFieldRef();
    const auto& UOldOld = model->URef().prevTimeRef().prevTimeRef().stkFieldRef();
    const auto& rhoF = model->rhoRef().stkFieldRef();
    const auto& muEffF = model->muEffRef().stkFieldRef();
    const auto& muF = model->muRef().stkFieldRef();
    const auto& gradUF = model->URef().gradRef().stkFieldRef();
    const auto& betaF = model->URef().blendingFactorRef().stkFieldRef();
    const auto& gradPF = model->pRef().gradRef().stkFieldRef();
    const auto& mDotF = model->mDotRef().stkFieldRef();
    const auto& divF = model->mDotRef().divRef().stkFieldRef();
    const auto& coordF = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);
    const auto& volF = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);
    const auto* FF =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto& areaF = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);
    const auto& coeffF = model->uWallCoeffsRef().stkFieldRef();
    const auto& sideUF = model->URef().sideFieldRef().stkFieldRef();

    const bool isUShifted = model->URef().isShifted();
    const bool isUGradientShifted = model->URef().isGradientShifted();
    const auto& controls = model->controlsRef();
    const scalar dt = controls.getTimestep();
    const auto c = BDF2::coeff(dt, controls.getTimestep(-1));
    const std::vector<scalar> msrc = domain->generalMomentumSourceRef().value_;

    scalar centralFrac = 0.0;
    bool centralEdge = false;      // same choice as the assembler
    {
        const YAML::Node solverControl =
            model->realmRef().simulationRef().getYAMLSolverControlNode();
        const YAML::Node expert = solverControl["expert_parameters"];
        if (expert && expert["momentum_advection_central_fraction"])
        {
            centralFrac = expert["momentum_advection_central_fraction"]
                              .template as<scalar>();
        }
        if (expert && expert["momentum_advection_central_interpolation"])
        {
            centralEdge = expert["momentum_advection_central_interpolation"]
                              .template as<std::string>() == "edge";
        }
    }

    const stk::mesh::PartVector& parts = domain->zonePtr()->interiorParts();
    const stk::mesh::Selector selOwned =
        metaData.locally_owned_part() & stk::mesh::selectUnion(parts);

    // global wall-parallel node levels (y), built once
    static std::vector<double> levels;
    if (levels.empty())
    {
        std::vector<double> loc;
        for (const stk::mesh::Bucket* b :
             bulkData.get_buckets(stk::topology::NODE_RANK, selOwned))
        {
            const scalar* x = stk::mesh::field_data(coordF, *b);
            for (size_t k = 0; k < b->size(); ++k)
            {
                loc.push_back(std::round(x[nd * k + 1] * 1e8) / 1e8);
            }
        }
        std::sort(loc.begin(), loc.end());
        loc.erase(std::unique(loc.begin(), loc.end()), loc.end());
        int n = static_cast<int>(loc.size());
        const int np = messager::nProcs();
        std::vector<int> counts(np), displs(np, 0);
        MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, messager::comm());
        for (int r = 1; r < np; ++r)
        {
            displs[r] = displs[r - 1] + counts[r - 1];
        }
        std::vector<double> all(displs[np - 1] + counts[np - 1]);
        MPI_Allgatherv(loc.data(), n, MPI_DOUBLE, all.data(), counts.data(),
                       displs.data(), MPI_DOUBLE, messager::comm());
        std::sort(all.begin(), all.end());
        all.erase(std::unique(all.begin(), all.end()), all.end());
        levels = all;
    }
    const label NL = static_cast<label>(levels.size());
    const label NP = NL - 1;
    auto levelOf = [&](scalar y)
    {
        const double yr = std::round(y * 1e8) / 1e8;
        auto it = std::lower_bound(levels.begin(), levels.end(), yr - 1e-7);
        return static_cast<label>(it - levels.begin());
    };

    // accumulators
    enum LTerm { L_TIME, L_PRESS, L_SRC, L_DIV, L_WALL, L_ADV, L_MOL, L_SGS, NLT };
    enum PTerm { P_ADV, P_MOL, P_SGS, P_MDOT, P_IPPROD, P_NODEPROD, P_AREA, NPT };
    enum GTerm { G_ADV, G_MOL, G_SGS, G_PRESS, G_SRC, G_DIV, G_WALL, G_TIME,
                 G_KE, G_KEOLD, G_MOM, G_VOL, NGT };
    std::vector<double> acc(NL * NLT + NP * NPT + NGT, 0.0);
    double* Lv = acc.data();
    double* Pv = acc.data() + NL * NLT;
    double* Gv = acc.data() + NL * NLT + NP * NPT;

    // ---- node terms (owned nodes) ----
    for (const stk::mesh::Bucket* b :
         bulkData.get_buckets(stk::topology::NODE_RANK, selOwned))
    {
        const scalar* U = stk::mesh::field_data(Uf, *b);
        const scalar* Uo = stk::mesh::field_data(UOld, *b);
        const scalar* Uoo = stk::mesh::field_data(UOldOld, *b);
        const scalar* rho = stk::mesh::field_data(rhoF, *b);
        const scalar* vol = stk::mesh::field_data(volF, *b);
        const scalar* gp = stk::mesh::field_data(gradPF, *b);
        const scalar* dv = stk::mesh::field_data(divF, *b);
        const scalar* x = stk::mesh::field_data(coordF, *b);
        const scalar* F = FF ? stk::mesh::field_data(*FF, *b) : nullptr;
        for (size_t k = 0; k < b->size(); ++k)
        {
            const label l = levelOf(x[nd * k + 1]);
            const scalar V = vol[k];
            double ke = 0, keo = 0;
            for (label i = 0; i < nd; ++i)
            {
                const scalar ui = U[nd * k + i];
                const scalar tTerm = -(c[0] * U[nd * k + i] + c[1] * Uo[nd * k + i] +
                                       c[2] * Uoo[nd * k + i]) *
                                     rho[k] * V / dt;
                const scalar pTerm = -gp[nd * k + i] * V;
                const scalar sTerm = ((F ? F[nd * k + i] : 0.0) + msrc[i]) * V;
                const scalar dTerm = dv[k] * ui;
                if (i == 0)
                {
                    Lv[l * NLT + L_TIME] += tTerm;
                    Lv[l * NLT + L_PRESS] += pTerm;
                    Lv[l * NLT + L_SRC] += sTerm;
                    Lv[l * NLT + L_DIV] += dTerm;
                    Gv[G_MOM] += rho[k] * V * ui;
                }
                Gv[G_TIME] += ui * tTerm;
                Gv[G_PRESS] += ui * pTerm;
                Gv[G_SRC] += ui * sTerm;
                Gv[G_DIV] += ui * dTerm;
                ke += ui * ui;
                keo += Uo[nd * k + i] * Uo[nd * k + i];
            }
            Gv[G_KE] += 0.5 * rho[k] * V * ke;
            Gv[G_KEOLD] += 0.5 * rho[k] * V * keo;
            Gv[G_VOL] += V;
        }
    }

    // ---- element (SCS) terms (owned elements) ----
    std::vector<scalar> wU, wBeta, wX, wGrad, wRho, wMuEff, wMu, area, dndx, deriv,
        detj, sfU, sfX;
    for (const stk::mesh::Bucket* bp :
         bulkData.get_buckets(stk::topology::ELEMENT_RANK, selOwned))
    {
        const stk::mesh::Bucket& b = *bp;
        MasterElement* me =
            MasterElementRepo::get_surface_master_element(b.topology());
        const label npe = me->nodesPerElement_;
        const label nip = me->numIntPoints_;
        const label* lr = me->adjacentNodes();
        wU.resize(npe * nd); wBeta.resize(npe * nd); wX.resize(npe * nd);
        wGrad.resize(npe * nd * nd); wRho.resize(npe); wMuEff.resize(npe);
        wMu.resize(npe); area.resize(nip * nd); dndx.resize(nd * nip * npe);
        deriv.resize(nd * nip * npe); detj.resize(nip); sfU.resize(nip * npe);
        sfX.resize(nip * npe);
        if (isUShifted)
        {
            me->shifted_shape_fcn(sfU.data());
        }
        else
        {
            me->shape_fcn(sfU.data());
        }
        me->shape_fcn(sfX.data());
        for (size_t e = 0; e < b.size(); ++e)
        {
            const scalar* mDot = stk::mesh::field_data(mDotF, b, e);
            stk::mesh::Entity const* nodes = b.begin_nodes(e);
            for (label n = 0; n < npe; ++n)
            {
                const scalar* U = stk::mesh::field_data(Uf, nodes[n]);
                const scalar* be = stk::mesh::field_data(betaF, nodes[n]);
                const scalar* x = stk::mesh::field_data(coordF, nodes[n]);
                const scalar* g = stk::mesh::field_data(gradUF, nodes[n]);
                wRho[n] = *stk::mesh::field_data(rhoF, nodes[n]);
                wMuEff[n] = *stk::mesh::field_data(muEffF, nodes[n]);
                wMu[n] = *stk::mesh::field_data(muF, nodes[n]);
                for (label i = 0; i < nd; ++i)
                {
                    wU[n * nd + i] = U[i];
                    wBeta[n * nd + i] = be[i];
                    wX[n * nd + i] = x[i];
                    for (label j = 0; j < nd; ++j)
                    {
                        wGrad[n * nd * nd + i * nd + j] = g[i * nd + j];
                    }
                }
            }
            scalar err = 0.0;
            me->determinant(1, wX.data(), area.data(), &err);
            if (isUGradientShifted)
            {
                me->shifted_grad_op(1, wX.data(), dndx.data(), deriv.data(),
                                    detj.data(), &err);
            }
            else
            {
                me->grad_op(1, wX.data(), dndx.data(), deriv.data(), detj.data(),
                            &err);
            }
            for (label ip = 0; ip < nip; ++ip)
            {
                const label il = lr[2 * ip], ir = lr[2 * ip + 1];
                const scalar md = mDot[ip];
                scalar Uip[SPATIAL_DIM] = {0}, Xip[SPATIAL_DIM] = {0};
                scalar muEffIp = 0, muIp = 0, rhoIp = 0;
                for (label n = 0; n < npe; ++n)
                {
                    const scalar r = sfU[ip * npe + n];
                    muEffIp += r * wMuEff[n];
                    muIp += r * wMu[n];
                    rhoIp += r * wRho[n];
                    for (label i = 0; i < nd; ++i)
                    {
                        Uip[i] += r * wU[n * nd + i];
                        Xip[i] += sfX[ip * npe + n] * wX[n * nd + i];
                    }
                }
                // fluxes leaving L towards R (rhs_L -= F, rhs_R += F)
                scalar Fadv[SPATIAL_DIM], Fmol[SPATIAL_DIM], Fsgs[SPATIAL_DIM];
                for (label i = 0; i < nd; ++i)
                {
                    const label up = md > 0 ? il : ir;
                    scalar dcorr = 0.0;
                    for (label j = 0; j < nd; ++j)
                    {
                        dcorr += wBeta[up * nd + i] * (Xip[j] - wX[up * nd + j]) *
                                 wGrad[up * nd * nd + i * nd + j];
                    }
                    const scalar uUp = wU[up * nd + i];
                    const scalar uC = centralEdge
                                          ? 0.5 * (wU[il * nd + i] + wU[ir * nd + i])
                                          : Uip[i];
                    Fadv[i] = centralFrac > 0.0
                                  ? md * ((1.0 - centralFrac) * (uUp + dcorr) +
                                          centralFrac * uC)
                                  : md * (uUp + dcorr);
                    scalar gsum = 0.0;
                    for (label n = 0; n < npe; ++n)
                    {
                        const label off = nd * npe * ip + n * nd;
                        for (label j = 0; j < nd; ++j)
                        {
                            gsum += -(dndx[off + j] * wU[n * nd + i] +
                                      dndx[off + i] * wU[n * nd + j]) *
                                    area[ip * nd + j];
                        }
                    }
                    Fmol[i] = muIp * gsum;
                    Fsgs[i] = (muEffIp - muIp) * gsum;
                }
                const label lL = levelOf(wX[il * nd + 1]);
                const label lR = levelOf(wX[ir * nd + 1]);
                Lv[lL * NLT + L_ADV] -= Fadv[0];
                Lv[lR * NLT + L_ADV] += Fadv[0];
                Lv[lL * NLT + L_MOL] -= Fmol[0];
                Lv[lR * NLT + L_MOL] += Fmol[0];
                Lv[lL * NLT + L_SGS] -= Fsgs[0];
                Lv[lR * NLT + L_SGS] += Fsgs[0];
                for (label i = 0; i < nd; ++i)
                {
                    const scalar du = wU[il * nd + i] - wU[ir * nd + i];
                    Gv[G_ADV] -= Fadv[i] * du;
                    Gv[G_MOL] -= Fmol[i] * du;
                    Gv[G_SGS] -= Fsgs[i] * du;
                }
                if (lL != lR)
                {
                    // y-plane between levels: express as flux in +y
                    const scalar s = lR > lL ? 1.0 : -1.0;
                    const label p = std::min(lL, lR);
                    const scalar Ay = s * area[ip * nd + 1];
                    scalar uDotA = 0.0;
                    for (label j = 0; j < nd; ++j)
                    {
                        uDotA += Uip[j] * area[ip * nd + j];
                    }
                    Pv[p * NPT + P_ADV] += s * Fadv[0];
                    Pv[p * NPT + P_MOL] += s * Fmol[0];
                    Pv[p * NPT + P_SGS] += s * Fsgs[0];
                    Pv[p * NPT + P_MDOT] += s * md;
                    Pv[p * NPT + P_IPPROD] += s * rhoIp * uDotA * Uip[0];
                    Pv[p * NPT + P_NODEPROD] +=
                        Ay * 0.5 * (wRho[il] * wU[il * nd + 1] * wU[il * nd] +
                                    wRho[ir] * wU[ir * nd + 1] * wU[ir * nd]);
                    Pv[p * NPT + P_AREA] += Ay;
                }
            }
        }
    }

    // ---- no-slip wall flux (owned sides) ----
    for (label iB = 0; iB < domain->zonePtr()->nBoundaries(); iB++)
    {
        if (domain->zonePtr()->boundaryRef(iB).type() != boundaryPhysicalType::wall ||
            model->URef().boundaryConditionRef(domain->index(), iB).type() !=
                boundaryConditionType::noSlip)
        {
            continue;
        }
        const stk::mesh::Selector selSides =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(domain->zonePtr()->boundaryRef(iB).parts());
        for (const stk::mesh::Bucket* bp :
             bulkData.get_buckets(metaData.side_rank(), selSides))
        {
            const stk::mesh::Bucket& b = *bp;
            MasterElement* meFC =
                MasterElementRepo::get_surface_master_element(b.topology());
            const label nps = meFC->nodesPerElement_;
            const label nbip = meFC->numIntPoints_;
            const label* ipNode = meFC->ipNodeMap();
            std::vector<scalar> sf(nbip * nps);
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(sf.data());
            }
            else
            {
                meFC->shape_fcn(sf.data());
            }
            for (size_t s = 0; s < b.size(); ++s)
            {
                stk::mesh::Entity side = b[s];
                stk::mesh::Entity const* nodes = bulkData.begin_nodes(side);
                const scalar* av = stk::mesh::field_data(areaF, side);
                const scalar* cf = stk::mesh::field_data(coeffF, side);
                const scalar* ubc = stk::mesh::field_data(sideUF, side);
                for (label ip = 0; ip < nbip; ++ip)
                {
                    scalar aMag = 0, nx[SPATIAL_DIM], ub[SPATIAL_DIM] = {0};
                    for (label j = 0; j < nd; ++j)
                    {
                        aMag += av[ip * nd + j] * av[ip * nd + j];
                    }
                    aMag = std::sqrt(aMag);
                    for (label j = 0; j < nd; ++j)
                    {
                        nx[j] = av[ip * nd + j] / aMag;
                    }
                    for (label n = 0; n < nps; ++n)
                    {
                        const scalar* U = stk::mesh::field_data(Uf, nodes[n]);
                        for (label j = 0; j < nd; ++j)
                        {
                            ub[j] += sf[ip * nps + n] * U[j];
                        }
                    }
                    stk::mesh::Entity nn = nodes[ipNode[ip]];
                    const scalar* Un = stk::mesh::field_data(Uf, nn);
                    const label l = levelOf(stk::mesh::field_data(coordF, nn)[1]);
                    for (label i = 0; i < nd; ++i)
                    {
                        scalar uT = 0, bT = 0;
                        for (label j = 0; j < nd; ++j)
                        {
                            const scalar w = (i == j ? 1.0 : 0.0) - nx[i] * nx[j];
                            uT += w * ub[j];
                            bT += w * ubc[ip * nd + j];
                        }
                        const scalar wTerm = -cf[ip] * (uT - bT);
                        if (i == 0)
                        {
                            Lv[l * NLT + L_WALL] += wTerm;
                        }
                        Gv[G_WALL] += Un[i] * wTerm;
                    }
                }
            }
        }
    }

    messager::sumReduce(acc);
    if (messager::master())
    {
        static bool first = true;
        std::ofstream out(path, std::ios::app);
        out << std::setprecision(10);
        if (first)
        {
            out << "# discrete budget; x-momentum per level: time press src div "
                   "wall adv mol sgs | per plane (between level p and p+1, "
                   "flux in +y): adv mol sgs mdot ipprod nodeprod area | global "
                   "u.rhs: adv mol sgs press src div wall time KE KEold Mx vol"
                << "\n# gamma " << centralFrac << " edge " << (centralEdge ? 1 : 0)
                << " levels " << NL;
            for (const double y : levels)
            {
                out << " " << y;
            }
            out << "\n";
            first = false;
        }
        out << std::scientific << std::setprecision(10) << controls.time << " "
            << dt;
        for (const double v : acc)
        {
            out << " " << v;
        }
        out << "\n";
    }
}

} // namespace

void navierStokesEquation::preTimeStep()
{
    // read-only discrete budget diagnostic, off unless requested (see above)
    if (const char* budgetPath = std::getenv("OPENACCEL_BUDGET_FILE"))
    {
        for (const auto& dom : this->domainVector_)
        {
            diagnoseDiscreteBudget(model_, dom, budgetPath);
        }
    }
    FOREACH_DOMAIN(model_->updateVelocityPrevTimeField);
    FOREACH_DOMAIN(model_->updateDensityPrevTimeField);
    FOREACH_DOMAIN(updateBulkVelocityForce_);
}

// Bulk-velocity controller (optional, per domain). Called once per time step,
// before the first outer iteration, with U holding the converged solution of
// the previous step. The volume-averaged velocity along the flow direction,
//
//   Ub = sum_owned(V_i U_i . d) / sum_owned(V_i),   V_i = dual nodal volume,
//
// is compared with the target and the uniform body force F (per unit volume)
// is updated by a PI law with closed-loop response time T:
//
//   e  = Ub_target - Ub
//   FI = FI + e dt / T^2        (integral part: converges to the mean drag)
//   F  = FI + 2 e / T           (proportional part: damps the response)
//
// For dUb/dt = F - drag this gives Ub'' + (2/T) Ub' + Ub/T^2 = ..., i.e. a
// critically damped return to the target over a few T. F is held fixed over
// all outer iterations of the step (it enters the BDF2 momentum equations
// like the constant source did) and replaces the input source value, which
// only initializes FI (unless FI is restored from a restart database).
void navierStokesEquation::updateBulkVelocityForce_(
    const std::shared_ptr<domain> domain)
{
    auto& source = domain->generalMomentumSourceRef();
    auto& ctrl = source.bulkVelocityControl_;
    if (!ctrl.active_)
    {
        return;
    }

    auto& controls = model_->controlsRef();
    stk::mesh::MetaData& metaData = model_->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = model_->meshRef().bulkDataRef();

    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

    // volume-weighted sums over locally owned nodes (each node counted once)
    std::vector<double> sums(2, 0.0);
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    for (const stk::mesh::Bucket* bucket :
         bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes))
    {
        const scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, *bucket);
        const scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, *bucket);
        for (size_t iNode = 0; iNode < bucket->size(); ++iNode)
        {
            scalar Ud = 0.0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                Ud += Ub[SPATIAL_DIM * iNode + i] * ctrl.direction_[i];
            }
            sums[0] += volb[iNode] * Ud;
            sums[1] += volb[iNode];
        }
    }
    if (messager::parallel())
    {
        messager::sumReduce(sums);
    }
    const scalar bulkVelocity = sums[0] / sums[1];

    if (!ctrl.initialized_)
    {
        const bool restored = controls.bulkForceIntegralRestored();
        scalar initial = 0.0;
        for (label i = 0; i < SPATIAL_DIM; i++)
        {
            initial += source.value_[i] * ctrl.direction_[i];
        }
        ctrl.integralForce_ = restored ? controls.bulkForceIntegral() : initial;
        ctrl.totalVolume_ = sums[1];
        ctrl.initialized_ = true;

        if (messager::master())
        {
            std::cout << "Bulk velocity control (domain " << domain->name()
                      << "): target Ub = " << std::scientific
                      << std::setprecision(10) << ctrl.targetBulkVelocity_
                      << ", response time T = " << ctrl.responseTime_
                      << ", control volume = " << sums[1]
                      << ", integral force initialized to "
                      << ctrl.integralForce_
                      << (restored ? " (restored from restart database)"
                                   : " (from input source value)")
                      << std::defaultfloat << std::endl;
        }
    }

    const scalar dt = controls.getTimestep();
    const scalar T = ctrl.responseTime_;
    const scalar error = ctrl.targetBulkVelocity_ - bulkVelocity;
    ctrl.integralForce_ += error * dt / (T * T);
    const scalar force = ctrl.integralForce_ + 2.0 * error / T;

    for (label i = 0; i < SPATIAL_DIM; i++)
    {
        source.value_[i] = force * ctrl.direction_[i];
    }
    ctrl.force_ = force;
    ctrl.bulkVelocity_ = bulkVelocity;
    controls.setBulkForceIntegral(ctrl.integralForce_);

    if (messager::master())
    {
        std::cout << "Bulk velocity control: Ub = " << std::scientific
                  << std::setprecision(10) << bulkVelocity
                  << ", error = " << error << ", force = " << force
                  << " (integral " << ctrl.integralForce_ << ")"
                  << std::defaultfloat << std::endl;

        const fs::path file =
            controls.getPostProcessingDirectory() / "bulk_velocity_control";
        const bool newFile = !fs::exists(file);
        std::ofstream out(file, std::ios::app);
        if (newFile)
        {
            out << "time\tUb_start_of_step\terror\tforce\tintegral_force\n";
        }
        out << std::scientific << std::setprecision(12) << controls.time
            << "\t" << bulkVelocity << "\t" << error << "\t" << force << "\t"
            << ctrl.integralForce_ << "\n";
    }
}

void navierStokesEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    const scalar U_scale_inv = 1.0 / (model_->URef().scale() + ::accel::SMALL);
    for (int i = 0; i < SPATIAL_DIM; i++)
    {
        residual_scales[i] = U_scale_inv;
    }
}

void navierStokesEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->URef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->URef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->URef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->URef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
