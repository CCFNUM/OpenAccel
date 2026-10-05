// File       : massTransferTests.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Unit tests for the mass-transfer models (parameter validation
//              and the Rayleigh-Plesset cavitation kernel). No mesh or MPI
//              is needed; the configuration is the fluidPairModel::massTransfer
//              struct of domain.h.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "rayleighPlessetCavitationMassTransfer.h"
#include "cavitationModel.h"

// std
#include <cmath>
#include <functional>
#include <iostream>

using namespace accel;

namespace
{

int nFailed = 0;
int nChecks = 0;

void check(bool ok, const std::string& what)
{
    ++nChecks;
    if (!ok)
    {
        ++nFailed;
        std::cerr << "FAIL: " << what << std::endl;
    }
}

void checkNear(double a, double b, double relTol, const std::string& what)
{
    const double scale = std::max(std::abs(a), std::abs(b));
    check(scale == 0.0 || std::abs(a - b) <= relTol * scale,
          what + " (" + std::to_string(a) + " vs " + std::to_string(b) + ")");
}

// run f; check that it throws a runtime_error whose message contains `needle`
void checkThrows(const std::function<void()>& f,
                 const std::string& needle,
                 const std::string& what)
{
    ++nChecks;
    try
    {
        f();
    }
    catch (const std::runtime_error& e)
    {
        if (std::string(e.what()).find(needle) == std::string::npos)
        {
            ++nFailed;
            std::cerr << "FAIL: " << what << ": message `" << e.what()
                      << "` lacks `" << needle << "`" << std::endl;
        }
        return;
    }
    ++nFailed;
    std::cerr << "FAIL: " << what << ": no exception" << std::endl;
}

const std::string path = "fluid_pair_models[0].mass_transfer";

// reference configuration (CFX defaults, water at ~25 C)
fluidPairModel::massTransfer refConfig()
{
    fluidPairModel::massTransfer cfg;
    cfg.option_ = massTransferModelOption::cavitation;
    cfg.cavitationModel_ = cavitationModelOption::rayleighPlesset;
    cfg.liquidPhase_ = "water";
    cfg.vaporPhase_ = "water_vapor";
    cfg.saturationPressure_ = 3169.0;
    return cfg;
}

// analytic test values from the specification
constexpr double pv = 3169.0, rhoL = 997.0, rhoV = 0.023, Rnuc = 1e-6,
                 rnuc = 5e-4, Fvap = 50.0, Fcond = 0.01;

double analyticVap(double p, double av)
{
    return Fvap * 3.0 * rnuc * (1.0 - av) * rhoV / Rnuc *
           std::sqrt(2.0 / 3.0 * (pv - p) / rhoL);
}

double analyticCond(double p, double av)
{
    return -Fcond * 3.0 * av * rhoV / Rnuc *
           std::sqrt(2.0 / 3.0 * (p - pv) / rhoL);
}

void testDefaults()
{
    const fluidPairModel::massTransfer cfg = refConfig();
    check(cfg.nucleationSiteVolumeFraction_ == 5.0e-4, "default r_nuc");
    check(cfg.nucleationSiteRadius_ == 1.0e-6, "default R_nuc");
    check(cfg.vaporizationCoefficient_ == 50.0, "default F_vap");
    check(cfg.condensationCoefficient_ == 0.01, "default F_cond");
    check(cfg.underRelaxation_ == 0.25, "default omega");
    check(!cfg.pressureClippingForRate_, "clipping default false");
    check(cfg.includeContinuitySource_, "continuity source default true");

    fluidPairModel::massTransfer none;
    check(createMassTransferModel(none) == nullptr,
          "factory returns null for none");
    const auto model = createMassTransferModel(cfg);
    check(model && model->name() == "RayleighPlessetCavitationMassTransfer",
          "factory creates RP model");
    check(dynamic_cast<cavitationModel*>(model.get()) != nullptr,
          "RP model is a cavitationModel");
}

void testValidation()
{
    auto mutated = [&](const std::function<void(fluidPairModel::massTransfer&)>&
                           f) {
        auto cfg = refConfig();
        f(cfg);
        return cfg;
    };
    auto validate = [&](const fluidPairModel::massTransfer& cfg) {
        cavitationModel::validate(cfg, path);
    };

    validate(refConfig()); // valid: must not throw

    checkThrows([&] {
        validate(mutated([](auto& c) { c.saturationPressure_ = -1; }));
    }, path + ".saturation_pressure", "negative saturation pressure");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.saturationPressure_ = 0; }));
    }, path + ".saturation_pressure", "zero saturation pressure");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.nucleationSiteRadius_ = 0; }));
    }, path + ".nucleation_site_radius", "zero radius");
    checkThrows([&] {
        validate(
            mutated([](auto& c) { c.nucleationSiteVolumeFraction_ = 1.5; }));
    }, path + ".nucleation_site_volume_fraction", "r_nuc > 1");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.vaporizationCoefficient_ = -2; }));
    }, path + ".vaporization_coefficient", "negative F_vap");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.condensationCoefficient_ = -2; }));
    }, path + ".condensation_coefficient", "negative F_cond");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.underRelaxation_ = 0; }));
    }, path + ".under_relaxation", "omega = 0");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.underRelaxation_ = 1.5; }));
    }, path + ".under_relaxation", "omega > 1");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.vaporPhase_ = c.liquidPhase_; }));
    }, path + ".vapor_phase", "same phase for both roles");
    checkThrows([&] {
        validate(mutated([](auto& c) { c.liquidPhase_.clear(); }));
    }, "liquid_phase and vapor_phase", "missing role");

    // the model constructor rejects invalid parameters as well
    checkThrows([&] {
        createMassTransferModel(
            mutated([](auto& c) { c.underRelaxation_ = 2.0; }));
    }, "under_relaxation", "factory validates");
}

void testRates()
{
    const auto model = createMassTransferModel(refConfig());

    // vaporization: p < p_v -> mdot > 0, matches the analytic value
    const double mv = model->rate({2000.0, 0.1, rhoL, rhoV});
    check(mv > 0.0, "vaporization sign");
    checkNear(mv, analyticVap(2000.0, 0.1), 1e-12, "vaporization value");

    // condensation: p > p_v -> mdot < 0
    const double mc = model->rate({5000.0, 0.1, rhoL, rhoV});
    check(mc < 0.0, "condensation sign");
    checkNear(mc, analyticCond(5000.0, 0.1), 1e-12, "condensation value");

    // equilibrium
    check(model->rate({pv, 0.1, rhoL, rhoV}) == 0.0, "equilibrium");

    // vaporization decreases to zero as alpha_v -> 1
    double prev = model->rate({2000.0, 0.0, rhoL, rhoV});
    bool monotone = true;
    for (double av : {0.25, 0.5, 0.75, 1.0})
    {
        const double m = model->rate({2000.0, av, rhoL, rhoV});
        monotone = monotone && m < prev;
        prev = m;
    }
    check(monotone, "vaporization decreases with alpha_v");
    check(model->rate({2000.0, 1.0, rhoL, rhoV}) == 0.0,
          "vaporization vanishes at alpha_v = 1");

    // condensation vanishes as alpha_v -> 0 and grows with alpha_v
    check(model->rate({5000.0, 0.0, rhoL, rhoV}) == 0.0,
          "condensation vanishes at alpha_v = 0");
    check(model->rate({5000.0, 0.5, rhoL, rhoV}) <
              model->rate({5000.0, 0.1, rhoL, rhoV}),
          "condensation magnitude grows with alpha_v");

    // alpha_v is clamped to [0,1]
    checkNear(model->rate({2000.0, -0.3, rhoL, rhoV}),
              model->rate({2000.0, 0.0, rhoL, rhoV}),
              1e-14,
              "alpha_v < 0 clamped");
    check(model->rate({5000.0, 1.7, rhoL, rhoV}) ==
              model->rate({5000.0, 1.0, rhoL, rhoV}),
          "alpha_v > 1 clamped");

    // invalid densities -> no rate, no NaN
    check(model->rate({2000.0, 0.1, 0.0, rhoV}) == 0.0, "rho_l = 0");
    check(model->rate({2000.0, 0.1, rhoL, 0.0}) == 0.0, "rho_v = 0");

    // negative absolute pressure: true pressure by default, clipped on request
    check(model->rate({-1000.0, 0.1, rhoL, rhoV}) >
              model->rate({0.0, 0.1, rhoL, rhoV}),
          "true pressure used without clipping");
    auto cfg = refConfig();
    cfg.pressureClippingForRate_ = true;
    RayleighPlessetCavitationMassTransfer clipped(cfg);
    checkNear(clipped.rate({-1000.0, 0.1, rhoL, rhoV}),
              clipped.rate({0.0, 0.1, rhoL, rhoV}),
              1e-14,
              "pressure clipping for rate");

    // array interface
    const double p[3] = {2000.0, pv, 5000.0};
    const double av[3] = {0.1, 0.1, 0.1};
    const double rl[3] = {rhoL, rhoL, rhoL};
    const double rv[3] = {rhoV, rhoV, rhoV};
    double out[3];
    model->compute(p, av, rl, rv, out);
    check(out[0] > 0.0 && out[1] == 0.0 && out[2] < 0.0,
          "compute() fills mdot");
}

void testCouplingHelpers()
{
    const auto model = createMassTransferModel(refConfig());

    // continuity source: div(u) = mdot (1/rho_v - 1/rho_l)
    const double mv = model->rate({2000.0, 0.1, rhoL, rhoV});
    const double Dv = massTransferModel::volumeSource(mv, rhoL, rhoV);
    check(Dv > 0.0, "expansion on vaporization (rho_l > rho_v)");
    checkNear(Dv, mv * (1.0 / rhoV - 1.0 / rhoL), 1e-14, "continuity source");
    const double mc = model->rate({5000.0, 0.1, rhoL, rhoV});
    check(massTransferModel::volumeSource(mc, rhoL, rhoV) < 0.0,
          "contraction on condensation");
    check(massTransferModel::volumeSource(mv, 0.0, rhoV) == 0.0,
          "volumeSource guards invalid density");

    // phase sources: S_alpha_v = mdot/rho_v, S_alpha_l = -mdot/rho_l; their sum
    // is the volumetric source
    checkNear(mv / rhoV - mv / rhoL, Dv, 1e-14, "sum of alpha sources = D");

    // under-relaxation: omega*raw + (1-omega)*prev
    checkNear(model->relax(8.0, 4.0), 0.25 * 8.0 + 0.75 * 4.0, 1e-14, "relax");
    check(model->includeContinuitySource(), "continuity switch from config");
}

} // namespace

int main()
{
    testDefaults();
    testValidation();
    testRates();
    testCouplingHelpers();

    std::cout << nChecks - nFailed << "/" << nChecks << " checks passed"
              << std::endl;
    return nFailed == 0 ? 0 : 1;
}
