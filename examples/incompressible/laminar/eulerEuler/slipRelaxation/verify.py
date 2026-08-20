#!/usr/bin/env python3
"""Acceptance checks for the Euler-Euler slip-relaxation verification case.

Run the case first, then:  python3 verify.py [results.e]

Every field stays spatially uniform in this case, so with m_k = alpha_k rho_k
the phasic momentum equations reduce to two coupled ODEs whose solution is

    u_r(t) = u_r(0) exp(-t/tau),   tau = 1 / (K (1/m_1 + 1/m_2))
    U_1(t) = U_inf + (m_2/(m_1+m_2)) u_r(t)
    U_2(t) = U_inf - (m_1/(m_1+m_2)) u_r(t)
    U_inf  = (m_1 U_1(0) + m_2 U_2(0)) / (m_1 + m_2)

and the mixture momentum m_1 U_1 + m_2 U_2 is conserved exactly.
"""
import glob
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

# --- case parameters, must match input.i -----------------------------------
ALPHA1, RHO1 = 0.4, 1000.0
ALPHA2, RHO2 = 0.6, 100.0
U1_0, U2_0 = 1.0, 0.0
K = 500.0

M1, M2 = ALPHA1 * RHO1, ALPHA2 * RHO2
TAU = 1.0 / (K * (1.0 / M1 + 1.0 / M2))
U_INF = (M1 * U1_0 + M2 * U2_0) / (M1 + M2)
MOM0 = M1 * U1_0 + M2 * U2_0

# Tolerances. The time integration is first-order backward Euler, so the
# velocity error is O(dt/tau); UNIFORM_TOL and SUM_TOL are discretization-free
# and therefore tight.
DECAY_TOL = 2.0e-2   # relative, vs analytic exponential
MOM_TOL = 1.0e-6     # relative drift of mixture momentum
# Spatial spread of velocity and volume fraction. The exact solution is
# uniform; the outer loop leaves a small residual pressure whose velocity
# correction imprints a spread at the linear-solver tolerance (observed ~7e-9
# at rtol 1e-4), far below any physical signal.
UNIFORM_TOL = 1.0e-7
SUM_TOL = 1.0e-10    # |sum(alpha) - 1|
ALPHA_TOL = 1.0e-9   # drift of alpha from its initial value
DRAG_TOL = 1.0e-12   # |sum_k M_k| / |M|, exact cancellation expected

# Pressure and mass divergence are dimensional, so they are judged against the
# physical scales of this case rather than against absolute numbers:
#   P_SCALE    ~ rho U^2   dynamic pressure of the initial slip   (1000 Pa)
#   MDOT_SCALE ~ (alpha rho) U h   a face mass flow rate          (40 kg/s)
P_SCALE = max(RHO1, RHO2) * (U1_0 - U2_0) ** 2
CELL_SIZE = 0.1                                          # 10 elements over 1 m
MDOT_SCALE = max(M1, M2) * abs(U1_0 - U2_0) * CELL_SIZE
PRESSURE_TOL = 1.0e-6  # |p|max / P_SCALE
DIV_TOL = 1.0e-9       # |div|max / MDOT_SCALE


def open_results(path="results.e"):
    """Open the run output, serial or decomposed.

    A parallel run writes one Exodus file per rank (`results.e.<n>.<r>`)
    instead of a single `results.e`, so both layouts are accepted and the
    per-rank nodal values are concatenated. Nodes shared between ranks appear
    more than once, which is harmless: every check here is a min/max/mean over
    a field that is spatially uniform in the exact solution.
    """
    parts = sorted(glob.glob(f"{path}.*.*")) if not Path(path).exists() else [path]
    if not parts:
        raise SystemExit(f"no results found: neither {path} nor {path}.<n>.<rank>")
    datasets = [nc.Dataset(p) for p in parts]
    for d in datasets:
        # Without this, netCDF4 returns masked arrays and silently replaces NaN
        # with the fill value, hiding a diverged run behind clean zeros.
        d.set_auto_mask(False)
    return datasets, len(parts)


def main(path="results.e"):
    datasets, n_ranks = open_results(path)
    d = datasets[0]
    names = [n.strip() for n in nc.chartostring(d.variables["name_nod_var"][:])]
    t = np.asarray(d.variables["time_whole"][:])
    if n_ranks > 1:
        print(f"reading {n_ranks} decomposed result files\n")

    def var(index, step):
        return np.concatenate([
            np.asarray(ds.variables[f"vals_nod_var{index + 1}"][step, :])
            for ds in datasets
        ])

    def field(name, step):
        return var(names.index(name), step)

    # The Exodus name limit (32 chars) truncates every
    # `interphase_momentum_source.<phase>_<comp>` variable to the same string,
    # so these are located by position, not name. They are written phase-major:
    # [phase1_x, phase1_y, ..., phase2_x, phase2_y, ...]. Verified against the
    # physics: phase 1 (losing momentum) carries the negative x-component.
    drag_idx = [i for i, n in enumerate(names)
                if n.startswith("interphase_momentum_source")]
    n_phases, n_comp = 2, len(drag_idx) // 2 if drag_idx else 0
    if drag_idx and len(drag_idx) != n_phases * n_comp:
        raise SystemExit(f"unexpected interphase_momentum_source layout: "
                         f"{len(drag_idx)} components found")
    div_names = [n for n in names if n.startswith("divergence.")]

    checks = []      # (label, ok, detail)
    worst_uniform = 0.0
    worst_sum = 0.0
    worst_alpha = 0.0
    worst_mom = 0.0
    worst_decay = 0.0
    worst_drag = 0.0
    worst_div = 0.0
    worst_p = 0.0
    saw_nan = False

    print(f"tau = {TAU:.6f} s    U_inf = {U_INF:.6f} m/s    "
          f"mixture momentum = {MOM0:.6f}\n")
    print(f"{'t [s]':>10} {'U1 (num)':>12} {'U1 (exact)':>12} "
          f"{'U2 (num)':>12} {'U2 (exact)':>12} {'|p|max':>11}")

    for step in range(len(t)):
        u1 = field("velocity.fluid1_x", step)
        u2 = field("velocity.fluid2_x", step)
        a1 = field("volume_fraction.fluid1", step)
        a2 = field("volume_fraction.fluid2", step)
        p = field("pressure", step)

        if not all(np.all(np.isfinite(v)) for v in (u1, u2, a1, a2, p)):
            saw_nan = True
            print(f"{t[step]:10.4f}  *** non-finite values ***")
            break

        # uniformity (peak-to-peak spread); pressure is judged separately
        # below against its physical scale, since it is dimensional
        worst_p = max(worst_p, float(np.max(np.abs(p))))
        for v in (u1, u2, a1, a2):
            worst_uniform = max(worst_uniform, float(np.ptp(v)))

        # equal-and-opposite drag: sum of the interphase momentum source over
        # all phases must vanish, componentwise, at every node
        if drag_idx:
            raw = [var(i, step) for i in drag_idx]
            magnitude = max(float(np.max(np.abs(v))) for v in raw)
            for c in range(n_comp):
                total = sum(raw[ph * n_comp + c] for ph in range(n_phases))
                if magnitude > 0.0:
                    worst_drag = max(worst_drag,
                                     float(np.max(np.abs(total))) / magnitude)

        # phasic mass-divergence residuals
        for nm in div_names:
            worst_div = max(worst_div, float(np.max(np.abs(field(nm, step)))))

        worst_sum = max(worst_sum, float(np.max(np.abs(a1 + a2 - 1.0))))
        worst_alpha = max(worst_alpha,
                          float(np.max(np.abs(a1 - ALPHA1))),
                          float(np.max(np.abs(a2 - ALPHA2))))

        u1m, u2m = float(np.mean(u1)), float(np.mean(u2))
        mom = M1 * u1m + M2 * u2m
        worst_mom = max(worst_mom, abs(mom - MOM0) / abs(MOM0))

        ur = (U1_0 - U2_0) * np.exp(-t[step] / TAU)
        u1e, u2e = U_INF + (M2 / (M1 + M2)) * ur, U_INF - (M1 / (M1 + M2)) * ur
        scale = max(abs(U1_0 - U2_0), 1.0e-12)
        worst_decay = max(worst_decay,
                          abs(u1m - u1e) / scale, abs(u2m - u2e) / scale)

        print(f"{t[step]:10.4f} {u1m:12.6f} {u1e:12.6f} "
              f"{u2m:12.6f} {u2e:12.6f} {np.max(np.abs(p)):11.3e}")

    print()
    # A run that went non-finite is stopped at that point, so the remaining
    # checks have only seen the steps before the blow-up. Reporting those as
    # PASS would be misleading -- they are inconclusive, not satisfied.
    evaluated = "" if not saw_nan else "  (inconclusive: run went non-finite)"

    checks.append(("all values finite", not saw_nan,
                   "NaN/Inf encountered" if saw_nan else "ok"))
    checks.append(("exponential slip decay matches analytic",
                   not saw_nan and worst_decay < DECAY_TOL,
                   f"max err {worst_decay:.3e}{evaluated}"))
    checks.append(("mixture momentum conserved",
                   not saw_nan and worst_mom < MOM_TOL,
                   f"max rel drift {worst_mom:.3e}{evaluated}"))
    checks.append(("volume fractions constant",
                   not saw_nan and worst_alpha < ALPHA_TOL,
                   f"max drift {worst_alpha:.3e}{evaluated}"))
    checks.append(("sum(alpha_k) == 1",
                   not saw_nan and worst_sum < SUM_TOL,
                   f"max |sum-1| {worst_sum:.3e}{evaluated}"))
    checks.append(("velocity/alpha fields spatially uniform",
                   not saw_nan and worst_uniform < UNIFORM_TOL,
                   f"max spread {worst_uniform:.3e}{evaluated}"))
    checks.append(("pressure uniform (zero to within rho*U^2)",
                   not saw_nan and worst_p / P_SCALE < PRESSURE_TOL,
                   f"max |p|/(rho U^2) {worst_p / P_SCALE:.3e}{evaluated}"))
    checks.append(("drag forces equal and opposite",
                   not saw_nan and bool(drag_idx) and worst_drag < DRAG_TOL,
                   ("no interphase_momentum_source in output"
                    if not drag_idx
                    else f"max |sum|/|M| {worst_drag:.3e}{evaluated}")))
    checks.append(("phasic mass divergence ~ 0",
                   not saw_nan and bool(div_names) and worst_div / MDOT_SCALE < DIV_TOL,
                   ("no divergence.* in output" if not div_names
                    else f"max |div|/mdot {worst_div / MDOT_SCALE:.3e}"
                         f"{evaluated}")))

    width = max(len(c[0]) for c in checks)
    failed = 0
    for label, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label:<{width}}  {detail}")
        failed += 0 if ok else 1

    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
