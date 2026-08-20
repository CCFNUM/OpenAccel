#!/usr/bin/env python3
"""Acceptance checks for the Euler-Euler batch settling column.

Run the case first, then:  python3 verify.py [results.e]

The column is closed, so at steady state the net volumetric flux vanishes,
alpha_d U_d + alpha_c U_c = 0. Dividing each phasic momentum balance by its own
alpha and subtracting eliminates grad(p):

    (rho_d - rho_c) g = K s (1/alpha_d + 1/alpha_c),   s = U_d - U_c

With Schiller-Naumann, K = 0.75 C_D rho_c alpha_c alpha_d |s| / d, every alpha
cancels and the slip is exactly the single-particle terminal velocity:

    0.75 C_D(Re) rho_c s^2 / d = (rho_d - rho_c) g,   Re = rho_c |s| d / mu

    U_d = -(1 - alpha_d) u_t     U_c = +alpha_d u_t
"""
import glob
import math
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

# --- case parameters, must match input.i -----------------------------------
RHO_D, RHO_C = 2500.0, 1000.0
MU_C = 1.0e-3
DIAM = 5.0e-5
G = 9.81
ALPHA_D0 = 0.05
HEIGHT = 0.05

# bulk sampling band: away from the walls and from the settling fronts
BAND = (0.018, 0.032)

SLIP_TOL = 2.0e-2    # relative, vs the implicit Schiller-Naumann terminal vel.
FLUX_TOL = 1.0e-6    # |alpha_d U_d + alpha_c U_c| / (alpha_d u_t)


def terminal_velocity():
    """Solve 0.75 C_D(Re) rho_c s^2 / d = (rho_d - rho_c) g for s."""
    s = (RHO_D - RHO_C) * G * DIAM ** 2 / (18.0 * MU_C)   # Stokes start
    for _ in range(200):
        re = RHO_C * s * DIAM / MU_C
        cd = 24.0 / re * (1.0 + 0.15 * re ** 0.687) if re <= 1000.0 else 0.44
        s = math.sqrt((RHO_D - RHO_C) * G * DIAM / (0.75 * cd * RHO_C))
    return s, RHO_C * s * DIAM / MU_C


def open_results(path="results.e"):
    parts = sorted(glob.glob(f"{path}.*.*")) if not Path(path).exists() else [path]
    if not parts:
        raise SystemExit(f"no results found: neither {path} nor {path}.<n>.<rank>")
    datasets = [nc.Dataset(p) for p in parts]
    for d in datasets:
        d.set_auto_mask(False)
    return datasets


def main(path="results.e"):
    datasets = open_results(path)
    d = datasets[0]
    names = [n.strip() for n in nc.chartostring(d.variables["name_nod_var"][:])]
    t = np.asarray(d.variables["time_whole"][:])

    def cat(key):
        return np.concatenate([np.asarray(ds.variables[key][:]) for ds in datasets])

    def field(name, step):
        i = names.index(name) + 1
        return np.concatenate([
            np.asarray(ds.variables[f"vals_nod_var{i}"][step, :]) for ds in datasets
        ])

    y = cat("coordy")
    band = (y > BAND[0]) & (y < BAND[1])
    if not band.any():
        raise SystemExit("bulk sampling band selected no nodes")

    u_t, re = terminal_velocity()
    ud_exact, uc_exact = -(1.0 - ALPHA_D0) * u_t, ALPHA_D0 * u_t
    print(f"u_t = {u_t:.6e} m/s   Re = {re:.4f}")
    print(f"exact:  U_d = {ud_exact:.6e}   U_c = {uc_exact:.6e}   "
          f"slip = {-u_t:.6e}\n")
    print(f"{'t [s]':>8} {'U_d bulk':>13} {'U_c bulk':>13} {'slip':>13} "
          f"{'alpha_d':>9} {'net flux':>11}")

    worst_slip = float("inf")
    worst_flux = 0.0
    saw_nan = False
    for step in range(len(t)):
        ud = field("velocity.particles_y", step)[band]
        uc = field("velocity.water_y", step)[band]
        ad = field("volume_fraction.particles", step)[band]
        if not (np.all(np.isfinite(ud)) and np.all(np.isfinite(uc))):
            saw_nan = True
            print(f"{t[step]:8.3f}  *** non-finite values ***")
            break
        udm, ucm, adm = float(ud.mean()), float(uc.mean()), float(ad.mean())
        slip = udm - ucm
        flux = adm * udm + (1.0 - adm) * ucm
        print(f"{t[step]:8.3f} {udm:13.6e} {ucm:13.6e} {slip:13.6e} "
              f"{adm:9.5f} {flux:11.3e}")
        # the slip is only meaningful once the drag transient has passed
        if t[step] > 0.05:
            worst_slip = min(worst_slip, abs(abs(slip) - u_t) / u_t)
            worst_flux = max(worst_flux, abs(flux) / (ALPHA_D0 * u_t))

    print()
    note = "" if not saw_nan else "  (inconclusive: run went non-finite)"
    checks = [
        ("all values finite", not saw_nan,
         "NaN/Inf encountered" if saw_nan else "ok"),
        ("slip equals Schiller-Naumann terminal velocity",
         not saw_nan and worst_slip < SLIP_TOL,
         f"best rel err {worst_slip:.3e}{note}"),
        ("zero net volumetric flux (closed column)",
         not saw_nan and worst_flux < FLUX_TOL,
         f"max |flux| {worst_flux:.3e}{note}"),
    ]
    width = max(len(c[0]) for c in checks)
    failed = sum(0 if ok else 1 for _, ok, _ in checks)
    for label, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label:<{width}}  {detail}")
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
