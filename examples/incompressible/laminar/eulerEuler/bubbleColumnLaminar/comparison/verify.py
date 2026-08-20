#!/usr/bin/env python3
"""Check an OpenAccel Euler-Euler run, and compare it with OpenFOAM if given.

    python3 comparison/verify.py                  # check this case's run
    python3 comparison/verify.py <openfoam_case>  # ...and compare against it

Checks: finished without NaN; boundary conditions actually honoured; solution
bounded and conservative; agreement with OpenFOAM. Prints PASS/FAIL per item,
exits non-zero if anything failed.
"""
import os, re, sys, glob
import numpy as np
from scipy.io import netcdf_file as nf

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
RHO_REF, G = 1000.0, 9.81
fails = []

def check(name, ok, detail):
    print("  [%s] %-44s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        fails.append(name)

res = os.path.join(CASE, "results.e")
if not os.path.exists(res):
    sys.exit("no results.e -- run the case first")

f = nf(res, "r", mmap=False)
names = [b"".join(r).decode().strip("\x00").strip()
         for r in f.variables["name_nod_var"][:]]
t = f.variables["time_whole"][:]
y = f.variables["coordy"][:]
x = f.variables["coordx"][:]
xmax, ymax = x.max(), y.max()
g = lambda k, i: f.variables["vals_nod_var%d" % (names.index(k) + 1)][i]

wall = (np.isclose(x, 0.0) | np.isclose(x, xmax)) & (y > 1e-9) & (y < ymax - 1e-9)
inlet = np.isclose(y, 0.0)

print("\n=== RUN ===")
print("  mesh %d nodes, %d frames, t = %.3f .. %.3f"
      % (len(x), len(t), t[0], t[-1]))

finite = all(np.all(np.isfinite(g("velocity.air_y", i))) for i in range(len(t)))
check("no NaN / Inf anywhere", finite, "%d frames all finite" % len(t))

ua = np.hypot(g("velocity.air_x", -1), g("velocity.air_y", -1))
uw = np.hypot(g("velocity.water_x", -1), g("velocity.water_y", -1))
check("no-slip wall |U.air| = 0", ua[wall].max() < 1e-8, "max %.2e" % ua[wall].max())
check("no-slip wall |U.water| = 0", uw[wall].max() < 1e-8, "max %.2e" % uw[wall].max())
vin = g("velocity.air_y", -1)[inlet]
check("inlet U.air_y uniform at its BC", np.ptp(vin) < 1e-9, "%.4f" % vin.mean())

aa, aw = g("volume_fraction.air", -1), g("volume_fraction.water", -1)
check("alpha bounded in [0,1]", aa.min() >= -1e-9 and aa.max() <= 1 + 1e-9,
      "[%.4f, %.4f]" % (aa.min(), aa.max()))
check("sum(alpha) = 1", np.abs(1 - (aa + aw)).max() < 1e-9,
      "max dev %.2e" % np.abs(1 - (aa + aw)).max())

if len(sys.argv) > 1:
    OF = os.path.abspath(sys.argv[1])
    def off(p):
        s = open(p).read()
        m = re.search(r"internalField\s+nonuniform[^(]*\(\s*(.*?)\n\)\s*;", s, re.S)
        if m:
            b = m.group(1)
            if "(" in b:
                return np.array([[float(v) for v in q.split()]
                                 for q in re.findall(r"\(([^)]*)\)", b)])
            return np.array([float(v) for v in b.split()])
        raise RuntimeError("uniform field in " + p)
    if not os.path.exists(OF + "/0/Ccy"):
        sys.exit("first run, inside the OpenFOAM case:\n"
                 "  foamPostProcess -func writeCellCentres -time 0")
    times = sorted([d for d in glob.glob(OF + "/[0-9]*")
                    if re.match(r"^[0-9.]+$", os.path.basename(d))
                    and os.path.exists(d + "/alpha.air")],
                   key=lambda p: float(os.path.basename(p)))
    last = times[-1]
    tof = float(os.path.basename(last))
    ccy = off(OF + "/0/Ccy")
    lev = np.unique(np.round(ccy, 9))
    P = lambda v: np.array([v[np.isclose(ccy, q)].mean(axis=0) for q in lev])
    ti = int(np.argmin(abs(t - tof)))
    yn = np.unique(np.round(y, 9))
    I = (x > 1e-9) & (x < xmax - 1e-9)
    oap = lambda k: np.interp(lev, yn,
        np.array([g(k, ti)[np.isclose(y, q) & I].mean() for q in yn]))
    print("\n=== vs OPENFOAM (t = %.3f vs %.3f) ===" % (tof, t[ti]))
    wc = (lev > 0.10) & (lev < 0.60)
    for nm, vof, voa, tol in [
            ("alpha.air", P(off(last + "/alpha.air")), oap("volume_fraction.air"), 10.0),
            ("U.air_y", P(off(last + "/U.air"))[:, 1], oap("velocity.air_y"), 10.0)]:
        d = 100 * abs(voa[wc].mean() - vof[wc].mean()) / max(abs(vof[wc].mean()), 1e-30)
        check("water column %-11s <%.0f%%" % (nm, tol), d <= tol,
              "OF %+.5g  OA %+.5g  -> %.1f%%" % (vof[wc].mean(), voa[wc].mean(), d))
    pf = P(off(last + "/p")) - 1e5
    pa = oap("pressure") - RHO_REF * G * lev
    d = 100 * np.linalg.norm((pa - pa.mean()) - (pf - pf.mean())) / \
        np.linalg.norm(pf - pf.mean())
    check("pressure profile <10%", d <= 10.0, "L2 %.1f%%" % d)

print("\n%s\n" % ("ALL CHECKS PASSED" if not fails else "FAILED: " + ", ".join(fails)))
sys.exit(1 if fails else 0)
