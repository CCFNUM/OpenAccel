# Euler-Euler batch settling column

A dilute suspension of heavy particles (`alpha = 0.05`, `rho = 2500`,
`d = 50 um`) in water settles under gravity in a closed 4 mm x 50 mm column
(4 x 100 QUAD4). Only drag couples the phases.

This is the first case that exercises **gravity, buoyancy, the `alpha*grad(p)`
split, Schiller-Naumann drag and genuine volume-fraction transport** — none of
which `slipRelaxation` touches.

## Analytic solution

The column is closed, so the net volumetric flux vanishes:
`alpha_d U_d + alpha_c U_c = 0`. Dividing each phasic momentum balance by its
own `alpha` and subtracting eliminates `grad(p)`:

```
(rho_d - rho_c) g = K s (1/alpha_d + 1/alpha_c),    s = U_d - U_c
```

With Schiller-Naumann, `K = 0.75 C_D rho_c alpha_c alpha_d |s| / d`, so every
`alpha` cancels and the slip is exactly the single-particle terminal velocity:

```
0.75 C_D(Re) rho_c s^2 / d = (rho_d - rho_c) g,     Re = rho_c |s| d / mu
U_d = -(1 - alpha_d) u_t        U_c = +alpha_d u_t
```

For these numbers `u_t = 1.982957e-3 m/s` at `Re = 0.0991`, **2.97 % below the
Stokes value** — so the case genuinely tests `C_D(Re)`, not merely Stokes drag.
The clear-water interface descends at `|U_d| = 1.883809e-3 m/s`; full settling
takes ~26.5 s.

## Important: this is *not* a Richardson-Zaki validation

The model as implemented yields `U_d / u_t = (1 - alpha_d)^1` — a hindered
settling exponent of exactly **1**. The empirical Richardson-Zaki exponent is
`n ~ 4.65` in the Stokes regime.

That difference is a **missing swarm / hindrance correction to the drag law**,
not a solver defect: plain Schiller-Naumann is a single-particle correlation.
Reproducing Richardson-Zaki requires a hindered-drag model (Wen-Yu, Gidaspow,
or an explicit `(1-alpha)^-n` factor), which is not implemented.

So treat this case as **verification against the analytic solution of the
implemented model**, not as validation against Richardson-Zaki. Comparing
directly to `n = 4.65` will look like a large error when nothing is wrong with
the solver. Adding a hindered-drag model is the prerequisite for genuine
Richardson-Zaki validation.

## Running

```sh
python3 make_mesh.py                       # regenerate mesh.e (optional)
../../../../../build-2d/openaccel-2D.exe -i input.i
python3 verify.py
```

## Acceptance checks

`verify.py` samples a bulk band (`0.018 < y < 0.032`), away from the walls and
the settling fronts, and checks:

- slip velocity equals the implicitly-solved Schiller-Naumann terminal velocity
- zero net volumetric flux, `alpha_d U_d + alpha_c U_c = 0`
- all values finite

## Status: FAILING — two bugs fixed, root cause NOT found

### Fixed and verified

1. **Buoyancy reference density.** OpenAccel solves for the reduced pressure
   (`flowModel::computeBodyForces` uses `(rho - rho_ref) g`). The Euler-Euler
   assembler applied the full `alpha_k rho_k g`. Corrected in
   `phasicNavierStokesAssemblerNodeTerms.cpp`. Verified: with equal phase
   densities the body force is exactly zero and the case is now exactly at rest
   (`0.000e+00` every step); before, it diverged on step 1.

2. **Body force never reached momentum.** `eulerEulerModel::updatePhaseBodyForces`
   seeded `body_forces`, but `flowModel::redistributeBodyForces` zeroes that
   field and gathers from `body_forces_original`. The seed must go into
   FOriginal. Verified: V-Mom RHS is now `4.20e-04`, matching
   `alpha (rho-rho_ref) g V`; it was `1.09e-36` (i.e. absent).

### Case-setup finding

The default `pressure: 0` is the *absolute* pressure, and `initializePressure`
adds `rho_ref |g| y` to form the stored reduced pressure. A uniform absolute
pressure is therefore NOT hydrostatic equilibrium -- it puts the whole column in
free fall (measured Uy = -0.0168 = 9.81*dt at iteration 1). The equilibrium
input is `pressure: "-rho_mix*|g|*y"`, i.e. `-10545.75*y` here. With that, outer
iterations 1-3 converge correctly toward the analytic slip.

### The remaining failure

From iteration ~4 the outer loop blows up geometrically at a **fixed gain of
~2.8 per outer iteration**, velocities oscillating symmetrically, until a phase
hits the residual-alpha floor and its momentum equation degenerates.

Ruled out by direct experiment -- do not re-test these:

| hypothesis | evidence against |
| --- | --- |
| singular pressure matrix | pin verified at row 148; `rowsNotSummingToZero=1`; diag ratio 4; RHS compatibility 3.8e-14 |
| linear solver | all systems solve to Drop ~1e-16 |
| direct LU as a fix | LU is BROKEN in this build -- corrupts even the trivial case |
| timestep / CFL | 100x reduction: identical |
| advection scheme | upwind vs high_resolution: identical |
| under-relaxation | 7x range (0.7/0.3 down to 0.1/0.03): identical |
| drag coupling / partial elimination | K/a_P from 84.5% to 0.1%: identical |
| wall vs symmetry boundaries | identical |
| initial pressure | equilibrium IC: converges 3 iterations then same blow-up |
| body-force magnitude | excess density 1500 -> 1 kg/m3: identical |
| body-force redistribution | `body_force_redistribution: false`: identical |

The decisive pair of observations:

- `gravity: [0, 0]` with the buoyant path active -> **survives**, exactly at rest
- any non-zero gravity -> diverges identically, *independent of magnitude*

A loop gain that is invariant to the forcing magnitude but requires the forcing
to be non-zero is the signature of a **discretely inconsistent operator** that
gravity switches on -- not a stiffness or stability-margin problem. The error
grows geometrically from whatever seed exists.

### Comparison with OpenFOAM (`multiphaseEuler/cellPressureCorrector.C`)

Two structural differences were identified:

1. OpenFOAM adds the body force to the **face flux**
   (`phase.URef() = HbyADs + reconstruct(alphaByADfs*mSfGradp - FgByADfs)`),
   with `Fgfs` built on faces from `ghf*snGrad(rho)*magSf` and
   `interpolate(rho_k - rho_mix)*(g & Sf)`. OpenAccel's
   `updatePhaseMassFluxInterior_` has only the Rhie-Chow pressure term and no
   body-force term. **Untested** -- the analytic justification is unclear for
   constant phase density, and disabling redistribution (which would restore
   exact cancellation) did not help.
2. OpenFOAM's `invADV` is a phase-**coupled** matrix; OpenAccel uses a per-phase
   scalar `du_k`. **Ruled out** by the drag sweep above.

Difference 1 is the only untested lead.
