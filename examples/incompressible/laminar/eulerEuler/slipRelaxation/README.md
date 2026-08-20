# Euler-Euler slip-relaxation verification case

Smallest useful verification case for the inhomogeneous Euler-Euler model. Two
incompressible phases fill a 1 m x 1 m periodic box (10 x 10 QUAD4 elements) with
uniform constant volume fractions, no gravity, no pressure gradient and no mass
transfer. The phases start at different uniform velocities and are coupled only
by a constant interphase drag coefficient `K`.

Every field stays spatially uniform, so the discretization contributes no error.
Any deviation from the analytic solution is an assembly error, which is what
makes this case diagnostic rather than merely indicative.

## Analytic solution

With `m_k = alpha_k * rho_k`, the phasic momentum equations reduce to

```
m_1 dU_1/dt = -K (U_1 - U_2)
m_2 dU_2/dt = +K (U_1 - U_2)
```

so the mixture momentum is conserved exactly and the slip decays exponentially:

```
tau    = 1 / ( K * (1/m_1 + 1/m_2) )
u_r(t) = u_r(0) * exp(-t/tau)
U_inf  = (m_1 U_1(0) + m_2 U_2(0)) / (m_1 + m_2)
U_1(t) = U_inf + (m_2/(m_1+m_2)) u_r(t)
U_2(t) = U_inf - (m_1/(m_1+m_2)) u_r(t)
```

This case uses `alpha_1 = 0.4, rho_1 = 1000` (m_1 = 400),
`alpha_2 = 0.6, rho_2 = 100` (m_2 = 60), `U_1(0) = (1,0)`, `U_2(0) = (0,0)`,
`K = 500`, giving `tau = 0.10434783 s` and `U_inf = 0.86956522 m/s`.
Total time 0.6 s is about 5.75 tau.

## Running

```sh
python3 make_mesh.py                       # regenerate mesh.e (optional)
../../../../../build-2d/openaccel-2D.exe -i input.i
python3 verify.py                          # acceptance checks
```

Note the `-i` flag: a bare positional argument is ignored and the solver falls
back to its default `input.i`.

Via the regression runner (the runner defaults to `build2D`, so point it at
this build):

```sh
python3 tools/python/regression_tests/quick_test.py slipRelaxation \
        --exe-2d build-2d/openaccel-2D.exe
```

`run.log` reaches ~150 MB at the default `verbose: 1` because every outer
iteration of all 2400 timesteps is printed. Redirect it to /dev/null, or lower
`verbose`, if you do not need the solver trace.

## Acceptance checks

`verify.py` asserts every item the project plan lists for this case:

| Plan requirement | Check |
| --- | --- |
| Conservation of `alpha_1 rho_1 U_1 + alpha_2 rho_2 U_2` | relative drift of the mixture momentum |
| Exponential decay of `U_1 - U_2` | pointwise vs the analytic exponential |
| Equal and opposite drag forces | sum over phases of `interphase_momentum_source`, componentwise, per node |
| Constant volume fractions | drift of each `alpha_k` from its initial value |
| Uniform pressure | `max(abs(p)) / (rho U^2)` |
| Near-zero phasic mass-divergence residuals | `max(abs(divergence.<phase>)) / mdot_scale` |
| `sum(alpha_k) = 1` | `max(abs(alpha_1 + alpha_2 - 1))` |

Pressure and mass divergence are dimensional, so they are judged as
dimensionless ratios against this case's own scales (`rho U^2` = 1000 Pa and a
face mass flow rate = 40 kg/s) rather than against absolute numbers.

Serial/MPI equivalence is the one plan item not automated -- run under
`mpirun` and diff the `verify.py` table.

## Status: passing (9/9)

```
     t [s]     U1 (num)   U1 (exact)     U2 (num)   U2 (exact)      |p|max
    0.0000     1.000000     1.000000     0.000000     0.000000   0.000e+00
    0.0050     0.993905     0.993897     0.040637     0.040684   5.746e-05
    ...
    0.6000     0.869983     0.869980     0.866778     0.866798   1.622e-07

  [PASS] all values finite                          ok
  [PASS] exponential slip decay matches analytic    max err 3.828e-04
  [PASS] mixture momentum conserved                 max rel drift 2.418e-10
  [PASS] volume fractions constant                  max drift 4.111e-10
  [PASS] sum(alpha_k) == 1                          max |sum-1| 0.000e+00
  [PASS] velocity/alpha fields spatially uniform    max spread 7.144e-09
  [PASS] pressure uniform (zero to within rho*U^2)  max |p|/(rho U^2) 5.746e-08
  [PASS] drag forces equal and opposite             max |sum|/|M| 0.000e+00
  [PASS] phasic mass divergence ~ 0                 max |div|/mdot 4.940e-10
```

Timestep refinement confirms the residual slip-decay error is purely
first-order backward Euler discretization:

| dt | max abs(U1 - exact) | ratio |
| --- | --- | --- |
| 1.00e-3 | 2.2878e-04 | |
| 5.00e-4 | 1.1463e-04 | 1.996 |
| 2.50e-4 | 5.7371e-05 | 1.998 |
| 1.25e-4 | 2.8700e-05 | 1.999 |

Fitting the decay constant from the data gives `tau = 0.10447 s` against the
exact `0.10434783 s` (0.12 %, consistent with the dt above).

### Note on under-relaxation

This case **requires** SIMPLE under-relaxation. Without
`relaxation_parameters` the outer loop over-corrects and the pressure diverges
by a fixed factor of ~3600 per timestep -- independent of `timestep` and of the
drag coefficient -- while momentum, drag and the mass fluxes remain individually
correct. The failure looks like a solver bug but is purely a missing
under-relaxation setting, so do not remove that block when adapting this case.

### Note on Exodus name truncation

Exodus limits variable names to 32 characters, so all four components of
`interphase_momentum_source.<phase>` are written under the same truncated name.
`verify.py` therefore locates them by position (phase-major ordering), which it
validated against the physics: phase 1, which loses momentum, carries the
negative x-component.
