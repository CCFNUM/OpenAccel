# Euler-Euler bubble column (OpenFOAM cross-check case)

OpenAccel version of `$FOAM_TUTORIALS/multiphaseEuler/bubbleColumnLaminar`,
matched to the OpenFOAM case so the two codes can be compared directly.

## Matched to OpenFOAM

| quantity | value | OpenFOAM source |
| --- | --- | --- |
| geometry | 0.15 m x 1.0 m | `system/blockMeshDict` vertices |
| mesh | 12 x 80 QUAD4 (960 cells, 1053 nodes) | refined from tutorial's `blocks hex ... (3 20 1)` |
| timestep | 1e-3 s, fixed | `controlDict deltaT` |
| end time | 5 s | `controlDict endTime` |
| phases | air dispersed in water | `constant/phaseProperties` |
| bubble diameter | 3 mm | `diameterModel d0` |
| drag | Schiller-Naumann | `constant/momentumTransfer` |
| lift / wall lubrication / dispersion | none | empty dicts in `momentumTransfer` |
| turbulence | laminar, both phases | `momentumTransport.*` |
| gravity | (0, -9.81) | `constant/g` |
| initial alpha.air | `0.5 + 0.49*tanh((y-0.701)/0.02)` | `0/alpha.air` codeStream |
| initial U.air | (0, 0.02) | `0/U.air` |
| inlet (bottom) | alpha 0.5/0.5, U.air (0,0.02), U.water 0 | `0/alpha.*`, `0/U.*` fixedValue |
| outlet (top) | opening, total pressure 0 | `0/p_rgh` prghTotalPressure |
| outlet velocity | zero-gradient / pressure driven | `0/U.*` pressureInletOutletVelocity |
| outlet alpha | air 1, water 0 on reversed flow | `0/alpha.*` inletOutlet |
| walls (left/right) | no-slip | `0/U.*` fixedValue (0 0 0) |

The `opening` boundary type is the correct analogue of OpenFOAM's
`prghTotalPressure` + `pressureInletOutletVelocity` + `inletOutlet` set.
A plain `outlet` would *not* match: it prescribes static, not total, pressure.

## Timestep note

The timestep is stability limited, not accuracy limited. The air phase sees a
body force `alpha*(rho_air - rho_ref)*g` against a transient coefficient of only
`alpha*rho_air`, so its initial acceleration is large until drag balances it.
The matched fixed timestep is 1e-3 s.

## Deliberate differences from OpenFOAM

The following models are switched off in OpenFOAM as well to keep the
comparison fair:

- surface tension (`sigma`)
- virtual mass (`Cvm`)
- compressible / thermal equations of state
- turbulent dispersion, lift, and wall lubrication

The main numerical difference is advection. This input deliberately uses
bounded first-order upwind; OpenFOAM uses `vanLeer` for alpha and
`limitedLinearV 1` for velocity. OpenAccel's current high-resolution CVFEM
limiter produces isolated phasic-velocity spikes in this stiff problem under
mesh refinement. Upwind remained bounded through 5 s on all three tested
meshes (3x20, 6x40, and 12x80).

## Running

```sh
../../../../../run-openaccel.sh -i input.i > log.openaccel
```

Runs 5000 timesteps to t = 5 s and writes 501 frames to `results.e`
(every 0.01 s). The log is verbose, so redirect it.

The included mesh is the 12x80 refinement used for comparison. To regenerate
it, run:

```sh
python3 make_mesh.py 12 80 mesh.e
```

## Viewing in ParaView

Open `results.e`. Available nodal fields:
`pressure`, `velocity.air`, `velocity.water`,
`volume_fraction.air`, `volume_fraction.water`.

**Pressure needs a conversion before comparing with OpenFOAM.** OpenAccel
stores the *reduced* (piezometric) pressure, `p - rho_ref*g*(y - y_ref)`, so it
is flat through the water and rises at `rho_ref*g = 9810 Pa/m` through the gas
cap. OpenFOAM's `p` is absolute and falls with height. They look opposite but
both are correct. To plot the same variable, add a Calculator filter:

```
pressure - 9810*coordsY
```

The buoyancy reference height is y = 0, which is why the expression subtracts
`9810*coordsY` and not `9810*(coordsY - 0.5)`.

Add 1e5 Pa after this conversion to compare directly with OpenFOAM's absolute
`p` rather than gauge pressure.

`volume_fraction.*` and `velocity.*` are directly comparable with no conversion.

## Comparison data

The implementation was compared with a separately generated 12x80 OpenFOAM
case using the same geometry, initial condition, physical models, timestep,
and output times. At 5 s, in the bubbly water column (`0.1 < y < 0.6`), the
mean alpha.air difference is 8.5 %, the mean vertical air-velocity difference
is 3.4 %, and the centred physical-pressure profile L2 difference is 2.7 %.
All three metrics pass the 10 % comparison threshold.

The alpha solve uses bounded first-order upwind followed by a conservative,
bounded phase-inventory projection. The projection restores phase volume where
the segregated boundedness correction removed it, while retaining the old-time
maximum principle. In the 5 s result, all 501 output frames are finite,
`alpha.air` stays in `[0.0003, 0.9880]`, and the two phase fractions sum to one
to machine precision. The whole-domain mean `alpha.air` is 0.2903, compared
with 0.3031 in OpenFOAM.

The phase-momentum correction uses a relaxation factor of 0.5. This removes a
localized interface mode from the maximum-velocity history without clipping
the velocity field. From 0.5 to 5 s, the RMS frame-to-frame change of the
maximum air velocity falls from 0.0910 to 0.0056 m/s, while retaining the
water-column velocity agreement quoted above.

`comparison/verify.py` can be used to inspect an OpenAccel result against an
OpenFOAM case. Exact agreement is not expected because OpenAccel uses nodal
CVFEM fields and bounded upwind here, while OpenFOAM uses cell-centred FVM and
higher-order limited advection.
