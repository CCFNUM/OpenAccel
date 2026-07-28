#!/usr/bin/env python3
"""Sweep a cone through angle of attack and report the lift curve as a PDF.

The cone sits apex-upstream in a box of far field. One mesh is built once and
reused for every angle: the incidence is imposed by rotating the free-stream
vector, not the body, so the sweep costs one mesh and N solves.

    tools/python/cone_aoa_sweep.py --angles 0 4 8 12 16

Each angle runs as its own case directory under the work directory, and the
force monitor on the `cone` side set gives the force vector that becomes the
lift and drag coefficients.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import time

from dataclasses import asdict, dataclass, field
from pathlib import Path


# --------------------------------------------------------------------------
# case definition
# --------------------------------------------------------------------------

@dataclass
class Case:
    """Everything that defines the sweep."""

    # geometry (metres); the cone's apex is at the origin, axis along +x
    cone_length: float = 1.0
    cone_radius: float = 0.25

    # far field box, in cone lengths
    upstream: float = 4.0
    downstream: float = 10.0
    lateral: float = 5.0

    # mesh sizes (metres)
    h_cone: float = 0.04
    h_far: float = 1.2
    refine_start: float = 0.15
    refine_end: float = 4.0

    # flow
    velocity: float = 1.0
    density: float = 1.0
    reynolds: float = 1000.0

    # solver
    max_iterations: int = 400
    residual_target: float = 1.0e-5
    ranks: int = 1

    angles: list[float] = field(default_factory=lambda: [0.0, 4.0, 8.0, 12.0, 16.0])

    @property
    def viscosity(self) -> float:
        """Dynamic viscosity that realises the requested Reynolds number."""
        return self.density * self.velocity * self.cone_length / self.reynolds

    @property
    def reference_area(self) -> float:
        """Cone base area — the reference area for the coefficients."""
        return math.pi * self.cone_radius ** 2

    @property
    def dynamic_pressure(self) -> float:
        return 0.5 * self.density * self.velocity ** 2

    @property
    def half_angle(self) -> float:
        """Semi-apex angle in degrees."""
        return math.degrees(math.atan2(self.cone_radius, self.cone_length))

    @property
    def box(self) -> tuple[float, float, float, float, float, float]:
        lc = self.cone_length
        lat = self.lateral * lc
        return (-self.upstream * lc, self.downstream * lc, -lat, lat, -lat, lat)


# --------------------------------------------------------------------------
# mesh
# --------------------------------------------------------------------------

GEO_TEMPLATE = """SetFactory("OpenCASCADE");

Cone(1) = {{0, 0, 0, {Lc}, 0, 0, 0, {Rc}, 2*Pi}};
Box(2) = {{{xmin}, {ymin}, {zmin}, {dx}, {dy}, {dz}}};
BooleanDifference(3) = {{ Volume{{2}}; Delete; }}{{ Volume{{1}}; Delete; }};

eps = {eps};
sIn[]   = Surface In BoundingBox{{ {xmin}-eps, {ymin}-eps, {zmin}-eps, {xmin}+eps, {ymax}+eps, {zmax}+eps }};
sOut[]  = Surface In BoundingBox{{ {xmax}-eps, {ymin}-eps, {zmin}-eps, {xmax}+eps, {ymax}+eps, {zmax}+eps }};
sY0[]   = Surface In BoundingBox{{ {xmin}-eps, {ymin}-eps, {zmin}-eps, {xmax}+eps, {ymin}+eps, {zmax}+eps }};
sY1[]   = Surface In BoundingBox{{ {xmin}-eps, {ymax}-eps, {zmin}-eps, {xmax}+eps, {ymax}+eps, {zmax}+eps }};
sZ0[]   = Surface In BoundingBox{{ {xmin}-eps, {ymin}-eps, {zmin}-eps, {xmax}+eps, {ymax}+eps, {zmin}+eps }};
sZ1[]   = Surface In BoundingBox{{ {xmin}-eps, {ymin}-eps, {zmax}-eps, {xmax}+eps, {ymax}+eps, {zmax}+eps }};
sCone[] = Surface In BoundingBox{{ -eps, {mRc}-eps, {mRc}-eps, {Lc}+eps, {Rc}+eps, {Rc}+eps }};

Physical Volume("fluid")     = {{3}};
Physical Surface("inlet")    = {{sIn[]}};
Physical Surface("outlet")   = {{sOut[]}};
Physical Surface("farfield") = {{sY0[], sY1[], sZ0[], sZ1[]}};
Physical Surface("cone")     = {{sCone[]}};

Field[1] = Distance;
Field[1].SurfacesList = {{sCone[]}};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = {h_cone};
Field[2].SizeMax = {h_far};
Field[2].DistMin = {refine_start};
Field[2].DistMax = {refine_end};
Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm3D = 1;
"""


def write_geo(case: Case, path: Path) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = case.box
    path.write_text(GEO_TEMPLATE.format(
        Lc=case.cone_length, Rc=case.cone_radius, mRc=-case.cone_radius,
        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
        dx=xmax - xmin, dy=ymax - ymin, dz=zmax - zmin,
        eps=1e-3 * case.cone_length,
        h_cone=case.h_cone, h_far=case.h_far,
        refine_start=case.refine_start, refine_end=case.refine_end,
    ))


def build_mesh(case: Case, workdir: Path, converter: Path) -> Path:
    """Mesh the cone once; every angle reuses the result."""
    geo = workdir / "cone.geo"
    msh = workdir / "cone.msh"
    exo = workdir / "mesh.e"

    if exo.exists():
        print(f"  reusing existing mesh {exo}")
        return exo

    write_geo(case, geo)
    print(f"  gmsh -> {msh.name}")
    run(["gmsh", str(geo), "-3", "-format", "msh", "-o", str(msh), "-v", "2"],
        cwd=workdir, what="gmsh")

    print(f"  converting -> {exo.name}")
    mod = load_module(converter, "gmsh_to_openaccel_exodus")
    data = mod.build_mesh(str(msh), "fluid")
    mod.write_exodus(data, str(exo))

    cells = sum(len(b["conn"]) for b in data["blocks"])
    sets = {s["name"]: len(s["elem"]) for s in data["sidesets"]}
    print(f"  mesh: {len(data['nodes'])} nodes, {cells} cells, side sets {sets}")

    missing = {"cone", "inlet", "outlet", "farfield"} - set(sets)
    if missing:
        raise SystemExit(f"mesh is missing expected side sets: {sorted(missing)}")
    return exo


# --------------------------------------------------------------------------
# solver input
# --------------------------------------------------------------------------

INPUT_TEMPLATE = """# Cone at {alpha:g} deg angle of attack -- generated by cone_aoa_sweep.py
mesh:
    file_path: mesh.e{decomposition}
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
        - name: default_domain
          location: [fluid]
          materials: [fluid]
          type: fluid
          domain_models:
            reference_pressure: 0.0
          fluid_models:
            turbulence:
                option: laminar
          boundaries:
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: {u:.10g}
                    v: {v:.10g}
                    w: 0.0
          # Free-stream velocity on the lateral boundaries. An `opening` would
          # be the more permissive choice, but its pressure is a TOTAL pressure
          # and the resulting system diverges here; a velocity far field is
          # stable and, with the box this far out, barely constrains the body.
          - name: farfield
            type: inlet
            location: [farfield]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: {u:.10g}
                    v: {v:.10g}
                    w: 0.0
          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
                mass_and_momentum:
                    option: static_pressure
                    relative_pressure: 0.0
          - name: cone
            type: wall
            location: [cone]
          initialization:
            velocity:
                option: value
                velocity: [{u:.10g}, {v:.10g}, 0.0]
            pressure:
                option: value
                pressure: 0.0
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                convergence_controls:
                    min_iterations: 5
                    max_iterations: {max_iterations}
                    physical_timescale: {timescale:.10g}
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: {residual_target:g}
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 30
                        rtol: 1.0e-2
                        atol: 1.0e-14
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
                    pressure_correction:
                        family: HYPRE
                        min_iterations: 3
                        max_iterations: 30
                        rtol: 1.0e-3
                        atol: 1.0e-14
                        options:
                            solver:
                                type: GMRES
                            precond:
                                type: BoomerAMG
                                coarsen_type: 10
                                interp_type: 6
                                relax_type: 18
                                strong_threshold: 0.25
                                num_sweeps: 1
                                max_levels: 20
                                aggressive_levels: 1
                                trunc_factor: 0.3
        output_control:
            file_path: results.e
            output_frequency: {max_iterations}
            output_fields: [velocity, pressure]
            corrected_boundary_values: true
            post_process:
            - name: cone_forces
              type: force
              options:
                calculate_moment: false
              location: [cone]
              frequency: 1
              write_to_file: true
    material_library:
    - name: fluid
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: {density:.10g}
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: {viscosity:.10g}
"""


def write_input(case: Case, alpha: float, path: Path) -> None:
    rad = math.radians(alpha)
    ex, ey = math.cos(rad), math.sin(rad)
    decomposition = ("\n    automatic_decomposition_type: rcb"
                     if case.ranks > 1 else "")
    path.write_text(INPUT_TEMPLATE.format(
        alpha=alpha,
        u=case.velocity * ex, v=case.velocity * ey,
        decomposition=decomposition,
        max_iterations=case.max_iterations,
        residual_target=case.residual_target,
        timescale=case.cone_length / case.velocity,
        density=case.density, viscosity=case.viscosity,
    ))


# --------------------------------------------------------------------------
# running and reading back
# --------------------------------------------------------------------------

@dataclass
class Result:
    alpha: float
    fx: float
    fy: float
    fz: float
    iterations: int
    converged: bool
    seconds: float

    def lift(self, alpha_rad: float) -> float:
        """Force component normal to the free stream, in the x-y plane."""
        return -self.fx * math.sin(alpha_rad) + self.fy * math.cos(alpha_rad)

    def drag(self, alpha_rad: float) -> float:
        return self.fx * math.cos(alpha_rad) + self.fy * math.sin(alpha_rad)


def read_forces(case_dir: Path) -> tuple[float, float, float, int]:
    """Last row of the force monitor: iter, force_x, force_y, force_z."""
    monitor = case_dir / "postProcessing" / "cone_forces"
    if not monitor.exists():
        raise FileNotFoundError(f"no force monitor at {monitor}")

    rows = [ln.split() for ln in monitor.read_text().splitlines() if ln.strip()]
    if len(rows) < 2:
        raise ValueError(f"force monitor {monitor} has no data rows")

    header, last = rows[0], rows[-1]
    idx = {name: i for i, name in enumerate(header)}
    for key in ("force_x", "force_y", "force_z"):
        if key not in idx:
            raise ValueError(f"force monitor is missing '{key}'; header={header}")
    return (float(last[idx["force_x"]]), float(last[idx["force_y"]]),
            float(last[idx["force_z"]]), int(float(last[idx[header[0]]])))


def run_angle(case: Case, alpha: float, workdir: Path, mesh: Path,
              solver: Path) -> Result:
    case_dir = workdir / f"aoa_{alpha:g}".replace(".", "p")
    case_dir.mkdir(parents=True, exist_ok=True)

    # A copy, not a symlink: the solver may write decomposition siblings.
    local_mesh = case_dir / "mesh.e"
    if not local_mesh.exists():
        shutil.copy2(mesh, local_mesh)

    write_input(case, alpha, case_dir / "input.i")
    shutil.rmtree(case_dir / "postProcessing", ignore_errors=True)

    cmd = [str(solver), "-i", "input.i"]
    if case.ranks > 1:
        cmd = ["mpirun", "-np", str(case.ranks), "--oversubscribe"] + cmd

    started = time.time()
    log = case_dir / "solver.log"
    with log.open("w") as handle:
        proc = subprocess.run(cmd, cwd=case_dir, stdout=handle,
                              stderr=subprocess.STDOUT)
    elapsed = time.time() - started

    if proc.returncode != 0:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-25:])
        raise SystemExit(
            f"solver failed for alpha={alpha} (rc={proc.returncode})\n"
            f"--- tail of {log} ---\n{tail}")

    fx, fy, fz, iterations = read_forces(case_dir)
    result = Result(alpha, fx, fy, fz, iterations,
                    iterations < case.max_iterations, elapsed)
    (case_dir / "summary.json").write_text(json.dumps(asdict(result), indent=2))
    return result


def load_result(case_dir: Path, alpha: float, max_iterations: int) -> Result:
    """Re-read a solved case so the report can be redrawn without re-solving."""
    cached = case_dir / "summary.json"
    if cached.exists():
        return Result(**json.loads(cached.read_text()))
    fx, fy, fz, iterations = read_forces(case_dir)
    return Result(alpha, fx, fy, fz, iterations,
                  iterations < max_iterations, float("nan"))


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def write_report(case: Case, results: list[Result], path: Path,
                 mesh_info: str) -> float | None:
    """Write the PDF and return the small-angle lift-curve slope, if fitted."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages

    alphas = np.array([r.alpha for r in results])
    rads = np.radians(alphas)
    q_s = case.dynamic_pressure * case.reference_area
    cl = np.array([r.lift(a) for r, a in zip(results, rads)]) / q_s
    cd = np.array([r.drag(a) for r, a in zip(results, rads)]) / q_s

    # Lift-curve slope from the small-angle end (<= 8 deg), if we have it.
    # Fitted through the origin: a cone is axisymmetric, so C_L(0) = 0.
    linear = alphas <= 8.0
    slope = None
    if linear.sum() >= 2 and np.any(rads[linear] != 0):
        slope = float(np.dot(rads[linear], cl[linear])
                      / np.dot(rads[linear], rads[linear]))

    # C_L at zero incidence must vanish by symmetry; whatever is left is noise.
    zero = np.flatnonzero(alphas == 0.0)
    zero_offset = float(cl[zero[0]]) if zero.size else None

    with PdfPages(path) as pdf:
        # ---- page 1: the lift curve -------------------------------------
        fig, (ax, txt) = plt.subplots(
            2, 1, figsize=(8.27, 11.69), gridspec_kw={"height_ratios": [3, 2]})
        fig.suptitle("Cone lift curve — OpenAccel", fontsize=15, y=0.97)

        ax.plot(alphas, cl, "o-", color="#1f77b4", lw=2, ms=6, label="$C_L$ (CFD)")
        if slope is not None:
            ax.plot(alphas, slope * rads, "--", color="#888888", lw=1.2,
                    label=rf"linear fit, $dC_L/d\alpha$ = {slope:.3f} /rad")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xlabel(r"angle of attack $\alpha$  [deg]")
        ax.set_ylabel(r"lift coefficient  $C_L$")
        ax.set_title(f"Re = {case.reynolds:g} (on cone length), laminar")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", frameon=False)

        for a, c in zip(alphas, cl):
            ax.annotate(f"{c:.3f}", (a, c), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7, color="#555555")

        txt.axis("off")
        rows = [["alpha [deg]", "C_L", "C_D", "L/D", "iters", "converged"]]
        for r, a, l, d in zip(results, alphas, cl, cd):
            rows.append([f"{a:g}", f"{l:+.4f}", f"{d:.4f}",
                         f"{l / d:+.3f}" if d else "-",
                         str(r.iterations), "yes" if r.converged else "NO"])
        table = txt.table(cellText=rows[1:], colLabels=rows[0],
                          loc="upper center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.4)
        pdf.savefig(fig)
        plt.close(fig)

        # ---- page 2: drag, polar, and the case definition ---------------
        fig, axes = plt.subplots(2, 2, figsize=(8.27, 11.69),
                                 gridspec_kw={"height_ratios": [1, 1.7]})
        fig.suptitle("Drag, polar, and case definition", fontsize=14, y=0.97)

        axes[0][0].plot(alphas, cd, "s-", color="#d62728", lw=2, ms=5)
        axes[0][0].set_xlabel(r"$\alpha$  [deg]")
        axes[0][0].set_ylabel(r"$C_D$")
        axes[0][0].set_title("Drag coefficient")
        axes[0][0].grid(alpha=0.3)

        axes[0][1].plot(cd, cl, "o-", color="#2ca02c", lw=2, ms=5)
        axes[0][1].set_xlabel(r"$C_D$")
        axes[0][1].set_ylabel(r"$C_L$")
        axes[0][1].set_title("Drag polar")
        axes[0][1].grid(alpha=0.3)

        for cell in (axes[1][0], axes[1][1]):
            cell.axis("off")

        xmin, xmax, ymin, ymax, zmin, zmax = case.box
        summary = [
            "GEOMETRY",
            f"  cone length          {case.cone_length:g} m",
            f"  base radius          {case.cone_radius:g} m",
            f"  semi-apex angle      {case.half_angle:.2f} deg",
            f"  reference area       {case.reference_area:.5g} m^2 (base area)",
            "",
            "DOMAIN",
            f"  x  {xmin:g} .. {xmax:g} m",
            f"  y  {ymin:g} .. {ymax:g} m",
            f"  z  {zmin:g} .. {zmax:g} m",
            "",
            "MESH",
            f"  {mesh_info}",
            f"  size on cone / far   {case.h_cone:g} / {case.h_far:g} m",
            "",
            "FLOW",
            f"  free-stream speed    {case.velocity:g} m/s",
            f"  density              {case.density:g} kg/m^3",
            f"  dynamic viscosity    {case.viscosity:.5g} Pa s",
            f"  Reynolds number      {case.reynolds:g}",
            f"  dynamic pressure     {case.dynamic_pressure:.5g} Pa",
            "",
            "SOLVER",
            "  steady state, laminar, high-resolution advection",
            f"  max iterations       {case.max_iterations}",
            f"  RMS residual target  {case.residual_target:g}",
            f"  MPI ranks            {case.ranks}",
            "",
            "METHOD",
            "  One mesh; incidence imposed by rotating the free-stream",
            "  vector, so the body never moves. Forces come from the",
            "  solver's force monitor on the 'cone' side set;",
            "  L = -Fx sin(a) + Fy cos(a), D = Fx cos(a) + Fy sin(a).",
            "",
            "CAVEATS",
            "  Lateral boundaries carry the free-stream velocity, so at",
            "  incidence they mildly constrain the cross-flow; the box is",
            "  sized to keep that small. Tetrahedra only -- no prism layers",
            "  -- so the boundary layer is resolved by size alone, which is",
            "  why the Reynolds number is kept low.",
        ]
        if zero_offset is not None:
            summary += [
                "",
                "  A cone is axisymmetric, so C_L(0) is exactly zero. This run",
                f"  gives {zero_offset:+.4f}, which is the asymmetry of the",
                "  unstructured mesh and sets the noise floor on every point.",
            ]
        axes[1][0].text(0, 1, "\n".join(summary), va="top", ha="left",
                        family="monospace", fontsize=7, transform=axes[1][0].transAxes)

        raw = ["RAW FORCES  [N]", "", f"{'alpha':>7} {'Fx':>13} {'Fy':>13} {'Fz':>13}"]
        for r in results:
            raw.append(f"{r.alpha:7g} {r.fx:13.5e} {r.fy:13.5e} {r.fz:13.5e}")
        total = sum(r.seconds for r in results)
        raw += ["", f"total solve time  {total:.0f} s" if math.isfinite(total)
                else "total solve time  (not recorded)"]
        axes[1][1].text(0, 1, "\n".join(raw), va="top", ha="left",
                        family="monospace", fontsize=8, transform=axes[1][1].transAxes)

        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)

    return slope


# --------------------------------------------------------------------------
# plumbing
# --------------------------------------------------------------------------

def run(cmd: list[str], cwd: Path, what: str) -> None:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        raise SystemExit(f"{what} failed (rc={proc.returncode})")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def find_solver(explicit: str | None, repo: Path) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get("OPENACCEL_BIN_DIR"):
        candidates.append(Path(os.environ["OPENACCEL_BIN_DIR"]) / "openaccel-3D.exe")
    candidates.append(repo / "build" / "openaccel-3D.exe")

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise SystemExit(
        "cannot find openaccel-3D.exe; pass --solver or set OPENACCEL_BIN_DIR.\n"
        "tried:\n  " + "\n  ".join(str(c) for c in candidates))


def find_converter(explicit: str | None, repo: Path) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    # the suite keeps the pre-processor beside the solver
    candidates.append(repo.parent / "OpenAccel-Pre" / "tools" / "gmsh_to_openaccel_exodus.py")

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise SystemExit(
        "cannot find gmsh_to_openaccel_exodus.py; pass --converter.\n"
        "tried:\n  " + "\n  ".join(str(c) for c in candidates))


def main() -> int:
    repo = Path(__file__).resolve().parents[2]

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--angles", type=float, nargs="+",
                   default=[0, 4, 8, 12, 16], metavar="DEG",
                   help="angles of attack in degrees (default: 0 4 8 12 16)")
    p.add_argument("--workdir", default=str(repo / "runs" / "cone_aoa_sweep"),
                   help="where meshes, cases and the report land")
    p.add_argument("--report", default=None,
                   help="PDF path (default: <workdir>/cone_lift_curve.pdf)")
    p.add_argument("--solver", default=None, help="path to openaccel-3D.exe")
    p.add_argument("--converter", default=None,
                   help="path to gmsh_to_openaccel_exodus.py")

    p.add_argument("--cone-length", type=float, default=1.0)
    p.add_argument("--cone-radius", type=float, default=0.25)
    p.add_argument("--velocity", type=float, default=1.0)
    p.add_argument("--density", type=float, default=1.0)
    p.add_argument("--reynolds", type=float, default=1000.0,
                   help="Reynolds number on cone length; sets the viscosity")
    p.add_argument("--h-cone", type=float, default=0.04,
                   help="target cell size on the cone surface")
    p.add_argument("--h-far", type=float, default=1.2,
                   help="target cell size at the far field")
    p.add_argument("--max-iterations", type=int, default=400)
    p.add_argument("--residual-target", type=float, default=1e-5)
    p.add_argument("--np", type=int, default=1, dest="ranks",
                   help="MPI ranks per solve (needs a parallel-netCDF Trilinos)")
    p.add_argument("--mesh-only", action="store_true",
                   help="build the mesh and stop")
    p.add_argument("--report-only", action="store_true",
                   help="redraw the PDF from cases already solved in --workdir")

    args = p.parse_args()

    case = Case(
        cone_length=args.cone_length, cone_radius=args.cone_radius,
        h_cone=args.h_cone, h_far=args.h_far,
        velocity=args.velocity, density=args.density, reynolds=args.reynolds,
        max_iterations=args.max_iterations,
        residual_target=args.residual_target,
        ranks=args.ranks, angles=list(args.angles),
    )

    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    if args.report_only:
        return report_only(case, workdir, args.report)

    solver = find_solver(args.solver, repo)
    converter = find_converter(args.converter, repo)

    print(f"work directory : {workdir}")
    print(f"solver         : {solver}")
    print(f"cone           : L={case.cone_length:g} m, R={case.cone_radius:g} m, "
          f"semi-apex {case.half_angle:.1f} deg")
    print(f"flow           : U={case.velocity:g} m/s, rho={case.density:g}, "
          f"mu={case.viscosity:.4g}, Re={case.reynolds:g}")
    print(f"angles         : {', '.join(f'{a:g}' for a in case.angles)} deg\n")

    print("[1/3] mesh")
    mesh = build_mesh(case, workdir, converter)
    mesh_info = mesh_summary(mesh)
    if args.mesh_only:
        print("\nmesh only — stopping here")
        return 0

    print(f"\n[2/3] solving {len(case.angles)} angles")
    results: list[Result] = []
    for i, alpha in enumerate(case.angles, start=1):
        print(f"  ({i}/{len(case.angles)}) alpha = {alpha:g} deg ...",
              end="", flush=True)
        result = run_angle(case, alpha, workdir, mesh, solver)
        results.append(result)
        rad = math.radians(alpha)
        cl = result.lift(rad) / (case.dynamic_pressure * case.reference_area)
        flag = "" if result.converged else "  [NOT CONVERGED]"
        print(f" CL = {cl:+.4f}  ({result.iterations} iters, "
              f"{result.seconds:.0f}s){flag}")

    report = Path(args.report) if args.report else workdir / "cone_lift_curve.pdf"
    print(f"\n[3/3] report -> {report}")
    slope = write_report(case, results, report, mesh_info)
    if slope is not None:
        print(f"  lift-curve slope (alpha <= 8 deg): {slope:.4f} per rad "
              f"({math.radians(1) * slope:.5f} per deg)")

    if any(not r.converged for r in results):
        bad = ", ".join(f"{r.alpha:g}" for r in results if not r.converged)
        print(f"\nWARNING: hit the iteration cap at alpha = {bad} deg; "
              f"raise --max-iterations before trusting those points.")
    return 0


def report_only(case: Case, workdir: Path, report_arg: str | None) -> int:
    """Redraw the PDF from case directories that were solved earlier."""
    results = []
    for alpha in case.angles:
        case_dir = workdir / f"aoa_{alpha:g}".replace(".", "p")
        if not case_dir.is_dir():
            raise SystemExit(f"no solved case at {case_dir}")
        results.append(load_result(case_dir, alpha, case.max_iterations))

    report = Path(report_arg) if report_arg else workdir / "cone_lift_curve.pdf"
    slope = write_report(case, results, report, mesh_summary(workdir / "mesh.e"))
    print(f"report -> {report}")
    if slope is not None:
        print(f"  lift-curve slope (alpha <= 8 deg): {slope:.4f} per rad")
    return 0


def mesh_summary(mesh: Path) -> str:
    try:
        import netCDF4 as nc
        with nc.Dataset(mesh) as ds:
            nodes = ds.dimensions["num_nodes"].size
            cells = ds.dimensions["num_elem"].size
        return f"{nodes} nodes, {cells} tetrahedra"
    except Exception:
        return "(mesh statistics unavailable)"


if __name__ == "__main__":
    sys.exit(main())
