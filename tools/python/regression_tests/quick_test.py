#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
"""
Quick (shallow) regression test runner for OpenAccel examples.

Runs a curated set of example cases with MPI using the correct 2D or 3D
binary.  The number of MPI ranks for each case is calculated from the mesh
size (number of elements) and the number of cores available on the machine.

After each run the log is parsed to report:
  - analysis type (steady-state or transient), read from input.i
  - number of iterations (and, for steady-state cases, the max-iteration
    cap set in input.i, shown as performed/max)
  - average CPU time per iteration
  - solver-reported total wall time (plus script-measured elapsed time)
  - whether the solver reported convergence ("Converged.")

The report is tagged with the current git commit and can optionally be
written to a file (the report only, not the per-case solver logs).

Each run's per-case metrics are stored as a JSON record under a results
directory (one file per commit, default tools/python/regression_tests/
results/).  On subsequent runs the most recent previous record is loaded as
a baseline and the report gains a "vs base" column that reports, per case,
whether the solver converged faster or slower than the baseline (fewer
solver iterations == faster).  Use --baseline to compare against a specific
stored record instead of the most recent one.

Binaries default to the in-tree build directories relative to the repository
root:
  - 3D: build/openaccel-3D.exe
  - 2D: build2D/openaccel-2D.exe
Passing --exe-3d / --exe-2d overrides those defaults.

A single case can be forced by passing its name as a positional argument.
Any directory under examples/ is accepted, not just the curated quick set,
and the report then covers only that case.  An unknown name is an error.

Usage examples:
    # Use the in-tree build/ and build2D/ binaries
    python3 tools/python/regression_tests/quick_test.py

    # Run only one example case
    python3 tools/python/regression_tests/quick_test.py pitzDaily

    # Override binary locations and save the report to a file
    python3 tools/python/regression_tests/quick_test.py \
        --exe-2d=/path/to/openaccel-2D.exe \
        --exe-3d=/path/to/openaccel-3D.exe \
        --report-file=run_report.txt

The script locates the repository (and hence the examples/ and build
directories) from its own location, so it can be run from any directory
inside the repository.
"""

import argparse
import concurrent.futures
import json
import os
import re
import shutil
import socket
import struct
import subprocess
import sys
import threading
import time
from pathlib import Path

# Lock to keep console output from concurrent cases from interleaving mid-line.
PRINT_LOCK = threading.Lock()


# ---------------------------------------------------------------------------
# Test matrix: human readable name -> directory under the repository root
#
# Note: OpenAccel stores examples as flat subdirectories under examples/,
#       unlike accel-stk which groups them by physics category.
# ---------------------------------------------------------------------------
SELECTED_CASES = [
    {"name": "pitzDaily",     "dir": "examples/pitzDaily"},
    {"name": "T106A",         "dir": "examples/T106A"},
    {"name": "T3A",           "dir": "examples/T3A"},
    {"name": "elbow",         "dir": "examples/elbow"},
    {"name": "slab",          "dir": "examples/slab"},
    {"name": "forwardStep",   "dir": "examples/forwardStep"},
    {"name": "circularArc",   "dir": "examples/circularArc"},
    {"name": "BernardCells",  "dir": "examples/BernardCells"},
    {"name": "bump2D",        "dir": "examples/bump2D"},
    {"name": "flange",        "dir": "examples/flange"},
    {"name": "plateHole",     "dir": "examples/plateHole"},
    {"name": "TaylorCouettePoiseuille", "dir": "examples/TaylorCouettePoiseuille"},
    {"name": "cavity",        "dir": "examples/cavity"},
]

EXAMPLES_DIR = "examples"


def available_example_cases(root: Path) -> list[str]:
    """Paths of every case directory under examples/, relative to it, sorted.

    Most examples are flat (``examples/cavity``), but some are grouped into
    sub-directories (``examples/incompressible/laminar/eulerEuler/
    slipRelaxation``), so nested directories holding an ``input.i`` are
    discovered too. Flat directory names are kept unconditionally so a case
    whose input file is named differently is still listed.
    """
    examples = root / EXAMPLES_DIR
    if not examples.is_dir():
        return []
    flat = {p.name for p in examples.iterdir() if p.is_dir()}
    nested = {
        str(p.parent.relative_to(examples))
        for p in examples.rglob("input.i")
        if p.is_file()
    }
    return sorted(flat | nested)


def resolve_single_case(root: Path, requested: str) -> dict:
    """Map a user-supplied case name to a test-matrix entry.

    Any directory under examples/ is accepted, not just the curated quick
    set, with or without an 'examples/' prefix and matched case-insensitively.
    Raises ValueError if no such case exists.
    """
    name = requested.strip().strip("/")
    if name.lower().startswith(EXAMPLES_DIR + "/"):
        name = name[len(EXAMPLES_DIR) + 1:].strip("/")

    if not name:
        raise ValueError("no case name given")

    names = available_example_cases(root)
    match = next((n for n in names if n == name), None)
    if match is None:
        match = next((n for n in names if n.lower() == name.lower()), None)
    if match is None:
        # Accept the leaf name of a nested case, e.g. `slipRelaxation` for
        # `incompressible/laminar/eulerEuler/slipRelaxation`.
        leaves = [n for n in names if Path(n).name.lower() == name.lower()]
        if len(leaves) > 1:
            raise ValueError(
                f"ambiguous case name '{requested}': matches "
                + ", ".join(leaves)
            )
        if leaves:
            match = leaves[0]
    if match is None:
        raise ValueError(f"unknown case: {requested}")

    for case in SELECTED_CASES:
        if case["name"] == match:
            return dict(case)
    # Cases outside the curated quick set are still runnable. Report under the
    # leaf name so nested cases stay readable in the results table.
    return {"name": Path(match).name, "dir": f"{EXAMPLES_DIR}/{match}"}


# Files are NetCDF classic/64-bit (CDF-1/2/5) or NetCDF-4/HDF5.
CDF_SIGNATURE = b"CDF"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------
def repo_root() -> Path:
    """Return the repository root (three levels above this script).

    The path is derived from the script's own location, so the examples/ and
    build directories can be found no matter which directory the script is
    invoked from inside the repository.
    """
    return Path(__file__).resolve().parents[3]


# Default in-tree binary locations, relative to the repository root.
DEFAULT_EXE_3D_REL = Path("build") / "openaccel-3D.exe"
DEFAULT_EXE_2D_REL = Path("build2D") / "openaccel-2D.exe"


def default_exe_3d(root: Path) -> Path:
    """Default 3D binary path (build/openaccel-3D.exe under the repo root)."""
    return root / DEFAULT_EXE_3D_REL


def default_exe_2d(root: Path) -> Path:
    """Default 2D binary path (build2D/openaccel-2D.exe under the repo root)."""
    return root / DEFAULT_EXE_2D_REL


def git_commit(root: Path) -> str | None:
    """Return the current git commit (short hash, with a '-dirty' suffix).

    Returns None if git is unavailable or the repo can't be queried.
    """
    if not shutil.which("git"):
        return None
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    commit = proc.stdout.strip()
    if not commit:
        return None

    try:
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        if dirty.returncode == 0 and dirty.stdout.strip():
            commit += "-dirty"
    except Exception:
        pass
    return commit


def available_cores() -> int:
    """Number of physical cores usable for MPI ranks (respects CPU affinity).

    Open MPI counts one slot per physical core by default, not per logical
    CPU / hyperthread.  Sizing ranks from logical CPUs would request more
    ranks than there are slots and the launch would be refused ("not enough
    slots available").  Falls back to the logical CPU count if the CPU
    topology cannot be read.
    """
    logical = allowed_cpu_ids()
    cores: set[tuple[str, str]] = set()
    for cpu in logical:
        topo = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        try:
            pkg = (topo / "physical_package_id").read_text().strip()
            core = (topo / "core_id").read_text().strip()
        except OSError:
            return len(logical)
        cores.add((pkg, core))
    return len(cores) if cores else len(logical)


def allowed_cpu_ids() -> list[int]:
    """Return the list of CPU IDs this process is allowed to use."""
    try:
        return sorted(os.sched_getaffinity(0))
    except AttributeError:
        return list(range(os.cpu_count() or 1))


def find_mpi_runner(preferred: str | None) -> str:
    """Pick a sensible MPI launcher."""
    if preferred:
        return preferred
    for candidate in ("mpiexec", "mpirun"):
        if shutil.which(candidate):
            return candidate
    raise RuntimeError(
        "No MPI launcher found in PATH (tried mpiexec, mpirun). "
        "Use --mpi-runner to provide one."
    )


def parse_input_dimension(case_dir: Path) -> str | None:
    """Read the '# This is a 2D/3D case ...' statement at the top of input.i."""
    input_file = case_dir / "input.i"
    if not input_file.is_file():
        return None
    # Only the first few kilobytes matter; the statement is always at the top.
    text = input_file.read_text(errors="ignore")[:4096]
    m = re.search(r"This\s+is\s+a\s+(\d)D\s+case", text, re.IGNORECASE)
    if m:
        return f"{m.group(1)}D"
    return None


def _regex_analysis_type(text: str) -> str | None:
    """Fallback: read analysis_type.option from raw input.i text."""
    m = re.search(r"analysis_type:\s*\n\s*option:\s*(\w+)", text)
    return m.group(1) if m else None


def _regex_max_iterations(text: str) -> int | None:
    """Fallback: read the outer-solver max_iterations under convergence_controls.

    Bounded to the convergence_controls block so it never picks up the
    unrelated max_iterations under linear_solver_settings.
    """
    m = re.search(r"convergence_controls:\s*\n(.*)", text, re.DOTALL)
    if not m:
        return None
    block = m.group(1)
    stop = re.search(r"\n\s*convergence_criteria:", block)
    if stop:
        block = block[:stop.start()]
    mm = re.search(r"max_iterations:\s*(\d+)", block)
    return int(mm.group(1)) if mm else None


def parse_case_config(case_dir: Path) -> tuple[str | None, int | None]:
    """Read (analysis_type, max_iterations) from a case's input.i (YAML).

    analysis_type is 'steady_state' or 'transient'.  max_iterations is the
    outer solver iteration cap under convergence_controls (not the
    linear-solver caps).  Uses PyYAML when available and falls back to a
    targeted regex parse otherwise.
    """
    input_file = case_dir / "input.i"
    if not input_file.is_file():
        return None, None
    text = input_file.read_text(errors="ignore")

    analysis = None
    max_iters = None
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
    except Exception:
        data = None
    if isinstance(data, dict):
        try:
            analysis = (
                data["simulation"]["physical_analysis"]["analysis_type"]["option"]
            )
        except (KeyError, TypeError):
            analysis = None
        try:
            max_iters = (
                data["simulation"]["solver"]["solver_control"]
                    ["basic_settings"]["convergence_controls"]["max_iterations"]
            )
        except (KeyError, TypeError):
            max_iters = None

    if analysis is None:
        analysis = _regex_analysis_type(text)
    if max_iters is None:
        max_iters = _regex_max_iterations(text)

    analysis = str(analysis) if analysis is not None else None
    try:
        max_iters = int(max_iters) if max_iters is not None else None
    except (TypeError, ValueError):
        max_iters = None
    return analysis, max_iters


def mesh_file_for_case(case_dir: Path) -> Path | None:
    """Return the mesh file path declared in input.i (defaults to mesh.e)."""
    input_file = case_dir / "input.i"
    if not input_file.is_file():
        return None
    text = input_file.read_text(errors="ignore")
    m = re.search(r"file_path:\s*(\S+)", text)
    mesh_name = m.group(1).strip("\"'") if m else "mesh.e"
    return case_dir / mesh_name


def _read_dimension_from_cdf(mesh_path: Path, dim_name: str) -> int | None:
    """Extract a dimension length from a classic/64-bit NetCDF (CDF-1/2/5) file."""
    with open(mesh_path, "rb") as f:
        data = f.read(min(262144, os.path.getsize(mesh_path)))

    if len(data) < 4 or data[:3] != CDF_SIGNATURE:
        return None

    version = data[3]
    if version in (1, 2):
        length_bytes, fmt = 4, ">I"
    elif version == 5:
        length_bytes, fmt = 8, ">Q"
    else:
        return None

    pos = data.find(dim_name.encode())
    if pos < 0:
        return None

    name_len = len(dim_name)
    pad = (4 - name_len % 4) % 4
    val_pos = pos + name_len + pad
    if val_pos + length_bytes > len(data):
        return None

    return struct.unpack(fmt, data[val_pos:val_pos + length_bytes])[0]


def _read_dimension_from_hdf5(mesh_path: Path, dim_name: str) -> int | None:
    """
    Extract a dimension length from a NetCDF-4/HDF5 Exodus file.

    Dimensions are stored as HDF5 datasets; the dimension length is the
    8-byte little-endian integer stored immediately after the dimension name
    in the file header.
    """
    with open(mesh_path, "rb") as f:
        data = f.read(min(1048576, os.path.getsize(mesh_path)))

    if not data.startswith(HDF5_SIGNATURE):
        return None

    pos = data.find(dim_name.encode())
    if pos < 0:
        return None

    val_pos = pos + len(dim_name)
    if val_pos + 8 > len(data):
        return None

    return struct.unpack("<Q", data[val_pos:val_pos + 8])[0]


def _sane_element_count(value: int | None, mesh_path: Path, min_bytes_per_elem: int = 8) -> bool:
    """Sanity-check an element count against the mesh file size."""
    if value is None or value <= 0:
        return False
    file_size = os.path.getsize(mesh_path)
    if value > file_size // min_bytes_per_elem:
        return False
    return True


def _count_elements_with_ncdump(mesh_path: Path) -> int | None:
    """Use ncdump -h to read num_elem (works for NetCDF classic and NetCDF-4)."""
    if not shutil.which("ncdump"):
        return None
    try:
        proc = subprocess.run(
            ["ncdump", "-h", str(mesh_path)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except Exception:
        return None
    m = re.search(r"^\s*num_elem\s*=\s*(\d+)\s*;", proc.stdout, re.MULTILINE)
    if m:
        return int(m.group(1))
    return None


def _count_elements_with_h5dump(mesh_path: Path) -> int | None:
    """Use h5dump -H to read the length of the num_elem dimension scale."""
    if not shutil.which("h5dump"):
        return None
    try:
        proc = subprocess.run(
            ["h5dump", "-H", "-A", str(mesh_path)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except Exception:
        return None
    # Look for DATASET "num_elem" and the size of its first dimension.
    m = re.search(
        r'DATASET\s+"num_elem"\s*\{[^}]*DATASPACE\s+SIMPLE\s*\{\s*\(\s*(\d+)',
        proc.stdout,
        re.DOTALL,
    )
    if m:
        return int(m.group(1))
    return None


def read_dimension(mesh_path: Path, dim_name: str) -> int | None:
    """Dispatch to the correct reader based on the mesh file format."""
    if not mesh_path.is_file():
        return None

    with open(mesh_path, "rb") as f:
        header = f.read(8)

    if header[:3] == CDF_SIGNATURE:
        return _read_dimension_from_cdf(mesh_path, dim_name)
    if header.startswith(HDF5_SIGNATURE):
        return _read_dimension_from_hdf5(mesh_path, dim_name)
    return None


def count_elements(mesh_path: Path) -> int | None:
    """Return the number of elements (num_elem) in an Exodus mesh file."""
    if not mesh_path.is_file():
        return None

    with open(mesh_path, "rb") as f:
        header = f.read(8)

    # CDF classic / 64-bit files: our internal parser is reliable.
    if header[:3] == CDF_SIGNATURE:
        value = _read_dimension_from_cdf(mesh_path, "num_elem")
        if _sane_element_count(value, mesh_path):
            return value
        # If it looks bogus, try ncdump as a fallback.
        value = _count_elements_with_ncdump(mesh_path)
        if _sane_element_count(value, mesh_path):
            return value
        return None

    # NetCDF-4 / HDF5 files: external tools are much more reliable than the
    # quick header heuristic, so try them first.
    if header.startswith(HDF5_SIGNATURE):
        value = _count_elements_with_ncdump(mesh_path)
        if _sane_element_count(value, mesh_path):
            return value
        value = _count_elements_with_h5dump(mesh_path)
        if _sane_element_count(value, mesh_path):
            return value
        value = _read_dimension_from_hdf5(mesh_path, "num_elem")
        if _sane_element_count(value, mesh_path):
            return value
        return None

    return None


def detect_dimension(case_dir: Path) -> str | None:
    """
    Determine whether a case is 2D or 3D.

    First try the explicit statement in input.i.  If that is missing, fall
    back to the num_dim dimension stored in the mesh file.
    """
    dim = parse_input_dimension(case_dir)
    if dim:
        return dim

    mesh_path = mesh_file_for_case(case_dir)
    if mesh_path and mesh_path.is_file():
        nd = read_dimension(mesh_path, "num_dim")
        if nd == 2:
            return "2D"
        if nd == 3:
            return "3D"
    return None


def compute_ranks(
    num_elem: int | None,
    dim: str,
    total_cores: int,
    target_2d: int,
    target_3d: int,
    max_ranks: int | None,
) -> int:
    """Compute a sensible MPI rank count for a case."""
    if num_elem is None or num_elem <= 0:
        return 1

    target = target_2d if dim == "2D" else target_3d
    ranks = (num_elem + target - 1) // target
    ranks = max(1, min(ranks, total_cores))
    if max_ranks is not None:
        ranks = min(ranks, max_ranks)
    return ranks


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------
_ITER_RE = re.compile(r"Iter\s*=\s*(\d+)")
_CPU_TIME_RE = re.compile(r"Iter CPU Time:\s*([0-9.eE+-]+)\s*(?:\[s\])?")
_WALL_TIME_RE = re.compile(r"Total wall time:\s*([0-9.eE+-]+)\s*\[s\]")
_CONVERGED_RE = re.compile(r"\bConverged\.")


def parse_log(log_file: Path) -> dict:
    """Parse iteration count, timing and convergence from a log."""
    if not log_file.is_file():
        return {
            "total_iterations": None,
            "avg_iter_cpu_time": None,
            "total_wall_time": None,
            "converged": None,
        }

    text = log_file.read_text(errors="ignore")
    iters = _ITER_RE.findall(text)
    cpu_times = _CPU_TIME_RE.findall(text)
    wall_times = _WALL_TIME_RE.findall(text)

    total_iterations = len(iters)
    avg_iter_cpu_time = (
        sum(float(t) for t in cpu_times) / len(cpu_times)
        if cpu_times else None
    )
    # Printed once in the solver footer; take the last match to be safe.
    total_wall_time = float(wall_times[-1]) if wall_times else None
    converged = bool(_CONVERGED_RE.search(text))

    return {
        "total_iterations": total_iterations,
        "avg_iter_cpu_time": avg_iter_cpu_time,
        "total_wall_time": total_wall_time,
        "converged": converged,
    }


def stream_output(proc: subprocess.Popen, log_fp, prefix: str) -> None:
    """Tee a running process's stdout to both the log file and the console."""
    for raw in proc.stdout:
        log_fp.write(raw)
        log_fp.flush()
        line = raw.rstrip("\n")
        with PRINT_LOCK:
            if line:
                print(f"[{prefix}] {line}", flush=True)
            else:
                print(flush=True)


# ---------------------------------------------------------------------------
# Execution helpers
# ---------------------------------------------------------------------------
class CorePool:
    """Simple counting pool so concurrent cases never use more than N cores."""

    def __init__(self, cores: int):
        self._free = cores
        self._cond = threading.Condition()

    def acquire(self, n: int) -> None:
        with self._cond:
            while self._free < n:
                self._cond.wait()
            self._free -= n

    def release(self, n: int) -> None:
        with self._cond:
            self._free += n
            self._cond.notify_all()


def format_cpu_list(cores: list[int]) -> str:
    """Format a list of CPU IDs for taskset -c (supports ranges)."""
    if not cores:
        return ""
    ranges = []
    start = prev = cores[0]
    for c in cores[1:]:
        if c == prev + 1:
            prev = c
        else:
            ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = prev = c
    ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ",".join(ranges)


class CpuAllocator:
    """
    Allocates dedicated physical CPU cores to running cases.

    This prevents concurrent cases from re-using the same cores, which is
    important when solvers spawn additional threads (Kokkos, OpenMP, etc.)
    beyond the requested MPI rank count.
    """

    def __init__(self, cores: list[int]):
        self._all = sorted(set(cores))
        self._free = set(self._all)
        self._cond = threading.Condition()

    def acquire(self, n: int) -> list[int]:
        with self._cond:
            while len(self._free) < n:
                self._cond.wait()
            chosen = sorted(self._free)[:n]
            self._free.difference_update(chosen)
            return chosen

    def release(self, cores: list[int]) -> None:
        with self._cond:
            self._free.update(cores)
            self._cond.notify_all()


def run_case(
    case: dict,
    exe_2d: Path,
    exe_3d: Path,
    mpi_runner: str,
    log_dir: Path,
    args: argparse.Namespace,
    core_pool: CorePool | None = None,
    cpu_allocator: CpuAllocator | None = None,
) -> dict:
    """Run one case and return a result dictionary."""
    name = case["name"]
    case_dir = case["case_dir"]
    dim = case["dim"]
    nprocs = case["nprocs"]
    log_file = log_dir / f"{name}.log"
    allocated_cores: list[int] = []
    acquired_count = 0

    result = {
        "name": name,
        "dim": dim,
        "nprocs": nprocs,
        "analysis_type": case.get("analysis_type"),
        "max_iters": case.get("max_iters"),
        "log": str(log_file),
        "returncode": None,
        "elapsed": None,
        "error": None,
        "total_iterations": None,
        "avg_iter_cpu_time": None,
        "total_wall_time": None,
        "converged": None,
        "status": "PENDING",
    }

    if case.get("missing"):
        result["status"] = "MISSING"
        result["error"] = case.get("missing_reason")
        return result

    exe = exe_2d if dim == "2D" else exe_3d
    cmd = [mpi_runner, "-n", str(nprocs), str(exe), "-i", "input.i"]
    result["cmd"] = " ".join(cmd)

    if args.dry_run:
        print(f"[DRY-RUN] {name}: {' '.join(cmd)}")
        result["returncode"] = 0
        result["status"] = "DRY-RUN"
        return result

    try:
        if cpu_allocator is not None:
            allocated_cores = cpu_allocator.acquire(nprocs)
            cmd = ["taskset", "-c", format_cpu_list(allocated_cores)] + cmd
            if args.verbose:
                print(f"[{name}] allocated CPUs: {format_cpu_list(allocated_cores)}")
        elif core_pool is not None:
            core_pool.acquire(nprocs)
            acquired_count = nprocs

        with open(log_file, "w") as log_fp:
            log_fp.write(f"# Case: {name}\n")
            log_fp.write(f"# Dimension: {dim}\n")
            log_fp.write(f"# MPI ranks: {nprocs}\n")
            log_fp.write(f"# Command: {' '.join(cmd)}\n")
            log_fp.write(f"# Working directory: {case_dir}\n")
            log_fp.write("#\n")
            log_fp.flush()

            start = time.time()
            proc = subprocess.Popen(
                cmd,
                cwd=case_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            stream_output(proc, log_fp, name)
            proc.wait()
            elapsed = time.time() - start

        result["returncode"] = proc.returncode
        result["elapsed"] = elapsed
        result["status"] = "PASS" if proc.returncode == 0 else "FAIL"

        metrics = parse_log(log_file)
        result.update(metrics)

        print(
            f"[{result['status']:6s}] {name:24s} {dim}  ranks={nprocs:2d}  "
            f"elapsed={elapsed:6.1f}s  log={log_file}"
        )
    except Exception as exc:
        result["returncode"] = -1
        result["error"] = str(exc)
        result["status"] = "ERROR"
        print(f"[ERROR] {name}: {exc}")
    finally:
        if allocated_cores:
            cpu_allocator.release(allocated_cores)
        elif acquired_count:
            core_pool.release(acquired_count)

    return result


def build_plan(
    args: argparse.Namespace,
    root: Path,
    cases: list[dict] | None = None,
) -> tuple[list[dict], list[str]]:
    """Inspect every selected case, warn about missing files, build a plan."""
    total_cores = args.max_total_cores or available_cores()
    plan = []
    warnings = []

    for case in (cases if cases is not None else SELECTED_CASES):
        name = case["name"]
        case_dir = root / case["dir"]
        input_file = case_dir / "input.i"
        mesh_path = None
        missing = False
        reason = None

        if not case_dir.is_dir():
            reason = f"directory not found: {case_dir}"
            warnings.append(f"[WARNING] {name}: {reason}")
            missing = True
        elif not input_file.is_file():
            reason = f"input file not found: {input_file}"
            warnings.append(f"[WARNING] {name}: {reason}")
            missing = True
        else:
            mesh_path = mesh_file_for_case(case_dir)
            if not mesh_path or not mesh_path.is_file():
                reason = f"mesh file not found: {mesh_path or case_dir / 'mesh.e'}"
                warnings.append(f"[WARNING] {name}: {reason}")
                missing = True

        if missing:
            plan.append({
                "name": name,
                "case_dir": case_dir,
                "dim": None,
                "num_elem": None,
                "nprocs": 0,
                "missing": True,
                "missing_reason": reason,
            })
            continue

        dim = detect_dimension(case_dir)
        if dim is None:
            reason = "cannot determine 2D/3D dimension"
            warnings.append(f"[WARNING] {name}: {reason}")
            plan.append({
                "name": name,
                "case_dir": case_dir,
                "dim": None,
                "num_elem": None,
                "nprocs": 0,
                "missing": True,
                "missing_reason": reason,
            })
            continue

        num_elem = count_elements(mesh_path) if mesh_path and mesh_path.is_file() else None
        nprocs = compute_ranks(
            num_elem,
            dim,
            total_cores,
            args.target_per_rank_2d,
            args.target_per_rank_3d,
            args.max_ranks,
        )
        analysis_type, max_iters = parse_case_config(case_dir)

        plan.append({
            "name": name,
            "case_dir": case_dir,
            "dim": dim,
            "num_elem": num_elem,
            "nprocs": nprocs,
            "analysis_type": analysis_type,
            "max_iters": max_iters,
            "missing": False,
            "missing_reason": None,
        })

    return plan, warnings


def print_plan(plan: list[dict]) -> None:
    """Print the execution plan to stdout."""
    print("\nExecution plan")
    print("-" * 70)
    print(f"{'Case':<24s} {'Dim':>3s} {'Elements':>10s} {'Ranks':>5s}")
    print("-" * 70)
    for c in plan:
        elem_str = str(c["num_elem"]) if c["num_elem"] is not None else "unknown"
        dim_str = c["dim"] if c["dim"] else "N/A"
        ranks_str = str(c["nprocs"]) if not c.get("missing") else "N/A"
        print(f"{c['name']:<24s} {dim_str:>3s} {elem_str:>10s} {ranks_str:>5s}")
    print("-" * 70)
    print()


# ---------------------------------------------------------------------------
# Result history / baselines
#
# Each run is stored as a JSON record (one file per commit) so that a later
# run can load the most recent previous record and report, per case, whether
# the solver converged faster or slower (measured by solver iteration count).
# ---------------------------------------------------------------------------
def commit_slug(commit: str | None) -> str:
    """Turn a commit label into a filesystem-safe file stem."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", commit or "unknown")


def build_run_record(
    results: list[dict],
    commit: str | None,
    mpi_runner: str,
    total_cores: int,
) -> dict:
    """Assemble a JSON-serialisable record of a run's per-case metrics."""
    cases = []
    for r in results:
        cases.append({
            "name": r.get("name"),
            "dim": r.get("dim"),
            "analysis_type": r.get("analysis_type"),
            "nprocs": r.get("nprocs"),
            "max_iters": r.get("max_iters"),
            "total_iterations": r.get("total_iterations"),
            "avg_iter_cpu_time": r.get("avg_iter_cpu_time"),
            "total_wall_time": r.get("total_wall_time"),
            "elapsed": r.get("elapsed"),
            "converged": r.get("converged"),
            "status": r.get("status"),
        })
    return {
        "commit": commit,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host": socket.gethostname(),
        "mpi_runner": mpi_runner,
        "usable_cores": total_cores,
        "cases": cases,
    }


def save_run_record(history_dir: Path, record: dict) -> Path:
    """Write a run record to <history_dir>/<commit>.json (pretty-printed)."""
    history_dir.mkdir(parents=True, exist_ok=True)
    path = history_dir / f"{commit_slug(record.get('commit'))}.json"
    with open(path, "w") as fp:
        json.dump(record, fp, indent=2)
        fp.write("\n")
    return path


def merge_run_record(history_dir: Path, record: dict) -> dict:
    """Fold a partial run into the record already stored for the same commit.

    Without this a single-case run would drop every other case from that
    commit's record, ruining it as a baseline for the next run.
    """
    path = history_dir / f"{commit_slug(record.get('commit'))}.json"
    existing = load_record(path) if path.is_file() else None
    if not existing:
        return record

    cases = [c for c in existing.get("cases", []) if isinstance(c, dict)]
    index = {c.get("name"): i for i, c in enumerate(cases)}
    for case in record.get("cases", []):
        i = index.get(case.get("name"))
        if i is None:
            cases.append(case)
        else:
            cases[i] = case

    merged = dict(record)
    merged["cases"] = cases
    return merged


def load_record(path: Path) -> dict | None:
    """Load a stored run record, returning None if it can't be read/parsed."""
    try:
        with open(path) as fp:
            data = json.load(fp)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def find_baseline_record(
    history_dir: Path, current_commit: str | None
) -> tuple[dict | None, Path | None]:
    """Return the most recent stored record for a *different* commit.

    The current commit's own file is skipped so that re-running the same
    commit compares against the previous distinct run rather than itself.
    """
    if not history_dir.is_dir():
        return None, None
    current_stem = commit_slug(current_commit)
    candidates = [
        p for p in history_dir.glob("*.json") if p.stem != current_stem
    ]
    if not candidates:
        return None, None
    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    return load_record(latest), latest


def cases_by_name(record: dict | None) -> dict[str, dict]:
    """Index a record's case list by case name."""
    if not record:
        return {}
    return {c.get("name"): c for c in record.get("cases", []) if c.get("name")}


def compare_iterations(cur: dict, base_case: dict | None) -> str:
    """Describe how a case's iteration count moved relative to the baseline.

    Fewer solver iterations means the case converged faster.  Returns a short
    label: "new" (no baseline), "n/a" (missing counts), "same", or
    "faster -N" / "slower +N".
    """
    if base_case is None:
        return "new"
    base_iters = base_case.get("total_iterations")
    cur_iters = cur.get("total_iterations")
    if base_iters is None or cur_iters is None:
        return "n/a"
    delta = cur_iters - base_iters
    if delta == 0:
        return "same"
    return f"{'faster' if delta < 0 else 'slower'} {delta:+d}"


_TYPE_LABELS = {"steady_state": "steady", "transient": "transient"}


def _format_report_row(
    r: dict,
    base_case: dict | None = None,
    show_compare: bool = False,
) -> str:
    """Format a single results row for the report table."""
    dim = r["dim"] if r["dim"] else "N/A"
    atype = r.get("analysis_type")
    type_str = _TYPE_LABELS.get(atype, "N/A")
    ranks = f"{r['nprocs']}" if r["nprocs"] is not None else "N/A"

    total_iters = r["total_iterations"]
    max_iters = r.get("max_iters")
    if total_iters is None:
        iters = "N/A"
    elif atype == "steady_state" and max_iters is not None:
        # Outer iterations performed out of the yaml cap (steady-state only).
        iters = f"{total_iters}/{max_iters}"
    else:
        iters = f"{total_iters}"

    avg = (
        f"{r['avg_iter_cpu_time']:.3e}"
        if r["avg_iter_cpu_time"] is not None else "N/A"
    )
    wall = (
        f"{r['total_wall_time']:.2f}"
        if r.get("total_wall_time") is not None else "N/A"
    )
    elapsed = f"{r['elapsed']:.2f}" if r["elapsed"] is not None else "N/A"
    conv = (
        "yes" if r["converged"] else "no"
        if r["converged"] is not None else "N/A"
    )
    status = r.get("status", "N/A")
    row = (
        f"{r['name']:<24s} {dim:>3s} {type_str:>9s} {ranks:>5s} {iters:>11s} "
        f"{avg:>14s} {wall:>11s} {elapsed:>11s} {conv:>10s} {status:>7s}"
    )
    if show_compare:
        row += f" {compare_iterations(r, base_case):>15s}"
    return row


def format_report(
    results: list[dict],
    commit: str | None = None,
    baseline_by_name: dict[str, dict] | None = None,
    baseline_commit: str | None = None,
) -> str:
    """Build the report table with timing and convergence information.

    Steady-state cases are listed first (iterations shown as performed/max),
    then transient cases below a horizontal divider.  The current git commit,
    when known, is shown in the report title.

    When ``baseline_by_name`` is supplied, a "vs base" column is appended that
    compares each case's solver-iteration count against the baseline (fewer
    iterations == converged faster).
    """
    show_compare = baseline_by_name is not None
    header = (
        f"{'Case':<24s} {'Dim':>3s} {'Type':>9s} {'Ranks':>5s} {'Iters/Max':>11s} "
        f"{'Avg CPU/iter':>14s} {'Wall [s]':>11s} {'Elapsed [s]':>11s} "
        f"{'Converged':>10s} {'Status':>7s}"
    )
    if show_compare:
        header += f" {'vs base':>15s}"
    rule = "-" * len(header)

    steady = [r for r in results if r.get("analysis_type") == "steady_state"]
    other = [r for r in results if r.get("analysis_type") != "steady_state"]

    bits = []
    if commit:
        bits.append(f"commit {commit}")
    if show_compare and baseline_commit:
        bits.append(f"vs baseline {baseline_commit}")
    title = "Results report"
    if bits:
        title += " (" + ", ".join(bits) + ")"

    def row(r: dict) -> str:
        base_case = baseline_by_name.get(r["name"]) if show_compare else None
        return _format_report_row(r, base_case, show_compare)

    lines = [title, rule, header, rule]
    for r in steady:
        lines.append(row(r))
    if steady and other:
        lines.append(rule)
    for r in other:
        lines.append(row(r))
    lines.append(rule)
    if show_compare:
        lines.append(
            "vs base: change in solver iterations relative to the baseline "
            "commit (fewer = converged faster)."
        )
    return "\n".join(lines)


def format_summary(results: list[dict]) -> str:
    """Build the pass/fail summary block (including per-category details)."""
    passed = [r for r in results if r["status"] == "PASS"]
    failed = [r for r in results if r["status"] == "FAIL"]
    missing = [r for r in results if r["status"] == "MISSING"]
    errors = [r for r in results if r["status"] == "ERROR"]

    lines = [
        "Summary",
        "-" * 70,
        f"Total cases:  {len(results)}",
        f"Passed:       {len(passed)}",
        f"Failed:       {len(failed)}",
        f"Missing:      {len(missing)}",
        f"Errors:       {len(errors)}",
    ]

    if failed:
        lines.append("")
        lines.append("Failed cases:")
        for r in failed:
            reason = r.get("error") or f"exit code {r['returncode']}"
            lines.append(f"  - {r['name']}: {reason}  (log: {r['log']})")
    if missing:
        lines.append("")
        lines.append("Missing cases:")
        for r in missing:
            lines.append(f"  - {r['name']}: {r.get('error')}")
    if errors:
        lines.append("")
        lines.append("Cases with execution errors:")
        for r in errors:
            lines.append(f"  - {r['name']}: {r.get('error')}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a quick set of OpenAccel example cases in parallel (MPI).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "case",
        nargs="?",
        default=None,
        help=(
            "Run only this example case (a directory name under examples/, "
            "with or without the 'examples/' prefix). Any case in examples/ "
            "is accepted, not just the curated quick set. Omit to run the "
            "full quick-test matrix."
        ),
    )
    parser.add_argument(
        "--exe-2d",
        default=None,
        help=(
            "Path to the 2D OpenAccel binary. Defaults to "
            "build2D/openaccel-2D.exe under the repository root."
        ),
    )
    parser.add_argument(
        "--exe-3d",
        default=None,
        help=(
            "Path to the 3D OpenAccel binary. Defaults to "
            "build/openaccel-3D.exe under the repository root."
        ),
    )
    parser.add_argument(
        "--mpi-runner",
        default=None,
        help="MPI launcher to use (default: auto-detect mpiexec, then mpirun).",
    )
    parser.add_argument(
        "--target-per-rank-2d",
        type=int,
        default=4000,
        help="Target number of 2D elements per MPI rank.",
    )
    parser.add_argument(
        "--target-per-rank-3d",
        type=int,
        default=4000,
        help="Target number of 3D elements per MPI rank.",
    )
    parser.add_argument(
        "--max-ranks",
        type=int,
        default=None,
        help="Hard cap on the number of ranks for any individual case.",
    )
    parser.add_argument(
        "--max-total-cores",
        type=int,
        default=None,
        help="Total number of cores the script may use (default: all available).",
    )
    parser.add_argument(
        "--parallel-cases",
        action="store_true",
        help="Allow multiple cases to run concurrently while respecting --max-total-cores.",
    )
    parser.add_argument(
        "--cpu-affinity",
        action="store_true",
        help=(
            "Pin each case to its own dedicated CPU cores using taskset. "
            "This prevents concurrent cases from re-using the same physical cores, "
            "which helps when solvers spawn extra threads beyond the MPI rank count."
        ),
    )
    parser.add_argument(
        "--log-dir",
        default="quick_test_logs",
        help="Directory where per-case log files are written.",
    )
    parser.add_argument(
        "--report-file",
        nargs="?",
        const="__DEFAULT__",
        default=None,
        help=(
            "Write the run report (results table + summary, not the solver "
            "logs) to a file. Give a path to choose the location, or pass the "
            "flag with no value to write <log-dir>/report.txt."
        ),
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help=(
            "Directory where per-commit JSON result records are stored and "
            "read back as baselines "
            "(default: tools/python/regression_tests/results)."
        ),
    )
    parser.add_argument(
        "--baseline",
        default=None,
        help=(
            "Path to a stored JSON result record to compare against. "
            "Overrides the automatic 'most recent previous commit' choice."
        ),
    )
    parser.add_argument(
        "--no-compare",
        action="store_true",
        help="Do not add the 'vs base' comparison column to the report.",
    )
    parser.add_argument(
        "--no-store",
        action="store_true",
        help="Do not store this run's results in the results directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would be run without executing them.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print extra diagnostic information.",
    )

    args = parser.parse_args(argv)

    root = repo_root()

    selected_cases = None
    if args.case is not None:
        try:
            selected_cases = [resolve_single_case(root, args.case)]
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            names = available_example_cases(root)
            if names:
                print(f"Available cases under {root / EXAMPLES_DIR}:",
                      file=sys.stderr)
                for n in names:
                    print(f"  - {n}", file=sys.stderr)
            else:
                print(f"ERROR: no case directories found under "
                      f"{root / EXAMPLES_DIR}", file=sys.stderr)
            return 2

    # Only the binary a single selected case needs has to be present.
    required_dims = {"2D", "3D"}
    if selected_cases:
        single_dir = root / selected_cases[0]["dir"]
        single_dim = detect_dimension(single_dir) if single_dir.is_dir() else None
        if single_dim:
            required_dims = {single_dim}

    # Binaries default to the in-tree build directories; explicit paths override.
    exe_3d = (
        Path(args.exe_3d).resolve() if args.exe_3d else default_exe_3d(root)
    )
    exe_2d = (
        Path(args.exe_2d).resolve() if args.exe_2d else default_exe_2d(root)
    )
    for label, exe, overridden, dim_key in (
        ("--exe-3d", exe_3d, bool(args.exe_3d), "3D"),
        ("--exe-2d", exe_2d, bool(args.exe_2d), "2D"),
    ):
        if dim_key not in required_dims:
            continue
        if not exe.is_file():
            hint = "" if overridden else " (build it, or pass an explicit path)"
            print(f"ERROR: {label} does not exist: {exe}{hint}", file=sys.stderr)
            return 1
        if not os.access(exe, os.X_OK):
            print(f"WARNING: {label} may not be executable: {exe}", file=sys.stderr)

    try:
        mpi_runner = find_mpi_runner(args.mpi_runner)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    commit = git_commit(root)
    print(f"Repository root: {root}")
    print(f"Commit:          {commit or 'unknown'}")
    if "3D" in required_dims:
        print(f"3D binary:       {exe_3d}")
    if "2D" in required_dims:
        print(f"2D binary:       {exe_2d}")
    if selected_cases:
        print(f"Single case:     {selected_cases[0]['name']} "
              f"({selected_cases[0]['dir']})")

    plan, warnings = build_plan(args, root, selected_cases)

    if warnings:
        print("\n".join(warnings), file=sys.stderr)
        print(file=sys.stderr)

    if not plan:
        print("ERROR: no cases could be prepared.", file=sys.stderr)
        return 1

    total_cores = args.max_total_cores or available_cores()
    print(
        f"Usable cores: {total_cores} physical (MPI slots); "
        f"{len(allowed_cpu_ids())} logical CPUs\n"
    )
    print_plan(plan)

    cpu_allocator = None
    core_pool = None
    if args.cpu_affinity:
        if not shutil.which("taskset"):
            print(
                "WARNING: --cpu-affinity requested but taskset is not available; "
                "falling back to core-count based scheduling.",
                file=sys.stderr,
            )
            core_pool = CorePool(total_cores) if args.parallel_cases else None
        else:
            cpus = allowed_cpu_ids()[:total_cores]
            cpu_allocator = CpuAllocator(cpus)
            print(
                f"CPU affinity enabled: {len(cpus)} cores available for exclusive allocation.\n"
            )
    elif args.parallel_cases:
        core_pool = CorePool(total_cores)

    results: list[dict] = []
    if args.parallel_cases:
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(plan)) as executor:
            futures = {
                executor.submit(
                    run_case,
                    case,
                    exe_2d,
                    exe_3d,
                    mpi_runner,
                    log_dir,
                    args,
                    core_pool,
                    cpu_allocator,
                ): case["name"]
                for case in plan
            }
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    else:
        for case in plan:
            results.append(
                run_case(
                    case,
                    exe_2d,
                    exe_3d,
                    mpi_runner,
                    log_dir,
                    args,
                    core_pool,
                    cpu_allocator,
                )
            )

    # Preserve the order defined in SELECTED_CASES
    order = {c["name"]: i for i, c in enumerate(SELECTED_CASES)}
    results.sort(key=lambda r: order.get(r["name"], 999))

    # Locate a baseline record to compare against (skipped for dry runs, which
    # carry no real metrics).
    results_dir = (
        Path(args.results_dir).resolve() if args.results_dir
        else Path(__file__).resolve().parent / "results"
    )
    baseline_by_name = None
    baseline_commit = None
    if not args.dry_run and not args.no_compare:
        if args.baseline:
            baseline_path = Path(args.baseline).resolve()
            baseline_record = load_record(baseline_path)
            if baseline_record is None:
                print(f"WARNING: could not read baseline: {baseline_path}",
                      file=sys.stderr)
        else:
            baseline_record, baseline_path = find_baseline_record(
                results_dir, commit
            )
        if baseline_record:
            baseline_by_name = cases_by_name(baseline_record)
            baseline_commit = baseline_record.get("commit")
            print(f"Comparing against baseline: {baseline_path} "
                  f"(commit {baseline_commit or 'unknown'})")

    report_text = format_report(
        results, commit, baseline_by_name, baseline_commit
    )
    summary_text = format_summary(results)

    print("\n" + report_text)
    print("\n" + summary_text)

    # Store this run's metrics for future comparisons (real runs only).
    if not args.dry_run and not args.no_store:
        record = build_run_record(results, commit, mpi_runner, total_cores)
        if selected_cases:
            record = merge_run_record(results_dir, record)
        try:
            stored_path = save_run_record(results_dir, record)
            print(f"\nResults stored: {stored_path}")
        except OSError as exc:
            print(f"WARNING: could not store results in {results_dir}: {exc}",
                  file=sys.stderr)

    if args.report_file is not None:
        if args.report_file == "__DEFAULT__":
            report_path = log_dir / "report.txt"
        else:
            report_path = Path(args.report_file).resolve()
        try:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            with open(report_path, "w") as fp:
                if commit:
                    fp.write(f"# Commit: {commit}\n")
                fp.write(report_text)
                fp.write("\n\n")
                fp.write(summary_text)
                fp.write("\n")
            print(f"\nReport written to: {report_path}")
        except OSError as exc:
            print(f"WARNING: could not write report to {report_path}: {exc}",
                  file=sys.stderr)

    failed = [r for r in results if r["status"] == "FAIL"]
    missing = [r for r in results if r["status"] == "MISSING"]
    errors = [r for r in results if r["status"] == "ERROR"]

    if failed or missing or errors:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
