#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
"""
Quick (shallow) regression test runner for OpenAccel examples.

Runs a curated set of example cases with MPI using the correct 2D or 3D
binary.  The number of MPI ranks for each case is calculated from the mesh
size (number of elements) and the number of cores available on the machine.

After each run the log is parsed to report:
  - total number of iterations
  - average CPU time per iteration
  - solver-reported total wall time (plus script-measured elapsed time)
  - whether the solver reported convergence ("Converged.")

Usage example:
    python3 tools/python/regression_tests/quick_test.py \
        --exe-2d=/path/to/openaccel-2D.exe \
        --exe-3d=/path/to/openaccel-3D.exe

The script is intended to be run from the repository root, but it will work
from any directory as long as the example paths are reachable from the
script's location.
"""

import argparse
import concurrent.futures
import os
import re
import shutil
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

# Files are NetCDF classic/64-bit (CDF-1/2/5) or NetCDF-4/HDF5.
CDF_SIGNATURE = b"CDF"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------
def repo_root() -> Path:
    """Return the repository root (three levels above this script)."""
    return Path(__file__).resolve().parents[3]


def available_cores() -> int:
    """Number of usable cores on this machine."""
    try:
        # If the process is pinned to a subset of cores, respect that.
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


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


def build_plan(args: argparse.Namespace, root: Path) -> tuple[list[dict], list[str]]:
    """Inspect every selected case, warn about missing files, build a plan."""
    total_cores = args.max_total_cores or available_cores()
    plan = []
    warnings = []

    for case in SELECTED_CASES:
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

        plan.append({
            "name": name,
            "case_dir": case_dir,
            "dim": dim,
            "num_elem": num_elem,
            "nprocs": nprocs,
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


def print_report(results: list[dict]) -> None:
    """Print a report table with timing and convergence information."""
    print("\nResults report")
    print(
        "-" * 117
    )
    print(
        f"{'Case':<24s} {'Dim':>3s} {'Ranks':>5s} {'Iters':>7s} "
        f"{'Avg CPU/iter':>14s} {'Wall [s]':>11s} {'Elapsed [s]':>11s} "
        f"{'Converged':>10s} {'Status':>7s}"
    )
    print(
        "-" * 117
    )
    for r in results:
        dim = r["dim"] if r["dim"] else "N/A"
        ranks = f"{r['nprocs']}" if r["nprocs"] is not None else "N/A"
        iters = f"{r['total_iterations']}" if r["total_iterations"] is not None else "N/A"
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
        print(
            f"{r['name']:<24s} {dim:>3s} {ranks:>5s} {iters:>7s} "
            f"{avg:>14s} {wall:>11s} {elapsed:>11s} {conv:>10s} {status:>7s}"
        )
    print(
        "-" * 117
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a quick set of OpenAccel example cases in parallel (MPI).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--exe-2d",
        required=True,
        help="Path to the 2D OpenAccel binary (e.g. /path/to/openaccel-2D.exe).",
    )
    parser.add_argument(
        "--exe-3d",
        required=True,
        help="Path to the 3D OpenAccel binary (e.g. /path/to/openaccel-3D.exe).",
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

    exe_2d = Path(args.exe_2d).resolve()
    exe_3d = Path(args.exe_3d).resolve()
    for label, exe in (("--exe-2d", exe_2d), ("--exe-3d", exe_3d)):
        if not exe.is_file():
            print(f"ERROR: {label} does not exist: {exe}", file=sys.stderr)
            return 1
        if not os.access(exe, os.X_OK):
            print(f"WARNING: {label} may not be executable: {exe}", file=sys.stderr)

    try:
        mpi_runner = find_mpi_runner(args.mpi_runner)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    root = repo_root()
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    plan, warnings = build_plan(args, root)

    if warnings:
        print("\n".join(warnings), file=sys.stderr)
        print(file=sys.stderr)

    if not plan:
        print("ERROR: no cases could be prepared.", file=sys.stderr)
        return 1

    total_cores = args.max_total_cores or available_cores()
    print(f"Available cores: {total_cores}\n")
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

    print_report(results)

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    failed = [r for r in results if r["status"] == "FAIL"]
    missing = [r for r in results if r["status"] == "MISSING"]
    errors = [r for r in results if r["status"] == "ERROR"]

    print("\nSummary")
    print("-" * 70)
    print(f"Total cases:  {len(results)}")
    print(f"Passed:       {len(passed)}")
    print(f"Failed:       {len(failed)}")
    print(f"Missing:      {len(missing)}")
    print(f"Errors:       {len(errors)}")

    if failed:
        print("\nFailed cases:")
        for r in failed:
            reason = r.get("error") or f"exit code {r['returncode']}"
            print(f"  - {r['name']}: {reason}  (log: {r['log']})")
    if missing:
        print("\nMissing cases:")
        for r in missing:
            print(f"  - {r['name']}: {r.get('error')}")
    if errors:
        print("\nCases with execution errors:")
        for r in errors:
            print(f"  - {r['name']}: {r.get('error')}")

    if failed or missing or errors:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
