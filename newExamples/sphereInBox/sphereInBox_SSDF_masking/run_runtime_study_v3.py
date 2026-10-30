#!/usr/bin/env python3
"""Run a stationary sphere/background refinement study. Requires PyYAML.

Run INSIDE Tell's uenv. Both cuBQL Release builds, gmsh and gmsh2exo are required.
No shell commands are constructed from user input. Existing results are preserved.
"""
import argparse
import copy
import csv
import hashlib
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
import tempfile
from datetime import datetime, timezone

try:
    import yaml
except ImportError:
    raise SystemExit("PyYAML is required: use a Python environment containing yaml.")

TIMER = re.compile(r"^SSDF_TOTAL_SECONDS=([0-9.eE+\-]+)\s*$", re.MULTILINE)
BOX = '''SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
Physical Volume("background", 1) = {1};
eps = 1e-6;
inletFaces[] = Surface In BoundingBox {-0.5-eps,-0.5-eps,-0.5-eps,-0.5+eps,0.5+eps,0.5+eps};
outletFaces[] = Surface In BoundingBox {0.5-eps,-0.5-eps,-0.5-eps,0.5+eps,0.5+eps,0.5+eps};
otherFaces[] = Boundary { Volume{1}; };
otherFaces[] -= inletFaces[];
otherFaces[] -= outletFaces[];
Physical Surface("inlet", 11) = {inletFaces[]};
Physical Surface("outlet", 12) = {outletFaces[]};
Physical Surface("symmetry", 13) = {otherFaces[]};
'''
SPHERE = '''SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.1675};
Physical Volume("immersed_body", 2) = {1};
surfaceTags[] = Boundary { Volume{1}; };
Physical Surface("sphere_wall", 21) = {surfaceTags[]};
'''

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def executable(value):
    found = shutil.which(str(value))
    if not found:
        raise ValueError(f"Executable not found or not executable: {value}")
    return str(Path(found).resolve())

def write_csv(path, rows):
    if not rows:
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)

def msh_counts(path):
    """Read counts from ASCII MSH4, without loading connectivity into memory."""
    nodes = triangles = tetrahedra = 0
    with path.open() as f:
        for line in f:
            if line.strip() == "$Nodes":
                nodes = int(next(f).split()[1])
            elif line.strip() == "$Elements":
                blocks = int(next(f).split()[0])
                for _ in range(blocks):
                    _, _, kind, count = map(int, next(f).split())
                    if kind == 2:
                        triangles += count
                    elif kind == 4:
                        tetrahedra += count
                    elif kind not in (1, 15):
                        raise ValueError(f"Unexpected non-linear/non-tet element type {kind}")
                    for _ in range(count):
                        next(f)
    if not nodes or not triangles or not tetrahedra:
        raise ValueError(f"Missing mesh counts in {path}")
    return nodes, triangles, tetrahedra

def case_input(template, box, sphere):
    data = copy.deepcopy(template)
    data["mesh"]["file_path"] = str(box)
    analysis = data["simulation"]["physical_analysis"]
    domains = analysis["domains"]
    if len(domains) != 1 or domains[0]["location"] != ["background"]:
        raise ValueError("Use a single fluid domain with location: [background].")
    regions = analysis["overrides"]["masking"]
    if len(regions) != 1 or regions[0]["shape"]["option"] != "mesh":
        raise ValueError("This study expects exactly one mesh-based masking region.")
    regions[0]["shape"]["path"] = str(sphere)
    solver = data["simulation"]["solver"]
    expert = solver["solver_control"].setdefault("expert_parameters", {})
    expert["disable_physics"] = True
    expert["force_wall_distance_calculation"] = True
    output = solver["output_control"]
    output["file_path"] = "results.e"
    output["output_frequency"] = 1
    fields = output.setdefault("output_fields", [])
    if "signed_distance" not in fields:
        fields.append("signed_distance")
    return yaml.safe_dump(data, sort_keys=False)


def checked(command, directory):
    log = directory / 'command.log'
    with log.open('w') as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        print(log.read_text(errors='replace')[-12000:], file=sys.stderr)
        raise RuntimeError(f'Command failed: {command[0]} ({result.returncode})')


def generate(args, directory, kind, h, target=None):
    geo, msh, exo = [directory / (kind + suffix) for suffix in ('.geo', '.msh', '.e')]
    for attempt in range(3):
        geo.write_text((BOX if kind == 'box' else SPHERE) +
            f'\nMesh.MeshSizeMin={h:.17g};\nMesh.MeshSizeMax={h:.17g};\n'
            'Mesh.ElementOrder=1;\nMesh.MshFileVersion=4.1;\nMesh.Binary=0;\n')
        print(f'Meshing {kind}: h={h:.6g}, target cells={target}, attempt={attempt+1}', flush=True)
        checked([args.gmsh, str(geo), '-3', '-format', 'msh4', '-nt', '1', '-o', str(msh)], directory)
        nodes, triangles, cells = msh_counts(msh)
        if target is None or abs(cells / target - 1) <= args.count_tolerance:
            break
        if attempt == 2:
            raise RuntimeError(f'Target {target}: achieved {cells} outside tolerance after 3 attempts')
        h *= (cells / target) ** (1/3)
    checked([args.converter, str(msh), str(exo)], directory)
    record = dict(h=h, nodes=nodes, triangles=triangles, cells=cells, sha256=digest(exo))
    msh.unlink()
    return exo, record


def check_build(path, gpu):
    cache = Path(path).parent / 'CMakeCache.txt'
    entries = {}
    for line in cache.read_text().splitlines():
        if ':' in line and '=' in line and not line.startswith('//'):
            key, value = line.split('=', 1)
            entries[key.split(':')[0]] = value
    enabled = lambda key: entries.get(key, '').upper() in ('ON', 'TRUE', '1', 'YES')
    if not enabled('SSDF_ENABLE_CUBIQL') or enabled('SSDF_ENABLE_CUDA') != gpu:
        raise ValueError(f'Wrong cuBQL/CUDA configuration in {cache}')
    if entries.get('CMAKE_BUILD_TYPE', '').lower() != 'release':
        raise ValueError(f'Expected Release build: {cache}')


def append_row(path, row):
    exists = path.exists()
    with path.open('a', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        if not exists:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def main():
    parser = argparse.ArgumentParser(description='CPU/GPU cuBQL masking study; only two timing CSVs retained.')
    here = Path(__file__).resolve().parent
    repo = next((p for p in here.parents if (p / '.git').exists()), here)
    parser.add_argument('--root', type=Path, default=here)
    parser.add_argument('--template', type=Path, default=here / 'input.i')
    parser.add_argument('--cpu-exe', default=str(repo / 'build-release-3D-cpu/openaccel-3D.exe'))
    parser.add_argument('--gpu-exe', default=str(repo / 'build-release-3D/openaccel-3D.exe'))
    parser.add_argument('--gmsh', default=str(Path.home() / 'HSLU/tools/meshing-env/bin/gmsh'))
    parser.add_argument('--converter', default=str(Path.home() / 'HSLU/tools/gmsh2exo/build/gmsh2exo'))
    parser.add_argument('--targets', type=int, nargs='+', default=[100000,1000000,2000000,4000000,8000000,16000000])
    parser.add_argument('--include-32m', action='store_true')
    parser.add_argument('--sphere-h', type=float, default=.004)
    parser.add_argument('--count-tolerance', type=float, default=.15)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--cpu-core', type=int, required=True, help='An allowed logical CPU; does not reserve it')
    parser.add_argument('--gpu-device', default='0')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    args.root = args.root.resolve()
    if args.include_32m:
        args.targets.append(32000000)
    args.targets = sorted(set(args.targets))
    if min(args.targets) <= 0 or args.sphere_h <= 0 or args.repeats < 1 or args.warmups < 0 or not 0 < args.count_tolerance < 1:
        parser.error('Invalid counts, sizing, repetitions or tolerance')
    for name in ('cpu_exe','gpu_exe','gmsh','converter'):
        setattr(args, name, executable(getattr(args, name)))
    executable('/usr/bin/time')
    check_build(args.cpu_exe, False)
    check_build(args.gpu_exe, True)
    if args.cpu_core not in os.sched_getaffinity(0):
        raise ValueError('Requested CPU is outside the allowed affinity set')
    template = yaml.safe_load(args.template.read_text())
    case_input(template, Path('/unused/box.e'), Path('/unused/sphere.e'))
    outputs = {mode:args.root / f'runtime_study_{mode}_v3' / 'runtime.csv' for mode in ('cpu','gpu')}
    if any(p.exists() for p in outputs.values()):
        raise ValueError('Existing timing CSV found. Use a fresh --root; no results will be overwritten.')
    print(f'Targets: {args.targets}; {len(args.targets)*args.repeats} recorded runs per backend; {args.warmups} warmup(s) per level/backend.')
    print('Both backends use cuBQL. GPU distance search; sign and masking remain on CPU.')
    if args.check_only:
        print('Preflight passed; no meshes generated or simulations started.')
        return
    os.sched_setaffinity(0, {args.cpu_core})
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close', OMP_PLACES='cores',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', SSDF_USE_BVH='1', SSDF_CELL_VALIDATE='0',
               CUDA_VISIBLE_DEVICES=args.gpu_device)
    for p in outputs.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    hashes = {mode:digest(getattr(args, mode+'_exe')) for mode in outputs}
    template_hash = digest(args.template)
    # Only this runner's temporary tree is removed; user meshes/results are never deleted.
    with tempfile.TemporaryDirectory(prefix='.runtime_v3_', dir=args.root) as scratch:
        scratch = Path(scratch)
        sphere, sr = generate(args, scratch, 'sphere', args.sphere_h)
        previous_h, previous_cells = .02, 561972
        for level,target in enumerate(args.targets,1):
            with tempfile.TemporaryDirectory(prefix='background_', dir=scratch) as work:
                work = Path(work)
                h = previous_h * (previous_cells / target)**(1/3)
                box, br = generate(args, work, 'box', h, target)
                previous_h, previous_cells = br['h'], br['cells']
                input_text = case_input(template, box, sphere)
                # Alternate backend order across repeats/levels to reduce fixed ordering bias.
                for phase,count in [('warmup',args.warmups), ('measured',args.repeats)]:
                    for repeat in range(1,count+1):
                        order = ('cpu','gpu') if (level+repeat)%2 == 0 else ('gpu','cpu')
                        for mode in order:
                            print(f'{phase}: {mode}, target={target}, actual={br["cells"]}, repeat={repeat}',flush=True)
                            with tempfile.TemporaryDirectory(prefix='run_',dir=work) as run:
                                run = Path(run)
                                (run/'input.i').write_text(input_text)
                                start = time.perf_counter()
                                with (run/'run.log').open('w') as stream:
                                    result = subprocess.run(['/usr/bin/time','-f','%M','-o',str(run/'rss.txt'),
                                        getattr(args,mode+'_exe'),'-i','input.i'],cwd=run,env=env,stdout=stream,stderr=subprocess.STDOUT)
                                wall = time.perf_counter()-start
                                log = (run/'run.log').read_text(errors='replace')
                                matches = TIMER.findall(log)
                                seconds = float(matches[0]) if len(matches)==1 else float('nan')
                                success = result.returncode==0 and math.isfinite(seconds) and seconds>0
                                rss = (run/'rss.txt').read_text().strip().splitlines()
                                row = dict(backend=mode,phase=phase,level=level,target_tetrahedra=target,
                                    background_nodes=br['nodes'],background_tetrahedra=br['cells'],box_h_m=br['h'],
                                    sphere_nodes=sr['nodes'],sphere_triangles=sr['triangles'],sphere_tetrahedra=sr['cells'],sphere_h_m=sr['h'],
                                    repeat=repeat,mpi_ranks=1,omp_threads=1,cpu_core=args.cpu_core,gpu_device=args.gpu_device if mode=='gpu' else '',
                                    immersed_total_seconds=seconds if math.isfinite(seconds) else '',process_seconds=wall,
                                    peak_process_rss_kib=rss[-1] if rss and rss[-1].isdigit() else '',
                                    status='ok' if success else 'failed',exit_code=result.returncode,
                                    timestamp_utc=datetime.now(timezone.utc).isoformat(),hostname=os.uname().nodename,
                                    load_average_1m=os.getloadavg()[0],box_sha256=br['sha256'],sphere_sha256=sr['sha256'],
                                    executable_sha256=hashes[mode],template_sha256=template_hash)
                                append_row(outputs[mode],row)
                                if not success:
                                    print(log[-12000:],file=sys.stderr)
                                    raise RuntimeError('Solver failed or timer invalid; failure recorded in CSV. Temporary outputs removed.')
    print('Finished. Only runtime_study_cpu_v3/runtime.csv and runtime_study_gpu_v3/runtime.csv retained.')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as error:
        sys.exit(f'ERROR: {error}')
