# CPU/GPU masking runtime study v3

Run `run_runtime_study_v3.py` inside Tell's uenv, using the meshing-env Python.
Requires PyYAML, gmsh, gmsh2exo, GNU time, and both built 3D Release executables.
Defaults locate tools under ~/HSLU/tools and builds relative to this repository.
Override --gmsh, --converter, --cpu-exe and --gpu-exe for other installations.

Geometry: unit cube centered at zero; sphere radius 0.1675 m centered at zero.
Fixed sphere mesh size 0.004 m. Background targets: 0.1, 1, 2, 4, 8, 16 million
linear tetrahedra, with +/-15% target tolerance and at most three meshing attempts.
Actual node, triangle and tetrahedron counts are recorded, never inferred from targets.
Use --include-32m only after assessing the 16-million case's time and memory.
The 32-million flag adds that level to a fresh study; it is not a resume option.

Both builds must enable SSDF_ENABLE_CUBIQL; CPU must disable SSDF_ENABLE_CUDA,
GPU must enable it. SSDF_USE_BVH=1 is explicit for both. The GPU path accelerates
distance BVH construction/query; classification and original masking remain CPU.
Build flags are checked from CMakeCache; rebuild stale binaries before running.

From this directory, inside uenv (do not nest uenv run):

    ~/HSLU/tools/meshing-env/bin/python run_runtime_study_v3.py --cpu-core 24 --check-only
    ~/HSLU/tools/meshing-env/bin/python -u run_runtime_study_v3.py --cpu-core 24 --gpu-device 0

Choose a currently idle core/GPU; 24 and 0 are examples, not reserved resources.
Defaults: 3 measured repetitions and 1 warmup per backend per background.
That is 18 measured + 6 warmup runs per backend. Warmup rows have phase=warmup;
filter phase=measured and status=ok for comparisons. Backend order alternates.

Only runtime_study_cpu_v3/runtime.csv and runtime_study_gpu_v3/runtime.csv remain
as result artifacts. CSVs are flushed after each run. Existing CSVs are never
overwritten; use a fresh --root directory for another study. Input/scripts remain.
One temporary sphere and one background at a time are shared across both backends.
Temporary meshes, input copies, solver field outputs and logs are removed on normal
exit, failure or handled interruption. On solver failure the log tail is printed,
a failed CSV row is saved, and execution stops. SIGKILL/power loss may leave scratch.
Mesh-generation failures print diagnostics and stop before recording a simulation.

immersed_total_seconds is SSDF_TOTAL_SECONDS from OpenAccel: masking initialization
plus immersed distance/sign calculation. It excludes external meshing and subsequent
ordinary wall-distance calculation/output. process_seconds includes whole solver
execution and output I/O. peak_process_rss_kib is GNU time maximum solver RSS,
not GPU memory or peak mesh-generation memory. No flow equations are solved.
Hashes, timestamps, CPU affinity and host load accompany the measurements.
Timing-only study: no numerical accuracy validation is performed by this runner.

To reproduce on Tell, pull SSDF and update submodules, then use newExamples paths.
Existing Tell directories under examples are not moved automatically by git because
their generated meshes/results are untracked. The v3 runner creates its own meshes.
