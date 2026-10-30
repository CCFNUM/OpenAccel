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

## Component timers (rebuild both executables after pulling)

The runner requires all component timers. Old binaries fail the first run with a
clear missing-timer error instead of producing an incomplete benchmark.

| CSV column | Measured work |
| --- | --- |
| masking_init_seconds | Entire existing masking initialization (inclusive subtotal) |
| surface_setup_seconds | Region initialization, including mesh loading, skinning and surface setup |
| graph_seconds | Background node-to-element graph construction |
| covered_classification_seconds | Existing markCovered calls, summed across bodies |
| layer_propagation_seconds | Existing computeLayers calls, including MPI work |
| reference_pairing_seconds | takeLayers: ring selection, closest-facet searches and reference pairing |
| masking_other_seconds | Remaining initialization, merging labels, field writes and collect |
| distance_seconds | Full edf_select calls, including BVH construction, searches, allocation, GPU transfers and cleanup |
| sign_seconds | Full points_inside_bvh calls, including CPU BVH construction and classification |
| distance_sign_2d_seconds | Combined 2D segment-distance and crossing loop; zero in 3D |
| field_seconds | Signed field writes and owned/shared/ghost communication |
| preparation_other_seconds | Query selection, surface/array packing, validation, merging and remaining immersed-distance work |

For one MPI rank, masking_init_seconds approximately equals its six subcomponents.
immersed_total_seconds approximately equals masking_init_seconds + distance_seconds
+ sign_seconds + distance_sign_2d_seconds + field_seconds + preparation_other_seconds.
Do not add the masking subtotal to its subcomponents again. Small clock overhead
accounts for rounding/reconciliation differences. All reductions and printing occur
after the total timer stops. Under MPI, each reported timer is the maximum over
ranks independently, so component maxima need not sum to the maximum total.

The current CUDA distance function synchronizes/copies results back before returning;
the enclosing distance timer includes completed device work. These are whole-call
timings, not separate GPU kernel/BVH-build/transfer timings. No additional distance
calculations or masking algorithm changes are introduced for instrumentation.

On Tell, inside uenv, from the repository root:

    cmake --build build-release-3D-cpu -j4
    cmake --build build-release-3D -j4

These are incremental builds of the existing configurations.
