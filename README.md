# OpenAccel

**OpenAccel** is a CPU-parallel, vertex-based finite volume (**CVFEM**) solver
built on the **Trilinos-STK** mesh infrastructure. It employs a
**pressure-based, segregated approach** to solve the governing equations, making
it well-suited for incompressible and low-Mach compressible flows. The solver
addresses a broad range of physics, including fluid flow, heat transfer,
turbulence, multiphase free-surface flows and solid mechanics--while providing a
framework for multiphysics applications, including Arbitrary Lagrangian-Eulerian
(ALE) methods and Conjugate Heat Transfer (CHT).

> **Note:** This project is under active development. API, features, and
> documentation may change without notice.

---


## Capabilities

* **Incompressible and compressible single-phase flows**
* **Incompressible multiphase flows** using the Volume of Fluid (VOF) method
* **Fluid-Structure Interaction (FSI)** via a partitioned ALE approach
* **Heat Transfer** including Conjugate Heat Transfer (CHT)
* **Turbulence Modeling** with support for RANS (e.g., k-ω SST and transition
  SST)


## Features

* **Modern Foundation:** Built using the **C++20** standard.
* **High Performance:** Optimized for massive parallelism via Trilinos-STK.


## Documentation

A **Theory Guide** covering the mathematical foundations of OpenAccel is
available in [`docs/theory/`](./docs/theory/) and bundled with releases.

The guide is automatically compiled to PDF on every commit via GitHub Actions.
To get the latest PDF:

1. Go to the **Actions** tab of this repository.
2. Open the most recent **Build Theory Guide** workflow run.
3. Download the `theory-guide-pdf` artifact.

Topics covered include governing equations, turbulence models (k-ε, k-ω SST,
Transition SST), heat transfer, mesh deformation (ALE/displacement diffusion),
free surface flow (VoF + FCT/cMULES), and discretisation schemes.

---


## Dependencies

### Required

* **Trilinos:** Must be built with **STK**. **Tpetra** and **Belos** are
  optionally enabled.
  > **⚠️ Note on Sierra migration:**
  > OpenAccel currently depends on a Trilinos build with
  > `-DSIERRA_MIGRATION:BOOL=ON` set.  This dependency will be removed in the
  > future.  A convenient way to install dependencies is to use
  > [Spack](https://spack.io/) which is further explained below.

  > **⚠️ Note on parallel netCDF (manual Trilinos builds):**
  > When building Trilinos manually, it is **essential** that it is linked
  > against the **PnetCDF-enabled netCDF**.  This is what allows SEACAS/Ioss to
  > read a serial Exodus mesh collectively, and therefore what enables the
  > **runtime automatic decomposition** used by STK (see [Parallel Execution
  > with MPI](#parallel-execution-with-mpi)).  Without it, meshes must be
  > decomposed offline with `decomp` instead.
  >
  > On Debian/Ubuntu the required flavour is provided by
  > `libnetcdf-pnetcdf-dev`:
  > ```sh
  > sudo apt-get install libnetcdf-pnetcdf-dev
  > ```
  > This installs a *parallel* netCDF **side by side** with the serial
  > `libnetcdf-dev`, under a separate prefix.  Trilinos' package prefind will
  > otherwise pick up the serial build, so the parallel one has to be pointed
  > to explicitly in the Trilinos configuration:
  > ```sh
  > -DNetcdf_ALLOW_PACKAGE_PREFIND:BOOL=OFF \
  > -DTPL_Netcdf_INCLUDE_DIRS:PATH=/usr/lib/x86_64-linux-gnu/netcdf/pnetcdf/include \
  > -DTPL_Netcdf_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/netcdf/pnetcdf/libnetcdf.so \
  > -DTPL_Netcdf_PARALLEL:BOOL=ON \
  > ```
  > Adjust the paths to match your distribution's multiarch triplet.  The
  > provided Spack environments already build a parallel netCDF, so this only
  > concerns manual builds.
* **YAML-cpp:** Required for input parsing and configuration.
* **MPI:** Essential for parallel execution and distributed memory
  communication.  [MPICH](https://www.mpich.org/) and
  [OpenMPI](https://www.open-mpi.org/) implementations are tested.


### Optional (Backend Solvers)

* **PETSc** (>= 3.18): Optional high-performance solver backend. (Some examples
  are configured with PETSc solvers. See
  [`examples/cavity/input.i`](./examples/cavity/input.i) and
  [`examples/elbow/input.i`](./examples/elbow/input.i).)
* **HYPRE** (>= 3.0): Optional multigrid preconditioner support. (See
  [`examples/airfoil/input.i`](./examples/airfoil/input.i) and
  [`examples/pitzDaily/input.i`](./examples/pitzDaily/input.i).)
* **gnuplot:** Optional live residual plotting.


### Submodules

The following libraries are provided as Git submodules and bundled into
releases:
* **[Eigen](https://gitlab.com/libeigen/eigen):** Template library for linear
  algebra.
* **[ExprTk](https://github.com/ArashPartow/exprtk):** Mathematical expression
  parsing and evaluation.
* **[nanoflann](https://github.com/jlblancoc/nanoflann):** KD-tree library for
  nearest-neighbour searches.
* **[gplotpp](https://github.com/ziotom78/gplotpp):** C++ interface for gnuplot.
* **[liblinsolve](https://github.com/CCFNUM/LibLinSolve.git):** Wrapper library
  for HYPRE, PETSc, and Trilinos linear solver utilities.

---

## Installation

### 0. Dependencies

The required and optional dependencies must be installed either manually or via
the system package manager.  The recommended way is to use the provided
[Spack](https://spack.io/) `mpich.yaml` or `openmpi.yaml` environments located
in [`./tools/spack`](./tools/spack).  See the Spack [getting
started](https://spack.readthedocs.io/en/latest/getting_started.html#getting-started)
section for setting up Spack on your system.  In order to install the correct
dependencies, the `ccfnum.openaccel` Spack repository must be added first
```sh
spack repo add ./tools/spack/repos/spack_repo/ccfnum/openaccel
```
Now you can add the desired environment (or both), for example using MPICH
```sh
spack env create mpich ./tools/spack/mpich.yaml
```
In order to use the environment it must be loaded using
```sh
spack env activate mpich
```
Before you can use the software stack after loading the environment for the
first time, the packages must be installed.  This can be done with
```sh
spack install [-j<ncores>]
```
As Spack is building the required dependencies from source, installation may
take a while.  Parallelism can be exploited by specifying `ncores` via the `-j`
option.


### 1. Cloning the Repository

Retrieve the source code and initialize submodules:

```bash
git clone https://github.com/CCFNUM/OpenAccel.git
cd OpenAccel
git submodule update --init --recursive
```


### 2. Configuration with CMake

To build OpenAccel, ensure the required dependencies are correctly installed or
activate one of the provided Spack environments.
Assuming one of the provided Spack environments is active, the code can built
using CMake with a Ninja backend:
```bash
cmake . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=Release
```
If dependencies are installed manually and CMake cannot find them automatically,
the install directories can be specified using `-DTrilinos_DIR=`, `-DYAML_DIR=`,
`-DPETSC_DIR=` and `-DHYPRE_DIR=`.  The spatial dimension can be specified with
`-DSPATIAL_DIM=` which must be either `2` or `3`.


### 3. Compilation

After configuring the build directory, the sources can be built with
```bash
cmake --build ./build [-j<ncores>]
```
Where `ncores` can be reduced in case the build process consumes too much memory
(by default all available cores are used).

---


## Running


### Serial Execution

To run OpenAccel in serial mode:

```bash
<path-to-executable> -i <path-to-input-file>
```

**Example:**
```bash
./build/openaccel-3D.exe -i case/input.yaml
```


### Parallel Execution with MPI

Running OpenAccel in parallel requires a domain decomposition.  This can either
be done manually using the `decomp` utility that ships with Trilinos or
automatically by specifying the `automatic_decomposition_type: <method>` under
the `mesh` field in the YAML input file (set in the examples to `rcb`).  The
decomposition algorithm `<method>` may be one of the following (defined by
Zoltan2):

* `hsfc`: Hilbert space-filling curve
* `rcb`: Recursive coordinate bisection (default)
* `rcb_ignore_z`: Recursive coordinate bisection, ignoring the z-coordinate
  dimension
* `rib`: Recursive intertial bisection

> [!IMPORTANT]
> The automatic decomposition is performed at runtime by STK/Ioss and requires
> Trilinos to have been linked against the **PnetCDF-enabled netCDF**.  If
> Trilinos was built manually against the serial netCDF, the automatic
> decomposition is not available and the mesh must be decomposed beforehand
> with `decomp` (see [Manual Mesh
> Decomposition](#manual-mesh-decomposition)).  See the
> [Dependencies](#required) section for how to configure this.

Execute OpenAccel in parallel via `mpirun`:

```bash
mpirun -np <num-procs> <path-to-executable> -i <path-to-input-file>
```

**Example:**
```bash
mpirun -np 8 ./build/openaccel-3D.exe -i case/input.yaml
```


#### Manual Mesh Decomposition

Should manual mesh decomposition be required, the `decomp` utility provided by
Trilinos can be used:

```bash
decomp --processors <num-procs> <path-to-mesh> --rcb --64 -V
```

**Parameters:**
- `--rcb` - Use Recursive Coordinate Bisection algorithm
- `--64` - Use 64-bit integers
- `-V` - Verbose output

**Example:**
```bash
decomp --processors 8 case/mesh.exo --rcb --64 -V
```

> **⚠️ Note on 64-bit meshes:**
> When decomposing 64-bit meshes, `nem_slice` may skip calling **ZOLTAN**, which
> is essential for producing contiguous grain-type partitions (as opposed to
> scattered partitions). It is recommended to convert the mesh to 32-bit format
> before decomposition using:
> ```bash
> ncdump <input-mesh> \
>   | sed 's/int64_status.*/int64_status = 0;/' \
>   | sed 's/\bint64\b/int/' \
>   | ncgen -5 -o <output-mesh>
> ```

---


## Examples

The repository includes a suite of test cases designed to demonstrate the
various numerical features and physical models available in the package. Below
is a categorized list of the available examples:


### Elementary Transport Phenomena

* **flange:** A pure thermal diffusion case in a solid flange part.


### Fluid Dynamics & Turbulence

* **cavity:** The classic lid-driven cavity problem; a closed domain with no
  inflow or outflow boundaries.
* **airfoil:** Incompressible flow over an airfoil featuring **non-conformal
  interfaces**, where the airfoil is discretized in a separate mesh zone.
* **elbow:** Incompressible laminar flow within a wedge-type mesh.
* **circularArc:** Compressible inviscid flow over a bump (circular arc) at
  subsonic speed.
* **pitzDaily:** Turbulent flow over a backward-facing step, used to benchmark
  turbulence models.
* **t106a:** Turbulent flow through a turbine cascade using periodic boundary
  conditions.
* **forwardStep**: 2D inviscid supersonic flow past a forward-facing step.
  > [!IMPORTANT]
  > This case requires the solver to be compiled for 2D by setting
  > `-DSPATIAL_DIM=2` during the CMake configuration.


### Multiphase & Free Surface Flow

* **damBreak:** A laminar multiphase free-surface flow simulation using the
  **Volume of Fluid (VOF)** method.
* **staticDroplet:** A stationary water droplet in an air domain. This VOF case
  isolates and tests the **Continuum Surface Force (CSF)** model for surface
  tension without gravity.


### Solid Mechanics

* **plateHole:** The Kirsch Problem benchmark uses a quarter-symmetry model to
  validate stress concentration at a circular hole against analytical solutions.


### Fluid-Structure Interaction (FSI) & Moving Meshes

* **flexibleDamBreak:** A water column collapse where the fluid impacts an
  elastic obstacle. This demonstrates the coupling of **VOF** with **ALE**
  formulations.
* **oscillatingBox:** A simple box in a closed cavity oscillating vertically,
  showcasing the **dynamic mesh** capabilities.
* **perpendicularFlap:** Incompressible flow against a flexible bar,
  highlighting the **ALE** approach for FSI.


### Heat Transfer & Buoyancy

* **BénardCells:** Buoyancy-driven flow (natural convection) utilizing the
  **Boussinesq approximation**.
* **slab:** Incompressible flow over a heated slab, demonstrating the **CHT**
  methodology.

---

## Acknowledgments

The development of OpenAccel is supported by the **Swiss National Science
Foundation** under the project *"Immersed Methods for
Fluid-Structure-Contact-Interaction Simulations and Complex Geometries"* (grant
nr. 215627), and by the **Platform for Advanced Scientific Computing (PASC)**
Program under the project *"XSES-FSI: towards eXtreme Scale Semi-Structured
discretizations for Fluid-Structure Interaction."*

---

## License

OpenAccel is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE) for
the full text.

## Contact

* **Project Coordinator:** Luca Mangani ([luca.mangani@hslu.ch](mailto:luca.mangani@hslu.ch))
* **Project Maintainer:** Lucian Hanimann ([lucian.hanimann@hslu.ch](mailto:lucian.hanimann@hslu.ch))
