# OpenAccel
**OpenAccel** is a CPU-parallel, vertex-based finite volume (**CVFEM**) solver built on the **Trilinos-STK** mesh infrastructure. It employs a **pressure-based, segregated approach** to solve the governing equations, making it well-suited for incompressible and low-Mach compressible flows. The solver addresses a broad range of physics, including fluid flow, heat transfer, turbulence, solid mechanics, and multiphase free-surface flows.

> **Note:** This project is under active development. APIs, features, and documentation may change without notice.

---

## 🔬 Capabilities
* **Incompressible and compressible single-phase flows**
* **Incompressible multiphase flows** using the Volume of Fluid (VOF) method
* **Fluid-Structure Interaction (FSI)** via a partitioned Arbitrary Lagrangian-Eulerian (ALE) approach
* **Heat Transfer** including Conjugate Heat Transfer (CHT)
* **Turbulence Modeling** with support for RANS (e.g., k-ω SST and transition SST)

## 🚀 Features
* **Modern Foundation:** Built using the **C++20** standard.
* **High Performance:** Optimized for massive parallelism via Trilinos-STK.

## 📖 Documentation

A **Theory Guide** covering the mathematical foundations of OpenAccel is available in [`tools/docs/theory/`](./tools/docs/theory/).

The guide is automatically compiled to PDF on every commit via GitHub Actions. To get the latest PDF:

1. Go to the **Actions** tab of this repository.
2. Open the most recent **Build Theory Guide** workflow run.
3. Download the `theory-guide-pdf` artifact.

Topics covered include governing equations, turbulence models (k-ε, k-ω SST, Transition SST), heat transfer, mesh deformation (ALE/displacement diffusion), free surface flow (VoF + FCT/cMULES), and discretisation schemes.

---

## 📚 Dependencies

### Required
* **Trilinos:** Must be built with **STK**. **Tpetra** and **Belos** are optionally enabled.
* **YAML-cpp:** Required for input parsing and configuration.
* **MPI:** Essential for parallel execution and distributed memory
  communication.  [MPICH](https://www.mpich.org/) and
  [OpenMPI](https://www.open-mpi.org/) implementations are tested.

### Optional (Backend Solvers)
* **PETSc:** Optional high-performance solver backend. (Some examples are
  configured with PETSc solvers. See [`examples/cavity/input.i`](./examples/cavity/input.i) and [`examples/elbow/input.i`](./examples/elbow/input.i).)
* **HYPRE:** Optional multigrid preconditioner support. (See [`examples/airfoil/input.i`](./examples/airfoil/input.i) and [`examples/pitzDaily/input.i`](./examples/pitzDaily/input.i).)
* **gnuplot:** Optional live residual plotting.

### Submodules
The following libraries are included as Git submodules:
* **[Eigen](https://gitlab.com/libeigen/eigen):** Template library for linear algebra.
* **[ExprTk](https://github.com/ArashPartow/exprtk):** Mathematical expression parsing and evaluation.
* **[nanoflann](https://github.com/jlblancoc/nanoflann):** KD-tree library for nearest-neighbour searches.
* **[gplotpp](https://github.com/ziotom78/gplotpp):** C++ interface for gnuplot.
* **[liblinsolve](https://gitlab.com/hslu_ccfnum/liblinsolve.git):** Wrapper library for HYPRE, PETSc, and Trilinos linear solver utilities.

---

## 📂 Installation

### 1. Cloning the Repository
Retrieve the source code and initialize submodules:

```bash
git clone https://github.com/CCFNUM/OpenAccel.git
cd OpenAccel
git submodule update --init --recursive
```

### 2. Configuration with CMake

Create a build directory and configure the project:

```bash
mkdir build
cd build

cmake -DCMAKE_CXX_EXTENSIONS=Off \
    -DTrilinos_DIR=<trilinos-install-cmake-directory> \
    -DYAML_DIR=<yaml-cpp-install-directory> \
    -DPETSC_DIR=<petsc-install-directory> \
    -DHYPRE_DIR=<hypre-install-directory> \
    -DGNUPLOT_DIR=<gnuplot-install-directory> \
    ..
```

### 3. Compilation

Once configuration is complete, compile the code using `make`:

```bash
make -j$(nproc)
```

The `-j$(nproc)` flag enables parallel compilation using all available CPU cores. Adjust the number if needed (e.g., `make -j4` for 4 cores).

---

## 🖥️ Running

### Serial Execution

To run OpenAccel in serial mode:

```bash
<path-to-executable> -i <path-to-input-file>
```

**Example:**
```bash
./build/accel-3D.exe -i case/input.yaml
```

### Parallel Execution with MPI

Running OpenAccel in parallel requires two steps: mesh decomposition and parallel execution.

#### Step 1: Mesh Decomposition

Use the `decomp` utility provided by Trilinos to decompose the mesh:

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
> When decomposing 64-bit meshes, `nem_slice` may skip calling **ZOLTAN**, which is essential for producing contiguous grain-type partitions (as opposed to scattered partitions). It is recommended to convert the mesh to 32-bit format before decomposition using:
> ```bash
> ncdump <input-mesh> \
>   | sed 's/int64_status.*/int64_status = 0;/' \
>   | sed 's/\bint64\b/int/' \
>   | ncgen -5 -o <output-mesh>
> ```

#### Step 2: Run with MPI

Execute OpenAccel using `mpirun`:

```bash
mpirun -np <num-procs> <path-to-executable> -i <path-to-input-file>
```

**Example:**
```bash
mpirun -np 8 ./build/accel-3D.exe -i case/input.yaml
```

---

## Examples

The repository includes a suite of test cases designed to demonstrate the various numerical features and physical models available in the package. Below is a categorized list of the available examples:

### Basic Transport Phenomena
* **flange:** A pure thermal diffusion case in a solid flange part.

### Fluid Dynamics & Turbulence
* **cavity:** The classic lid-driven cavity problem; a closed domain with no inflow or outflow boundaries.
* **airfoil:** Incompressible flow over an airfoil featuring **non-conformal interfaces**, where the airfoil is discretized in a separate mesh zone.
* **elbow:** Incompressible laminar flow within a wedge-type mesh.
* **circularArc:** Compressible inviscid flow over a bump (circular arc) at subsonic speed.
* **pitzDaily:** Incompressible turbulent internal flow over a backward-facing step. This case utilizes the **k-ω SST** turbulence model.

### Multiphase & Free Surface Flow
* **damBreak:** A laminar multiphase free-surface flow simulation using the **Volume of Fluid (VOF)** method.
* **staticDroplet:** A stationary water droplet in an air domain. This VOF case isolates and tests the **Continuum Surface Force (CSF)** model for surface tension without gravity.

### Fluid-Structure Interaction (FSI) & Moving Meshes
* **flexibleDamBreak:** A water column collapse where the fluid impacts an elastic obstacle. This demonstrates the coupling of **VOF** with **Arbitrary Lagrangian-Eulerian (ALE)** formulations.
* **oscillatingBox:** A simple box in a closed cavity oscillating vertically, showcasing the **dynamic mesh** capabilities.
* **perpendicularFlap:** Incompressible flow against a flexible bar, highlighting the **ALE** approach for FSI.

### Heat Transfer & Buoyancy
* **BénardCells:** Buoyancy-driven flow (natural convection) utilizing the **Boussinesq approximation**.
* **slab:** Incompressible flow over a heated slab, demonstrating the **Conjugate Heat Transfer (CHT)** methodology.

---

## Acknowledgments

The development of OpenAccel is supported by the **Swiss National Science Foundation** under the project *"Immersed Methods for Fluid-Structure-Contact-Interaction Simulations and Complex Geometries"* (grant nr. 215627), and by the **Platform for Advanced Scientific Computing (PASC)** Program under the project *"XSES-FSI: towards eXtreme Scale Semi-Structured discretizations for Fluid-Structure Interaction."*

---

## License

OpenAccel is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE) for the full text.

## Contact

* **Project Coordinator:** Luca Mangani ([luca.mangani@hslu.ch](mailto:luca.mangani@hslu.ch))
* **Project Maintainer:** Lucian Hanimann ([lucian.hanimann@hslu.ch](mailto:lucian.hanimann@hslu.ch))
