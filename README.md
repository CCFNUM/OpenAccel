# OpenAccel
**OpenAccel** is a parallel vertex-based finite volume fluid flow solver (**CVFEM**) powered by the **STK library**.

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
* **HYPRE:** Optional multigrid preconditioner support.
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

## License

OpenAccel is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE) for the full text.

## Contact

* **Project Coordinator:** Luca Mangani ([luca.mangani@hslu.ch](mailto:luca.mangani@hslu.ch))
* **Project Maintainer:** Lucian Hanimann ([lucian.hanimann@hslu.ch](mailto:lucian.hanimann@hslu.ch))
