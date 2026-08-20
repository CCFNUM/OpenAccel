#!/usr/bin/env bash
# Run an OpenAccel solver with the correct MPI on this machine.
#
# ~/.bashrc sources OpenFOAM, which puts the system OpenMPI 4 ahead of the
# spack OpenMPI 5 that OpenAccel is built against. Mixing them fails at load
# time with:
#     libmpi_mpifh.so.40: undefined symbol: ompi_instance_count
# This wrapper prepends the spack MPI for the child process only, so OpenFOAM
# in the same shell is unaffected.
#
# Usage:
#   ./run-openaccel.sh -i input.i                 # 2D solver (default)
#   EXE=build-3d/openaccel-3D.exe ./run-openaccel.sh -i input.i
#   NP=4 ./run-openaccel.sh -i input.i            # under mpirun
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPACK_MPI_LIB="$(dirname "$(ldd "$REPO/${EXE:-build-2d/openaccel-2D.exe}" \
    | awk '/libmpi_mpifh\.so/ {print $3}')")"

if [[ -z "$SPACK_MPI_LIB" || ! -d "$SPACK_MPI_LIB" ]]; then
    echo "could not locate the spack MPI lib dir" >&2
    exit 1
fi

export LD_LIBRARY_PATH="$SPACK_MPI_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

EXE_PATH="$REPO/${EXE:-build-2d/openaccel-2D.exe}"
if [[ -n "${NP:-}" ]]; then
    exec mpirun --oversubscribe -np "$NP" "$EXE_PATH" "$@"
else
    exec "$EXE_PATH" "$@"
fi
