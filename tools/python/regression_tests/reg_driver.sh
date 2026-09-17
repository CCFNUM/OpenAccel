#!/usr/bin/env bash
# File       : reg_driver.sh
# Created    : Thu Sep 17 2026 07:46:14 (+0200)
# Author     : Fabian Wermelinger
# Description: Regression driver wrapper
# Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

usage() {
    cat <<EOF
USAGE: $0 [-h|--help] [-t|-omp-num-threads <n>] -2D|-3D <exec> [-2D|-3D <exec>] [quick_test.py args]

SYNOPSIS: Executable wrapper around quick_test.py python script.  Either the
          -2D or -3D option (or both) must be specified.

OPTIONS:
    -h, --help
        Print this help message and exit.

    -t, --omp-num-threads
        Specify the number of OMP_NUM_THREADS to be used (default 1)

    -2D <exec>
        Path to 2D test executable.

    -3D <exec>
        Path to 3D test executable.

    [quick_test.py args]
        Remaining arguments not recognized by this script are passed to the
        quick_test.py python script.
EOF
    exit 1
}

exec_2d=""
exec_3d=""
n_threads=1
pass_through=()
while [ $# -gt 0 ]; do
    case $1 in
        -h|--help) usage ;;
        -t|--omp-num-threads) n_threads=$2; shift 2 ;;
        -2D) exec_2d=$2; shift 2 ;;
        -3D) exec_3d=$2; shift 2 ;;
        *) pass_through+=("$1"); shift ;;
    esac
done

if [ ! -x "${exec_2d}" ] && [ ! -x "${exec_3d}" ]; then
    echo "Either the -2D or -3D must be specified"
    exit 1
fi

args=()
if [ -x "${exec_2d}" ]; then
    args+=(--exe-2d "$(readlink -f "${exec_2d}")")
fi
if [ -x "${exec_3d}" ]; then
    args+=(--exe-3d "$(readlink -f "${exec_3d}")")
fi

SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
PROJECT_ROOT="${SCRIPT_DIR%/tools/python/regression_tests}"

export OMP_NUM_THREADS=${n_threads}
(cd "${PROJECT_ROOT}" && python3 "${SCRIPT_DIR}/quick_test.py" "${args[@]}" "${pass_through[@]}")
