#!/usr/bin/env python3
# File       : exodus_to_cgns.py
# Description: Convert a multi-block HEX8 Exodus mesh (.e / .exo) into a
#              CFX-/ParaView-importable unstructured CGNS file.
#
#   * one CGNS Unstructured Zone_t per Exodus element block, with local node
#     renumbering (blocks may have disjoint nodes, e.g. DG/non-conformal),
#   * one HEXA_8 volume section per zone,
#   * one QUAD_4 surface section per sideset, named after the sideset, attached
#     to the zone of the block it bounds. CFX-Pre reads each surface section as
#     a 2D region for boundary-condition assignment.
#
# The file is written by calling the system CGNS Mid-Level Library (libcgns)
# directly through ctypes, so the output is SIDS-compliant (passes cgnscheck).
#
#   usage:  exodus_to_cgns.py mesh.e [-o mesh.cgns]
#
# Note: mesh + boundary regions only (no solution fields are transferred).

import sys


# ---------------------------------------------------------------------------
# Dependency preflight -- check everything is available before doing any work
# ---------------------------------------------------------------------------
def check_dependencies():
    import importlib.util
    import ctypes

    missing = []

    for mod in ("numpy", "netCDF4"):
        if importlib.util.find_spec(mod) is None:
            missing.append(f"Python module '{mod}'")

    try:
        ctypes.CDLL("libcgns.so")
    except OSError:
        missing.append("system library 'libcgns.so' (CGNS Mid-Level Library)")

    if missing:
        sys.stderr.write("ERROR: missing required dependencies:\n")
        for m in missing:
            sys.stderr.write(f"  - {m}\n")
        sys.stderr.write(
            "\nInstall, e.g.:\n"
            "  pip install numpy netCDF4\n"
            "  apt-get install libcgns-dev      # Debian/Ubuntu (provides "
            "libcgns.so)\n"
            "  conda install -c conda-forge cgns\n"
        )
        sys.exit(1)


check_dependencies()

import argparse
import ctypes as C
import numpy as np
import netCDF4


# --- CGNS enums / constants (from cgnslib.h) -------------------------------
CG_MODE_WRITE = 1
Unstructured  = 3
RealDouble    = 4
QUAD_4        = 7
HEXA_8        = 17

# HEX8 face -> local node map (Exodus and CGNS share HEX8 vertex ordering)
FMAP = {1: [0, 1, 5, 4], 2: [1, 2, 6, 5], 3: [2, 3, 7, 6],
        4: [0, 4, 7, 3], 5: [0, 3, 2, 1], 6: [4, 5, 6, 7]}


def load_cgns():
    lib = C.CDLL("libcgns.so")
    ip = C.POINTER(C.c_int)
    sig = {
        "cg_open": [C.c_char_p, C.c_int, ip],
        "cg_base_write": [C.c_int, C.c_char_p, C.c_int, C.c_int, ip],
        "cg_zone_write": [C.c_int, C.c_int, C.c_char_p, C.c_void_p, C.c_int, ip],
        "cg_coord_write": [C.c_int, C.c_int, C.c_int, C.c_int, C.c_char_p,
                           C.c_void_p, ip],
        "cg_section_write": [C.c_int, C.c_int, C.c_int, C.c_char_p, C.c_int,
                             C.c_int, C.c_int, C.c_int, C.c_void_p, ip],
        "cg_close": [C.c_int],
    }
    for name, args in sig.items():
        f = getattr(lib, name)
        f.argtypes = args
        f.restype = C.c_int
    lib.cg_get_error.restype = C.c_char_p
    return lib


def decode_names(var):
    """Decode an Exodus name table (char matrix) to a list of strings."""
    out = []
    for row in var[:]:
        s = row.tobytes().decode("ascii", "ignore") if hasattr(row, "tobytes") \
            else "".join(np.asarray(row).astype(str))
        out.append(s.split("\x00", 1)[0].strip())
    return out


def cgns_name(s):
    """Sanitize to a valid 32-char CGNS node name."""
    return s[:32] if s else "unnamed"


def main():
    ap = argparse.ArgumentParser(
        description="Convert a multi-block HEX8 Exodus mesh to CGNS "
                    "(zones per block, sidesets as named surface regions).")
    ap.add_argument("input", help="input Exodus file (.e / .exo)")
    ap.add_argument("-o", "--output", default=None,
                    help="output CGNS file (default: <input>.cgns)")
    args = ap.parse_args()

    out = args.output or (args.input.rsplit(".", 1)[0] + ".cgns")

    lib = load_cgns()

    def chk(ierr, what):
        if ierr != 0:
            raise RuntimeError(f"{what}: {lib.cg_get_error().decode()}")

    d = netCDF4.Dataset(args.input)

    X = np.asarray(d.variables["coordx"][:], float)
    Y = np.asarray(d.variables["coordy"][:], float)
    Z = np.asarray(d.variables["coordz"][:], float)

    # element blocks
    blocks = []
    i = 1
    while f"connect{i}" in d.variables:
        conn = np.asarray(d.variables[f"connect{i}"][:], np.int64) - 1
        if conn.shape[1] != 8:
            sys.exit(f"ERROR: block {i} is not HEX8 (got {conn.shape[1]} "
                     "nodes/elem); only HEX8 is supported.")
        blocks.append(conn)
        i += 1
    if not blocks:
        sys.exit("ERROR: no element blocks (connect*) found.")

    blk_names = decode_names(d.variables["eb_names"]) \
        if "eb_names" in d.variables else [f"block{j+1}"
                                           for j in range(len(blocks))]
    blk_off = np.cumsum([0] + [b.shape[0] for b in blocks])   # global offsets
    allconn = np.vstack(blocks)

    # sidesets: faces (global node ids) + the block they bound
    ss_names = decode_names(d.variables["ss_names"]) \
        if "ss_names" in d.variables else []
    sidesets = {j: [] for j in range(len(blocks))}            # block -> [(name, faces)]
    for k, name in enumerate(ss_names, start=1):
        if f"elem_ss{k}" not in d.variables:
            continue
        elems = np.asarray(d.variables[f"elem_ss{k}"][:], np.int64) - 1
        sides = np.asarray(d.variables[f"side_ss{k}"][:], np.int64)
        faces = np.array([allconn[e][FMAP[int(s)]]
                          for e, s in zip(elems, sides)])
        blk = int(np.searchsorted(blk_off, elems[0], side="right") - 1)
        if not np.all((elems >= blk_off[blk]) & (elems < blk_off[blk + 1])):
            sys.stderr.write(f"  warning: sideset '{name}' spans >1 block; "
                             f"attaching to block '{blk_names[blk]}'\n")
        sidesets[blk].append((name, faces))

    # --- write CGNS --------------------------------------------------------
    fn = C.c_int()
    chk(lib.cg_open(out.encode(), CG_MODE_WRITE, C.byref(fn)), "cg_open")
    Bidx = C.c_int()
    chk(lib.cg_base_write(fn, b"Base", 3, 3, C.byref(Bidx)), "cg_base_write")

    for zi, blk in enumerate(blocks):
        ncell = blk.shape[0]
        local = np.unique(blk)
        g2l = -np.ones(X.shape[0], np.int64)
        g2l[local] = np.arange(local.size)
        nnode = local.size

        size = np.array([nnode, ncell, 0], np.int32)
        Zidx = C.c_int()
        chk(lib.cg_zone_write(fn, Bidx, cgns_name(blk_names[zi]).encode(),
                              size.ctypes.data, Unstructured, C.byref(Zidx)),
            "cg_zone_write")

        for cname, comp in [(b"CoordinateX", X), (b"CoordinateY", Y),
                            (b"CoordinateZ", Z)]:
            arr = np.ascontiguousarray(comp[local], np.float64)
            Cidx = C.c_int()
            chk(lib.cg_coord_write(fn, Bidx, Zidx, RealDouble, cname,
                                   arr.ctypes.data, C.byref(Cidx)),
                "cg_coord_write")

        hexc = np.ascontiguousarray((g2l[blk] + 1).astype(np.int32).ravel())
        Sidx = C.c_int()
        chk(lib.cg_section_write(fn, Bidx, Zidx,
                                 cgns_name(blk_names[zi] + "_cells").encode(),
                                 HEXA_8, 1, ncell, 0, hexc.ctypes.data,
                                 C.byref(Sidx)), "cg_section_write(hex)")

        start = ncell + 1
        written = []
        for name, faces in sidesets[zi]:
            nf = faces.shape[0]
            quad = np.ascontiguousarray((g2l[faces] + 1).astype(np.int32).ravel())
            end = start + nf - 1
            Sidx = C.c_int()
            chk(lib.cg_section_write(fn, Bidx, Zidx, cgns_name(name).encode(),
                                     QUAD_4, start, end, 0, quad.ctypes.data,
                                     C.byref(Sidx)), f"cg_section_write({name})")
            start = end + 1
            written.append(name)
        print(f"zone '{blk_names[zi]}': {nnode} nodes, {ncell} hex, "
              f"regions: {written}")

    chk(lib.cg_close(fn), "cg_close")
    d.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
