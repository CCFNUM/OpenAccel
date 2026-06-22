#!/usr/bin/env python3
# File       : split_sideset_by_normal.py
# Created    : Mon Jun 15 2026
# Author     : Mhamad Mahdi Alloush
# Description: Split Exodus side set(s) into two by outward-normal . reference
# Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.
"""
Split one or more Exodus side sets into two new side sets according to the sign
of (outward face normal) . (reference vector). Faces whose outward normal points
*against* the reference vector by more than a margin go to the first ("inflow")
set; all other faces go to the second ("outflow") set:

    n . v < -threshold   ->  inflow set
    otherwise            ->  outflow set

Typical use for external aerodynamics: split a closed/circular far field into a
windward INLET (where the free stream enters) and a static-pressure OUTLET
(top/bottom/leeward). With the inlet restricted to clearly-windward faces, the
imposed-direction inlet stays well-conditioned (its dir.n is bounded away from
zero), while the outlet -- which imposes no direction -- handles the grazing and
leeward faces without the tangential degeneracy or slip-wall clipping that a
single direction-imposed far field suffers at angle of attack.

Side sets not listed for splitting are copied through unchanged. Everything else
in the mesh (coordinates, blocks, names, ...) is preserved; only the side-set
section is rebuilt. 2D meshes (QUAD4 / TRI3 edges) are supported.

Usage:
    python3 split_sideset_by_normal.py INPUT.e -o OUTPUT.e \
        --sidesets farfield_in,farfield_out \
        --vector 0.99918744,0.04030615 \
        --names inlet,outlet --threshold 0.15
"""
import argparse
import sys
import numpy as np
from netCDF4 import Dataset

# local node ordinals (0-based) of each side, per element type
SIDE_NODES = {
    "QUAD4": {1: (0, 1), 2: (1, 2), 3: (2, 3), 4: (3, 0)},
    "QUAD": {1: (0, 1), 2: (1, 2), 3: (2, 3), 4: (3, 0)},
    "TRI3": {1: (0, 1), 2: (1, 2), 3: (2, 0)},
    "TRI": {1: (0, 1), 2: (1, 2), 3: (2, 0)},
}


def _arr(v):
    return np.ma.getdata(v[:])


def _name_list(d, var):
    return [b"".join(r[r != b""]).decode().strip("\x00") for r in _arr(d.variables[var])]


def read_blocks(d):
    """(conn0, lo, hi, elem_type) per block; conn0 is 0-based."""
    blocks, off, k = [], 0, 1
    while ("connect%d" % k) in d.variables:
        c = d.variables["connect%d" % k]
        conn0 = _arr(c).astype(np.int64) - 1
        et = getattr(c, "elem_type", "QUAD4").upper()
        blocks.append((conn0, off, off + conn0.shape[0], et))
        off += conn0.shape[0]
        k += 1
    return blocks


def outward_normal(coords, elem_nodes0, side):
    """Unit outward normal of a 2D element edge `side` (1-based)."""
    et_nodes = elem_nodes0
    npe = len(et_nodes)
    et = "QUAD4" if npe == 4 else "TRI3"
    ia, ib = SIDE_NODES[et][side]
    a = coords[et_nodes[ia]]
    b = coords[et_nodes[ib]]
    e = b - a
    n = np.array([e[1], -e[0]])              # one of the two edge normals
    centroid = coords[et_nodes].mean(axis=0)
    mid = 0.5 * (a + b)
    if np.dot(n, mid - centroid) < 0:        # flip to point outward
        n = -n
    ln = np.linalg.norm(n)
    return n / ln if ln > 0 else n


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--sidesets", required=True,
                    help="comma list of side-set names or ids to split")
    ap.add_argument("--vector", required=True,
                    help="reference vector 'x,y' (e.g. freestream unit dir)")
    ap.add_argument("--names", default="inflow,outflow",
                    help="names for the two new side sets (inflow,outflow)")
    ap.add_argument("--threshold", type=float, default=0.0,
                    help="margin: n.v < -threshold => inflow set")
    args = ap.parse_args()

    src = Dataset(args.infile)
    dim = src.dimensions["num_dim"].size
    if dim != 2:
        sys.exit("only 2D meshes supported (got num_dim=%d)" % dim)

    coords = np.column_stack([_arr(src.variables["coordx"]),
                              _arr(src.variables["coordy"])]).astype(float)
    blocks = read_blocks(src)
    v = np.array([float(t) for t in args.vector.split(",")])
    v = v / np.linalg.norm(v)
    new_names = args.names.split(",")

    ss_ids = list(_arr(src.variables["ss_prop1"]))
    ss_names = _name_list(src, "ss_names")
    name2idx = {nm: i for i, nm in enumerate(ss_names)}      # 0-based ss position
    id2idx = {int(i): i0 for i0, i in enumerate(ss_ids)}

    split_idx = set()
    for tok in args.sidesets.split(","):
        if tok in name2idx:
            split_idx.add(name2idx[tok])
        else:
            split_idx.add(id2idx[int(tok)])

    # gather (elem, side) of every face to split, classify by normal . v
    grpA, grpB = [], []   # (elem_1based, side)
    for i0 in split_idx:
        i1 = i0 + 1
        elems = _arr(src.variables["elem_ss%d" % i1]).astype(np.int64)
        sides = _arr(src.variables["side_ss%d" % i1]).astype(np.int64)
        for e1, s in zip(elems, sides):
            e0 = e1 - 1
            for conn0, lo, hi, et in blocks:
                if lo <= e0 < hi:
                    nodes0 = conn0[e0 - lo]
                    break
            n = outward_normal(coords, nodes0, int(s))
            if np.dot(n, v) < -args.threshold:
                grpA.append((int(e1), int(s)))
            else:
                grpB.append((int(e1), int(s)))
    print("split %d faces -> %s=%d (inflow)  %s=%d (outflow)  (vector=%s thr=%.3f)"
          % (len(grpA) + len(grpB), new_names[0], len(grpA),
             new_names[1], len(grpB), v.tolist(), args.threshold))

    # assemble the new side-set list: kept (unsplit) first, then the two new
    kept = [i0 for i0 in range(len(ss_names)) if i0 not in split_idx]
    next_id = (max(int(i) for i in ss_ids) + 1) if ss_ids else 1

    out_sets = []   # (name, id, elem_array, side_array)
    for i0 in kept:
        i1 = i0 + 1
        out_sets.append((ss_names[i0], int(ss_ids[i0]),
                         _arr(src.variables["elem_ss%d" % i1]).astype(np.int32),
                         _arr(src.variables["side_ss%d" % i1]).astype(np.int32)))
    for grp, nm in ((grpA, new_names[0]), (grpB, new_names[1])):
        if not grp:
            print("  WARNING: '%s' set is empty - skipping" % nm)
            continue
        ea = np.array([g[0] for g in grp], dtype=np.int32)
        sa = np.array([g[1] for g in grp], dtype=np.int32)
        out_sets.append((nm, next_id, ea, sa))
        next_id += 1

    # ----- write new file: copy everything, rebuild side-set section -----
    skip_dims = {"num_side_sets"}
    skip_dims |= {"num_side_ss%d" % (i + 1) for i in range(len(ss_names))}
    skip_dims |= {"num_df_ss%d" % (i + 1) for i in range(len(ss_names))}
    skip_vars = {"ss_status", "ss_prop1", "ss_names"}
    for i in range(len(ss_names)):
        skip_vars |= {"elem_ss%d" % (i + 1), "side_ss%d" % (i + 1),
                      "dist_fact_ss%d" % (i + 1)}

    dst = Dataset(args.out, "w", format=src.data_model)
    for a in src.ncattrs():
        dst.setncattr(a, getattr(src, a))
    for name, dimobj in src.dimensions.items():
        if name in skip_dims:
            continue
        dst.createDimension(name, None if dimobj.isunlimited() else dimobj.size)
    dst.createDimension("num_side_sets", len(out_sets))
    for k, (nm, sid, ea, sa) in enumerate(out_sets, start=1):
        dst.createDimension("num_side_ss%d" % k, ea.size)

    for name, var in src.variables.items():
        if name in skip_vars:
            continue
        nv = dst.createVariable(name, var.dtype, var.dimensions)
        for a in var.ncattrs():
            nv.setncattr(a, getattr(var, a))
        nv[:] = var[:]

    ln = src.dimensions["len_name"].size
    ss_status = dst.createVariable("ss_status", "i4", ("num_side_sets",))
    ss_prop1 = dst.createVariable("ss_prop1", "i4", ("num_side_sets",))
    ss_prop1.setncattr("name", "ID")
    ss_names_v = dst.createVariable("ss_names", "S1",
                                    ("num_side_sets", "len_name"))
    ss_status[:] = np.ones(len(out_sets), dtype=np.int32)
    ss_prop1[:] = np.array([s[1] for s in out_sets], dtype=np.int32)
    names_char = np.zeros((len(out_sets), ln), dtype="S1")
    for k, (nm, sid, ea, sa) in enumerate(out_sets):
        for j, ch in enumerate(nm.encode()[:ln - 1]):
            names_char[k, j] = bytes([ch])
        ev = dst.createVariable("elem_ss%d" % (k + 1), "i4",
                                ("num_side_ss%d" % (k + 1),))
        sv = dst.createVariable("side_ss%d" % (k + 1), "i4",
                                ("num_side_ss%d" % (k + 1),))
        ev[:] = ea
        sv[:] = sa
    ss_names_v[:] = names_char

    src.close()
    dst.close()
    print("wrote %s with %d side sets: %s"
          % (args.out, len(out_sets), [s[0] for s in out_sets]))


if __name__ == "__main__":
    main()
