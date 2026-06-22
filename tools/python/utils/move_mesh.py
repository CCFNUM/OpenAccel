#!/usr/bin/env python3
# File       : move_mesh.py
# Created    : Mon Jun 15 2026
# Author     : Mhamad Mahdi Alloush
# Description: Move mesh boundaries (translate / rotate) and IDW-warp the volume
# Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.
"""
Rigid-body move one or more boundaries (side sets) of an Exodus mesh by a
translation and/or a rotation, and deform the interior mesh accordingly using
Inverse-Distance-Weighting (IDW) warping -- the classic "IDWarp" strategy.

Every boundary node carries a prescribed displacement: the moving boundaries
get the displacement of their rigid transform; every *other* boundary is fixed
in space (zero displacement) by default. Each interior node then receives the
inverse-distance-weighted average of those boundary displacements,

    u_i = sum_b w_ib u_b / sum_b w_ib ,   w_ib = 1 / d_ib**power ,

summed over all boundary (control) nodes. Because the displacement field -- not
the rotation angle -- is interpolated, it decays smoothly to zero at the fixed
boundaries, the near-wall layers follow their wall almost rigidly (boundary
layer / y+ preserved), and the far field stays put. The transform is exact on
the boundaries themselves.

Works for 2D and 3D meshes. Element non-inversion is verified after warping.

A common use is to impose a geometric angle of attack on an airfoil mesh so the
free stream can be set mesh-aligned (alpha = 0 in the solver), e.g.

    python3 move_mesh.py mesh.e -o mesh_aoa.e --move wing:rotate:-2.31

Move specification (repeatable --move), target is a side-set name or integer id:

    NAME:translate:dx,dy[,dz]
    NAME:rotate:deg[,px,py[,pz][,ax,ay,az]]

For rotate, the pivot defaults to the centroid of the moving boundary and the
axis defaults to +z (the only choice in 2D). Sign: a positive angle is a
counter-clockwise (right-hand-rule about the axis) rotation; for a +x free
stream a NOSE-UP airfoil (positive angle of attack) is "rotate:-alpha" about
the spanwise +z axis.

Usage:
    python3 move_mesh.py INPUT.e [-o OUTPUT.e] --move SPEC [--move SPEC ...]
        [--power 4.0] [--free NAME ...] [--chunk 4000]
"""
import argparse
import shutil
import sys
import numpy as np
from netCDF4 import Dataset


# --------------------------------------------------------------------------- #
# Exodus helpers
# --------------------------------------------------------------------------- #
def _arr(v):
    return np.ma.getdata(v[:])


def read_sideset_table(d):
    """Return {name: ssid} and {ssid: index(1-based)} for the side sets."""
    ids = list(_arr(d.variables["ss_prop1"])) if "ss_prop1" in d.variables else []
    names = []
    if "ss_names" in d.variables:
        for row in _arr(d.variables["ss_names"]):
            names.append(b"".join(row[row != b""]).decode().strip("\x00"))
    name2id = {nm: int(i) for nm, i in zip(names, ids)}
    id2index = {int(i): k + 1 for k, i in enumerate(ids)}  # elem_ss<index>
    return name2id, id2index, names, [int(i) for i in ids]


def sideset_nodes(d, index, blocks0):
    """0-based node ids touched by faces of the side set at this 1-based index.

    Exodus side sets reference (elem, local-side); we take every node of the
    parent elements, which is exact for rigid-body boundary motion (all nodes
    of a boundary face move with the body)."""
    elems = _arr(d.variables["elem_ss%d" % index]).astype(np.int64) - 1
    nodes = []
    for conn0, lo, hi in blocks0:
        m = (elems >= lo) & (elems < hi)
        if np.any(m):
            nodes.append(conn0[elems[m] - lo].ravel())
    if not nodes:
        return np.empty(0, dtype=np.int64)
    return np.unique(np.concatenate(nodes))


def read_blocks(d):
    """List of (conn0, lo, hi): 0-based connectivity and global elem range."""
    blocks, off = [], 0
    k = 1
    while ("connect%d" % k) in d.variables:
        conn0 = _arr(d.variables["connect%d" % k]).astype(np.int64) - 1
        ne = conn0.shape[0]
        blocks.append((conn0, off, off + ne))
        off += ne
        k += 1
    return blocks


# --------------------------------------------------------------------------- #
# Geometry / transforms
# --------------------------------------------------------------------------- #
def rotation_matrix(deg, axis, dim):
    a = np.deg2rad(deg)
    c, s = np.cos(a), np.sin(a)
    if dim == 2:
        return np.array([[c, -s], [s, c]])
    ax = np.asarray(axis, float)
    ax = ax / np.linalg.norm(ax)
    x, y, z = ax
    C = 1.0 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def signed_measure(x, conn0, dim):
    """Per-element signed area (2D) or signed volume-Jacobian (3D).

    Returns None for topologies we do not check (caller skips them)."""
    npe = conn0.shape[1]
    P = x[conn0]  # (ne, npe, dim)
    if dim == 2:  # shoelace, exact for tri/quad
        xs, ys = P[:, :, 0], P[:, :, 1]
        xr, yr = np.roll(xs, -1, 1), np.roll(ys, -1, 1)
        return 0.5 * np.sum(xs * yr - xr * ys, axis=1)
    if dim == 3 and npe == 4:  # tet determinant
        a = P[:, 1] - P[:, 0]
        b = P[:, 2] - P[:, 0]
        c = P[:, 3] - P[:, 0]
        return np.einsum("ij,ij->i", np.cross(a, b), c) / 6.0
    if dim == 3 and npe == 8:  # hex centroid Jacobian sign (trilinear map)
        sx = np.array([-1, 1, 1, -1, -1, 1, 1, -1]) / 8.0
        se = np.array([-1, -1, 1, 1, -1, -1, 1, 1]) / 8.0
        sz = np.array([-1, -1, -1, -1, 1, 1, 1, 1]) / 8.0
        J0 = np.einsum("a,nad->nd", sx, P)
        J1 = np.einsum("a,nad->nd", se, P)
        J2 = np.einsum("a,nad->nd", sz, P)
        return np.einsum("ij,ij->i", np.cross(J0, J1), J2)
    return None


# --------------------------------------------------------------------------- #
# Move-spec parsing
# --------------------------------------------------------------------------- #
def parse_move(spec, dim):
    """'target:op:args' -> dict(target, op, params)."""
    parts = spec.split(":")
    if len(parts) < 3:
        raise ValueError("bad --move '%s' (need target:op:args)" % spec)
    target, op = parts[0], parts[1].lower()
    vals = [float(v) for v in parts[2].split(",")]
    if op == "translate":
        if len(vals) < dim:
            raise ValueError("translate needs %d components in '%s'" % (dim, spec))
        return dict(target=target, op="translate", vec=np.array(vals[:dim]))
    if op == "rotate":
        deg = vals[0]
        pivot = None
        axis = [0, 0, 1]
        rest = vals[1:]
        if dim == 2 and len(rest) >= 2:
            pivot = np.array(rest[:2])
        if dim == 3:
            if len(rest) >= 3:
                pivot = np.array(rest[:3])
            if len(rest) >= 6:
                axis = rest[3:6]
        return dict(target=target, op="rotate", deg=deg, pivot=pivot, axis=axis)
    raise ValueError("unknown op '%s' in '%s'" % (op, spec))


def resolve_target(target, name2id, valid_ids):
    if target in name2id:
        return name2id[target]
    try:
        ssid = int(target)
    except ValueError:
        raise ValueError("side set '%s' not found by name" % target)
    if ssid not in valid_ids:
        raise ValueError("side set id %d not in mesh" % ssid)
    return ssid


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("-o", "--out", dest="outfile", default=None)
    ap.add_argument("--move", action="append", default=[], required=True,
                    metavar="SPEC", help="target:op:args (repeatable)")
    ap.add_argument("--free", action="append", default=[], metavar="NAME",
                    help="boundary left unconstrained (slides; not a control set)")
    ap.add_argument("--power", type=float, default=4.0,
                    help="IDW exponent (higher => more local/rigid near wall)")
    ap.add_argument("--chunk", type=int, default=4000)
    args = ap.parse_args()

    out = args.outfile or args.infile.rsplit(".", 1)[0] + "_warped.e"
    shutil.copyfile(args.infile, out)
    d = Dataset(out, "r+")

    dim = d.dimensions["num_dim"].size
    coord_names = ["coordx", "coordy", "coordz"][:dim]
    X = np.column_stack([_arr(d.variables[c]).astype(np.float64)
                         for c in coord_names])
    nnode = X.shape[0]
    blocks0 = read_blocks(d)
    name2id, id2index, names, valid_ids = read_sideset_table(d)

    moves = [parse_move(s, dim) for s in args.move]
    moving_ids = {resolve_target(m["target"], name2id, valid_ids) for m in moves}
    free_ids = {resolve_target(t, name2id, valid_ids) for t in args.free}

    # prescribed boundary displacements
    U = np.zeros((nnode, dim))
    is_moving = np.zeros(nnode, dtype=bool)
    for m in moves:
        ssid = resolve_target(m["target"], name2id, valid_ids)
        nodes = sideset_nodes(d, id2index[ssid], blocks0)
        if m["op"] == "translate":
            U[nodes] = m["vec"]
        else:
            pivot = m["pivot"]
            if pivot is None:
                pivot = X[nodes].mean(axis=0)
            R = rotation_matrix(m["deg"], m["axis"], dim)
            U[nodes] = (X[nodes] - pivot) @ R.T - (X[nodes] - pivot)
        is_moving[nodes] = True

    # control set = all boundary nodes except the "free" ones; moving wins
    is_ctrl = np.zeros(nnode, dtype=bool)
    fixed_count = 0
    for ssid in valid_ids:
        if ssid in moving_ids or ssid in free_ids:
            continue
        nodes = sideset_nodes(d, id2index[ssid], blocks0)
        nodes = nodes[~is_moving[nodes]]          # shared node -> moving wins
        U[nodes] = 0.0
        is_ctrl[nodes] = True
        fixed_count += nodes.size
    is_ctrl[is_moving] = True                     # moving nodes are controls too

    ctrl = np.where(is_ctrl)[0]
    interior = np.where(~is_ctrl)[0]
    print("dim=%d nodes=%d  moving=%d fixed=%d free-excluded interior=%d"
          % (dim, nnode, int(is_moving.sum()), fixed_count, interior.size))
    for m in moves:
        print("  move %-14s %s %s" % (m["target"], m["op"],
              m.get("vec", "%.4f deg" % m.get("deg", 0))))

    # IDW interpolation of the displacement vector
    Xc = X[ctrl]
    Uc = U[ctrl]
    p = args.power
    eps = 1e-300
    for k in range(0, interior.size, args.chunk):
        idx = interior[k:k + args.chunk]
        diff = X[idx][:, None, :] - Xc[None, :, :]
        d2 = np.einsum("ncd,ncd->nc", diff, diff)
        w = 1.0 / (d2 ** (0.5 * p) + eps)
        U[idx] = (w @ Uc) / w.sum(axis=1)[:, None]

    Xnew = X + U

    # non-inversion check across all blocks
    bad = 0
    worst = 1.0
    for conn0, lo, hi in blocks0:
        m0 = signed_measure(X, conn0, dim)
        m1 = signed_measure(Xnew, conn0, dim)
        if m0 is None:
            print("  (validity check skipped for npe=%d)" % conn0.shape[1])
            continue
        bad += int(np.sum(np.sign(m0) != np.sign(m1)))
        r = np.abs(m1) / np.maximum(np.abs(m0), 1e-300)
        worst = min(worst, float(r.min()))
    print("element check: inverted=%d  min|new/old|=%.4f" % (bad, worst))
    if bad:
        print("ERROR: %d inverted elements - lower the move or raise --power"
              % bad, file=sys.stderr)
        d.close()
        sys.exit(1)

    for j, c in enumerate(coord_names):
        d.variables[c][:] = Xnew[:, j]
    d.close()
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
