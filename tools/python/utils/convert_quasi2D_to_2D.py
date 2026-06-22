#!/usr/bin/env python3
# File       : convert_quasi2D_to_2D.py
# Created    : Thu May 28 2026
# Author     : Mhamad Mahdi Alloush
# Description: Convert quasi-2D exodus to literal 2D
# Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.
"""
Convert a quasi-2D Exodus mesh (a 3D mesh that is exactly one element thick in
one direction) into a *literal* 2D Exodus mesh (num_dim = 2).

The "depth" direction is detected automatically as the coordinate axis that has
only two distinct values (i.e. a single element layer). It may be x, y or z.
The two remaining axes become the planar coordinates. Volume elements collapse
to their planar face:

    HEX8  / HEX    ->  QUAD4
    WEDGE6 / WEDGE ->  TRI3

Side sets are remapped: lateral faces become element edges; the two
depth-normal faces (front/back) collapse to the element itself and are dropped.
A side set that becomes empty (e.g. "front", "back", "frontAndBack") is removed
and reported. Node sets, block IDs, side-set IDs and all names are preserved.

Coincident nodes are merged by snapping each kept coordinate onto an absolute
per-axis grid of cell = (axis extent) * 10**(-sig-figs) (default sig-figs=10).
A per-axis grid (rather than one tolerance from the global bounding box) keeps
distinct fine-mesh nodes near a body (e.g. an airfoil boundary layer) apart,
while still floating-point "dust" on a symmetry centreline (an exact 0.0 paired
with its 1e-18 translate) onto the same grid point so the pair merges.

The conversion is fully vectorised (numpy), so meshes with millions of elements
convert in seconds.

This is a mesh-only conversion. After converting you must also adapt the
companion input.i: drop any boundary whose location points at a dropped
front/back side set, and use 2-component vectors.

Usage:
    python3 convert_quasi2D_to_2D.py INPUT.e [OUTPUT.e]
    python3 convert_quasi2D_to_2D.py INPUT.e -o OUTPUT.e \
            [--depth-axis x|y|z] [--sig-figs N] [--dry-run]

If OUTPUT is omitted it defaults to "<input-stem>_2D.e" next to the input.
"""

import argparse
import os

import numpy as np
import netCDF4


# ---------------------------------------------------------------------------
# Topology: exodus local side -> ordered list of local node indices (1-based).
# Reference: Exodus II side-set node ordering.
# ---------------------------------------------------------------------------
SIDE_NODES = {
    "HEX": {
        1: [1, 2, 6, 5], 2: [2, 3, 7, 6], 3: [3, 4, 8, 7],
        4: [4, 1, 5, 8], 5: [1, 4, 3, 2], 6: [5, 6, 7, 8],
    },
    "WEDGE": {
        1: [1, 2, 5, 4], 2: [2, 3, 6, 5], 3: [3, 1, 4, 6],
        4: [1, 3, 2], 5: [4, 5, 6],
    },
}
# 2D element edges as 0-based (i, j) endpoint pairs, in local-side order.
EDGE_IDX = {
    "QUAD4": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "TRI3":  [(0, 1), (1, 2), (2, 0)],
}
# All edges of the 3D source element as 0-based local-node pairs. Used to pair
# front/back nodes topologically: an edge whose two endpoints sit on opposite
# depth planes is a "depth edge", and its endpoints are the same 2D node.
ELEM_EDGES = {
    "HEX":   [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7)],
    "WEDGE": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3),
              (0, 3), (1, 4), (2, 5)],
}
# nodes-per-element of the source -> (canonical source key, 2D type, planar n)
SOURCE_INFO = {
    8: ("HEX",   "QUAD4", 4),
    6: ("WEDGE", "TRI3",  3),
}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def decode_name_rows(var):
    """Decode an Exodus char-array variable into a list of python strings."""
    out = []
    for row in np.ma.filled(var[:], b"\x00"):
        raw = b"".join(bytes(c) for c in row.tolist())
        out.append(raw.split(b"\x00", 1)[0].decode("ascii", "ignore").strip())
    return out


def write_name_rows(var, names):
    """Write a list of strings into an Exodus char-array variable."""
    nrow, ncol = var.shape
    buf = np.zeros((nrow, ncol), dtype="S1")
    for i, name in enumerate(names):
        for j, ch in enumerate(name.encode("ascii", "ignore")[: ncol - 1]):
            buf[i, j] = bytes([ch])
    var[:] = buf


def read_coords(ds):
    """Return (ndim, list-of-coordinate-arrays, list-of-axis-names)."""
    ndim = ds.dimensions["num_dim"].size
    letters = ["x", "y", "z"][:ndim]
    if "coord" in ds.variables:
        coord = np.asarray(ds.variables["coord"][:])
        arrs = [np.asarray(coord[i], dtype=float) for i in range(ndim)]
    else:
        arrs = [np.asarray(ds.variables["coord" + L][:], dtype=float)
                for L in letters]
    if "coor_names" in ds.variables:
        names = decode_name_rows(ds.variables["coor_names"])
    else:
        names = letters
    return ndim, arrs, names


def snap(a, sig):
    """Snap an axis onto an absolute grid of cell = span * 10**(-sig).

    An absolute (per-axis) grid is deliberate, not a significant-figure /
    relative tolerance:

      * Coincident nodes always merge. Front/back node pairs are a pure
        translation, so their kept coordinates agree to round-off. On a body
        centreline one partner is often an exact 0.0 while the other is
        round-off "dust" (e.g. 3e-18); a relative tolerance treats those as
        wildly different and fails to merge them, cracking the mesh. Snapping
        to an absolute cell floors that dust onto the same grid point.
      * Distinct nodes stay apart. With sig=10 the cell is 1e-10 of the axis
        extent, far finer than any real mesh spacing (e.g. an airfoil boundary
        layer), yet far coarser than floating-point noise.

    Because the snap is to a fixed grid, every point in a cell maps to the
    identical float, so np.unique groups them exactly and the output can never
    contain two coincident-but-separate nodes.
    """
    a = np.asarray(a, dtype=float)
    span = float(a.max() - a.min())
    if span <= 0.0:
        return a.copy()
    cell = span * 10.0 ** (-sig)
    return np.round(a / cell) * cell


def axis_is_two_layers(a):
    """True if axis values fall into exactly two thin clusters (a single layer).

    Robust to floating-point "dust": a translation-extruded mesh can carry
    round-off on the depth coordinate (e.g. the top plane written as
    0.001 +/- 1e-9), so the layer is two *clusters*, not two exact values.
    The axis qualifies when both a low and a high cluster exist, each cluster's
    internal spread is negligible (< 1% of the gap between them), i.e. the
    nodes really do live on two parallel planes.
    """
    a = np.asarray(a, dtype=float)
    lo, hi = float(a.min()), float(a.max())
    span = hi - lo
    if span <= 0.0:
        return False
    mid = 0.5 * (lo + hi)
    near_lo = a[a <= mid]
    near_hi = a[a > mid]
    if near_lo.size == 0 or near_hi.size == 0:
        return False
    spread = max(near_lo.max() - near_lo.min(), near_hi.max() - near_hi.min())
    gap = near_hi.min() - near_lo.max()
    return gap > 0.0 and spread < 0.01 * gap


def detect_depth_axis(snapped):
    """Index of the single-layer (two-plane) axis; -1 if none.

    If several qualify, the one with the smallest extent (thinnest layer) wins.
    """
    candidates = [(float(a.max() - a.min()), i)
                  for i, a in enumerate(snapped) if axis_is_two_layers(a)]
    if not candidates:
        return -1
    candidates.sort()
    return candidates[0][1]


# ---------------------------------------------------------------------------
# core conversion
# ---------------------------------------------------------------------------
def convert(in_path, out_path, forced_axis=None, sig=10, dry_run=False):
    ds = netCDF4.Dataset(in_path, "r")

    ndim, coords, coor_names = read_coords(ds)
    if ndim != 3:
        ds.close()
        raise SystemExit(f"  SKIP: num_dim={ndim}, not a 3D mesh.")

    snapped = [snap(c, sig) for c in coords]
    depth = ({"x": 0, "y": 1, "z": 2}[forced_axis] if forced_axis
             else detect_depth_axis(snapped))
    if depth < 0:
        ds.close()
        raise SystemExit("  SKIP: no axis resolves into two parallel planes; "
                         "this mesh is not quasi-2D. (Try --depth-axis.)")

    keep = [i for i in range(3) if i != depth]
    depth_letter = ["x", "y", "z"][depth]
    keep_names = [coor_names[i] for i in keep]
    ka, kb = coords[keep[0]], coords[keep[1]]
    sdepth = snapped[depth]

    # --- node de-duplication: collapse each depth edge to one planar node ---
    # Front/back node pairs are merged TOPOLOGICALLY -- by the vertical element
    # edges that join the two depth planes -- never by (x, y) coincidence.
    #
    # A coordinate merge is fragile in two independent ways:
    #   1. Warped source mesh. If the two planes are not a clean translation
    #      (their front/back (x, y) differ by far more than round-off -- common
    #      in meshes imported from Fluent etc.), a fine-grid coordinate merge
    #      never pairs those nodes: the back node is left ORPHAN (referenced by
    #      no 2D element) and every side-set face built on it is dropped. This
    #      is exactly the bug this rewrite removes.
    #   2. Non-conformal / DG interface. The two sides of such an interface are
    #      geometrically coincident but are *distinct* node ids in different
    #      blocks, and that separation is what the interface requires (each face
    #      owned by its own element so its area vector points outward). A
    #      coordinate merge welds the sides -> normals align -> per-side fluxes
    #      add instead of cancel (~50% imbalance, unphysical fields).
    #
    # The topological merge handles both for free: a depth edge always lies
    # inside one element (hence one block), so front/back pairs collapse while
    # the two sides of an interface -- which share no element and no depth edge
    # -- stay separate. The planar (x, y) is taken from each column's front
    # (low-depth) node as the reference plane.
    n_old_nodes = len(ka)
    dmid = 0.5 * (float(sdepth.min()) + float(sdepth.max()))
    plane = (sdepth > dmid).astype(np.int8)        # 0 = front (low), 1 = back

    # gather the unique depth-edge node pairs across all blocks
    pair_chunks = []
    _b = 1
    while f"connect{_b}" in ds.variables:
        conn0 = np.asarray(ds.variables[f"connect{_b}"][:], dtype=np.int64) - 1
        canon = SOURCE_INFO[conn0.shape[1]][0]
        for i, j in ELEM_EDGES[canon]:
            a, b = conn0[:, i], conn0[:, j]
            cross = plane[a] != plane[b]
            if cross.any():
                pair_chunks.append(np.stack([a[cross], b[cross]], axis=1))
        _b += 1
    pairs = (np.unique(np.sort(np.vstack(pair_chunks), axis=1), axis=0)
             if pair_chunks else np.empty((0, 2), dtype=np.int64))

    # union-find over old nodes, joining each depth-edge pair
    parent = np.arange(n_old_nodes, dtype=np.int64)

    def _find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:                   # path compression
            parent[i], i = root, parent[i]
        return root

    for u, v in pairs:
        ru, rv = _find(int(u)), _find(int(v))
        if ru != rv:
            parent[ru] = rv

    roots = np.array([_find(i) for i in range(n_old_nodes)], dtype=np.int64)
    _, node_map = np.unique(roots, return_inverse=True)      # old0 -> new0
    node_map = np.asarray(node_map, dtype=np.int64).ravel()
    n_new_nodes = int(node_map.max()) + 1

    # representative planar coordinate = the front (low-depth) node of each
    # column; assign back nodes first so the front node wins the slot.
    order = np.argsort(-plane, kind="stable")                # back (1) first
    new_x = np.empty(n_new_nodes); new_x[node_map[order]] = ka[order]
    new_y = np.empty(n_new_nodes); new_y[node_map[order]] = kb[order]

    # --- element blocks: collapse each cell to its planar cap face ----------
    eb_names_src = (decode_name_rows(ds.variables["eb_names"])
                    if "eb_names" in ds.variables else None)
    blocks = []          # one dict per element block
    belong = []          # global elem (0-based) -> block index
    local0 = []          # global elem (0-based) -> local index in its block
    blk = 1
    n_collapsed = 0
    while f"connect{blk}" in ds.variables:
        cvar = ds.variables[f"connect{blk}"]
        conn = np.asarray(cvar[:], dtype=np.int64)            # (nel, nnode) 1-based
        nel, nnode = conn.shape
        src_type = cvar.getncattr("elem_type").upper()
        if nnode not in SOURCE_INFO:
            ds.close()
            raise SystemExit(f"  Unsupported element type {src_type!r} "
                             f"({nnode} nodes/elem) in block {blk}.")
        canon, type2d, n2d = SOURCE_INFO[nnode]

        # plane label (0 = front/low, 1 = back/high) of every node of every
        # element. For each element pick its front cap face -- the side whose
        # nodes all lie in the SAME depth plane, preferring the front (label 0).
        #
        # Classification is by plane MEMBERSHIP, not by exact snapped-depth
        # equality. A plane is rarely perfectly flat: translation-extruded
        # meshes carry round-off "dust" so a single plane is written as 0 mixed
        # with 1e-8, or 0.005 mixed with 0.005001. That dust is negligible
        # against the front/back gap (axis_is_two_layers already vouched for
        # this) but can exceed the snap-grid cell, in which case an exact-depth
        # cap test splits an otherwise-flat face and leaves the element with no
        # cap. Bucketing each node into front/back by the midpoint is immune to
        # within-plane dust at any --sig-figs.
        pl = plane[conn - 1]                                   # (nel, nnode) 0/1
        best_side = np.zeros(nel, dtype=np.int64)
        best_depth = np.full(nel, 2, dtype=np.int64)           # > max plane label
        for s, locs in SIDE_NODES[canon].items():
            idx = np.asarray(locs) - 1
            sub = pl[:, idx]
            is_cap = sub.max(axis=1) == sub.min(axis=1)
            cap_d = sub[:, 0].astype(np.int64)                 # 0 front, 1 back
            better = is_cap & (cap_d < best_depth)
            best_side[better] = s
            best_depth[better] = cap_d[better]
        if np.any(best_side == 0):
            ds.close()
            raise SystemExit(f"  Block {blk}: {int((best_side == 0).sum())} "
                             f"element(s) have no depth-normal face; the mesh is "
                             f"not a single layer along '{depth_letter}'.")

        new_conn = np.empty((nel, n2d), dtype=np.int64)
        for s, locs in SIDE_NODES[canon].items():
            sel = best_side == s
            if not sel.any():
                continue
            idx = np.asarray(locs) - 1
            new_conn[sel] = node_map[conn[sel][:, idx] - 1] + 1

        # enforce counter-clockwise winding (positive area)
        px = new_x[new_conn - 1]
        py = new_y[new_conn - 1]
        area = 0.5 * np.sum(px * np.roll(py, -1, axis=1)
                            - np.roll(px, -1, axis=1) * py, axis=1)
        flip = area < 0.0
        new_conn[flip] = new_conn[flip][:, ::-1]

        # detect any element that collapsed to fewer than n2d distinct nodes
        srt = np.sort(new_conn, axis=1)
        dup = np.any(srt[:, 1:] == srt[:, :-1], axis=1)
        n_collapsed += int(dup.sum())

        gstart = sum(b["nel"] for b in blocks)
        belong.append(np.full(nel, len(blocks), dtype=np.int64))
        local0.append(np.arange(nel, dtype=np.int64))
        _bid = int(np.asarray(ds.variables["eb_prop1"][:])[blk - 1])
        _bname = (eb_names_src[blk - 1]
                  if eb_names_src and eb_names_src[blk - 1] else f"block_{_bid}")
        blocks.append({
            "id": _bid, "name": _bname,
            "type2d": type2d, "n2d": n2d, "canon": canon, "nel": nel,
            "gstart": gstart, "conn_old": conn, "conn_new": new_conn,
        })
        blk += 1
    n_elem = sum(b["nel"] for b in blocks)
    belong = np.concatenate(belong)
    local0 = np.concatenate(local0)

    if n_collapsed:
        ds.close()
        raise SystemExit(
            f"  {n_collapsed} element(s) collapsed to a degenerate 2D cell.\n"
            f"  Distinct nodes were merged: increase precision with a larger\n"
            f"  --sig-figs (currently {sig}).")

    # --- side sets: lateral faces -> edges; cap faces -> dropped -----------
    n_ss = (ds.dimensions["num_side_sets"].size
            if "num_side_sets" in ds.dimensions else 0)
    ss_ids = np.asarray(ds.variables["ss_prop1"][:]) if n_ss else []
    ss_names = (decode_name_rows(ds.variables["ss_names"])
                if "ss_names" in ds.variables else None)
    kept_ss, dropped_ss = [], []
    lost_faces = {}            # side-set name -> count of edge faces left untagged
    for s in range(n_ss):
        g0 = np.asarray(ds.variables[f"elem_ss{s + 1}"][:], dtype=np.int64) - 1
        sides = np.asarray(ds.variables[f"side_ss{s + 1}"][:], dtype=np.int64)
        out_e, out_s = [], []
        n_lost = 0
        bsel = belong[g0]
        for bi in np.unique(bsel):
            b = blocks[bi]
            for sd in np.unique(sides[bsel == bi]):
                m = (bsel == bi) & (sides == sd)
                loc = local0[g0[m]]
                locs = np.asarray(SIDE_NODES[b["canon"]][int(sd)]) - 1
                face2d = node_map[b["conn_old"][loc][:, locs] - 1] + 1
                # distinct-node count per face row
                srt = np.sort(face2d, axis=1)
                ndist = 1 + np.sum(srt[:, 1:] != srt[:, :-1], axis=1)
                edge = ndist == 2
                if not edge.any():
                    continue                       # whole group is a cap face
                fe = face2d[edge]
                A = fe[:, 0]
                B = fe[np.arange(len(fe)), np.argmax(fe != A[:, None], axis=1)]
                nc = b["conn_new"][loc[edge]]
                snew = np.zeros(len(nc), dtype=np.int64)
                for k, (i, j) in enumerate(EDGE_IDX[b["type2d"]]):
                    hit = (((nc[:, i] == A) & (nc[:, j] == B))
                           | ((nc[:, i] == B) & (nc[:, j] == A)))
                    snew[hit] = k + 1
                ok = snew > 0
                n_lost += int((~ok).sum())     # edge face matching no quad edge
                ge = (g0[m][edge] + 1)[ok]
                out_e.append(ge)
                out_s.append(snew[ok])
        if out_e:
            kept_ss.append((int(ss_ids[s]),
                            ss_names[s] if ss_names else "",
                            np.concatenate(out_e), np.concatenate(out_s)))
        else:
            dropped_ss.append((ss_names[s] if ss_names else f"#{int(ss_ids[s])}"))
        if n_lost:
            lost_faces[ss_names[s] if ss_names else f"#{int(ss_ids[s])}"] = n_lost

    if lost_faces:
        ds.close()
        raise SystemExit(
            "  side-set conversion lost lateral faces that matched no 2D edge "
            f"(a topology/merge inconsistency): {lost_faces}.\n"
            "  This must not happen with the topological node merge -- the mesh "
            "may not be a clean single layer along the depth axis.")

    # --- node sets ---------------------------------------------------------
    n_ns = (ds.dimensions["num_node_sets"].size
            if "num_node_sets" in ds.dimensions else 0)
    ns_ids = np.asarray(ds.variables["ns_prop1"][:]) if n_ns else []
    ns_names = (decode_name_rows(ds.variables["ns_names"])
                if "ns_names" in ds.variables else None)
    kept_ns = []
    for s in range(n_ns):
        nodes0 = np.asarray(ds.variables[f"node_ns{s + 1}"][:], dtype=np.int64) - 1
        remapped = np.unique(node_map[nodes0] + 1)
        kept_ns.append((int(ns_ids[s]),
                        ns_names[s] if ns_names else "", remapped))

    # --- safety: every output node must be used by some 2D element ---------
    referenced = np.zeros(n_new_nodes, dtype=bool)
    for b in blocks:
        referenced[b["conn_new"].ravel() - 1] = True
    n_orphan = int((~referenced).sum())
    if n_orphan:
        ds.close()
        raise SystemExit(
            f"  {n_orphan} of {n_new_nodes} output nodes are orphan (used by no "
            "element).\n  The front/back planes did not collapse cleanly -- the "
            "input is not a single layer along the depth axis.")

    title = ds.getncattr("title") if "title" in ds.ncattrs() else "mesh"
    ds.close()

    # --- report ------------------------------------------------------------
    print(f"  depth axis : {depth_letter}  (kept planar axes: "
          f"{keep_names[0]}, {keep_names[1]})  [sig-figs={sig}]")
    print(f"  nodes      : {len(node_map)} -> {n_new_nodes}")
    print(f"  elements   : {n_elem} (" +
          ", ".join(f"{b['type2d']}x{b['nel']}" for b in blocks) + ")")
    print(f"  side sets  : kept {[k[1] for k in kept_ss]}")
    if dropped_ss:
        print(f"               dropped (front/back) {dropped_ss}")
    if kept_ns:
        print(f"  node sets  : {[k[1] for k in kept_ns]}")
    if dry_run:
        print("  (dry-run: nothing written)")
        return

    write_exodus(out_path, title, keep_names, new_x, new_y,
                 blocks, kept_ss, kept_ns)
    print(f"  written    : {out_path}")


def write_exodus(out_path, title, keep_names, new_x, new_y,
                 blocks, kept_ss, kept_ns):
    out = netCDF4.Dataset(out_path, "w", format="NETCDF3_64BIT_OFFSET")
    out.api_version = np.float32(7.03)
    out.version = np.float32(7.03)
    out.floating_point_word_size = np.int32(8)
    out.file_size = np.int32(1)
    out.maximum_name_length = np.int32(32)
    out.int64_status = np.int32(0)
    out.title = (str(title) + " (literal 2D)")[:80]

    out.createDimension("len_string", 33)
    out.createDimension("len_line", 81)
    out.createDimension("four", 4)
    out.createDimension("len_name", 256)
    out.createDimension("time_step", None)
    out.createDimension("num_dim", 2)
    out.createDimension("num_nodes", len(new_x))
    out.createDimension("num_elem", sum(b["nel"] for b in blocks))
    out.createDimension("num_el_blk", len(blocks))

    out.createVariable("time_whole", "f8", ("time_step",))
    out.createVariable("eb_status", "i4", ("num_el_blk",))[:] = \
        np.ones(len(blocks), dtype=np.int32)
    ebp = out.createVariable("eb_prop1", "i4", ("num_el_blk",))
    ebp.setncattr("name", "ID")
    ebp[:] = np.array([b["id"] for b in blocks], dtype=np.int32)

    out.createVariable("coordx", "f8", ("num_nodes",))[:] = new_x
    out.createVariable("coordy", "f8", ("num_nodes",))[:] = new_y

    write_name_rows(out.createVariable("eb_names", "S1",
                    ("num_el_blk", "len_name")),
                    [b["name"] for b in blocks])
    write_name_rows(out.createVariable("coor_names", "S1",
                    ("num_dim", "len_name")), keep_names)

    for i, b in enumerate(blocks, start=1):
        out.createDimension(f"num_el_in_blk{i}", b["nel"])
        out.createDimension(f"num_nod_per_el{i}", b["n2d"])
        cv = out.createVariable(f"connect{i}", "i4",
                                (f"num_el_in_blk{i}", f"num_nod_per_el{i}"))
        cv.elem_type = b["type2d"]
        cv[:] = b["conn_new"].astype(np.int32)

    if kept_ss:
        out.createDimension("num_side_sets", len(kept_ss))
        out.createVariable("ss_status", "i4", ("num_side_sets",))[:] = \
            np.ones(len(kept_ss), dtype=np.int32)
        ssp = out.createVariable("ss_prop1", "i4", ("num_side_sets",))
        ssp.setncattr("name", "ID")
        ssp[:] = np.array([k[0] for k in kept_ss], dtype=np.int32)
        write_name_rows(out.createVariable("ss_names", "S1",
                        ("num_side_sets", "len_name")),
                        [k[1] for k in kept_ss])
        for i, (sid, name, e, sd) in enumerate(kept_ss, start=1):
            out.createDimension(f"num_side_ss{i}", len(e))
            out.createVariable(f"elem_ss{i}", "i4", (f"num_side_ss{i}",))[:] = \
                e.astype(np.int32)
            out.createVariable(f"side_ss{i}", "i4", (f"num_side_ss{i}",))[:] = \
                sd.astype(np.int32)

    if kept_ns:
        out.createDimension("num_node_sets", len(kept_ns))
        out.createVariable("ns_status", "i4", ("num_node_sets",))[:] = \
            np.ones(len(kept_ns), dtype=np.int32)
        nsp = out.createVariable("ns_prop1", "i4", ("num_node_sets",))
        nsp.setncattr("name", "ID")
        nsp[:] = np.array([k[0] for k in kept_ns], dtype=np.int32)
        write_name_rows(out.createVariable("ns_names", "S1",
                        ("num_node_sets", "len_name")),
                        [k[1] for k in kept_ns])
        for i, (nid, name, nodes) in enumerate(kept_ns, start=1):
            out.createDimension(f"num_nod_ns{i}", len(nodes))
            out.createVariable(f"node_ns{i}", "i4", (f"num_nod_ns{i}",))[:] = \
                nodes.astype(np.int32)

    out.close()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="input quasi-2D Exodus mesh (.e)")
    ap.add_argument("output", nargs="?", help="output 2D mesh (default <stem>_2D.e)")
    ap.add_argument("-o", "--out", dest="out_opt", help="alternative to positional output")
    ap.add_argument("--depth-axis", choices=["x", "y", "z"],
                    help="force the depth axis instead of auto-detecting")
    ap.add_argument("--sig-figs", type=int, default=10,
                    help="node-merge grid resolution: planar nodes are snapped "
                         "to a grid of (axis extent) * 10**(-sig-figs). Lower it "
                         "if distinct nodes merge (degenerate cells); raise it "
                         "if coincident nodes fail to merge. Default 10.")
    ap.add_argument("--dry-run", action="store_true",
                    help="analyse and report, but do not write a file")
    args = ap.parse_args()

    out_path = args.output or args.out_opt
    if out_path is None:
        stem, _ = os.path.splitext(args.input)
        out_path = stem + "_2D.e"

    print(f"{args.input}")
    convert(args.input, out_path, forced_axis=args.depth_axis,
            sig=args.sig_figs, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
