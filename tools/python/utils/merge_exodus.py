#!/usr/bin/env python3
# File       : merge_exodus.py
# Created    : Fri Sep 11 2026
# Author     : Mhamad Mahdi Alloush (drafted with Claude Code)
# Description: Merge several Exodus II meshes into one multi-block file
# Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.
"""
Merge any number of Exodus II meshes into a single file WITHOUT merging nodes:
every input keeps its own nodes and elements, which are appended with a global
offset. Overlapping or touching meshes are allowed (a rotor/stator pair stays
non-conformal at the common surface); coincident nodes are only reported.

Element blocks, side sets and node sets are carried over. IDs are re-assigned
sequentially (1..N) in order of appearance. Names are kept unless renamed and
MUST be unique across the merged file: the merge aborts with the list of
duplicates otherwise. Original id maps and distribution factors are dropped.

Usage (options after an --input apply to that input only):

    python3 merge_exodus.py -o mesh.e [--scale 0.001] [--title TEXT] \
        --input rotor.exo  --block rotor  --rename WALL_INT=ROTOR_WALL_INT ... \
        --input stator.exo --block stator --rename WALL_INT=STATOR_WALL_INT ... \
        --input casing.exo --prefix CASING_

    --block NEW          rename the (single) element block of this input
    --block OLD=NEW      rename one block of a multi-block input
    --rename OLD=NEW     rename a side set or node set of this input
    --prefix P           prepend P to every block/side-set/node-set name of this
                         input that is not renamed explicitly
    --scale S            multiply all coordinates by S (e.g. mm -> m)
    --tol D              coincident-node report tolerance (after scaling)
    --no-shared-check    skip the pairwise coincident-node report
"""
import os
import sys
import time
from datetime import datetime

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree

ELEM_TYPE = {(3, 8): "HEX8", (3, 4): "TET4", (3, 6): "WEDGE6", (3, 5): "PYRAMID5",
             (2, 4): "QUAD4", (2, 3): "TRI3", (3, 20): "HEX20", (3, 10): "TET10"}


def _arr(v):
    return np.ma.getdata(v[:])


def decode_names(var):
    out = []
    for row in np.ma.filled(var[:], b"\x00"):
        raw = b"".join(bytes(c) for c in row.tolist())
        out.append(raw.split(b"\x00", 1)[0].decode("ascii", "ignore").strip())
    return out


def write_names(var, names):
    nrow, ncol = var.shape
    buf = np.zeros((nrow, ncol), dtype="S1")
    for i, name in enumerate(names):
        for j, ch in enumerate(name.encode("ascii", "ignore")[: ncol - 1]):
            buf[i, j] = bytes([ch])
    var[:] = buf


def read_coords(ds):
    nd = ds.dimensions["num_dim"].size
    if "coord" in ds.variables:
        return np.asarray(_arr(ds.variables["coord"]), float).T[:, :nd]
    return np.stack([np.asarray(_arr(ds.variables[c]), float)
                     for c in ("coordx", "coordy", "coordz")[:nd]], axis=1)


def parse_args(argv):
    g = {"out": None, "scale": 1.0, "title": None, "tol": 1e-9, "inputs": [],
         "shared_check": True}
    cur = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-o", "--out"):
            g["out"] = argv[i + 1]; i += 2
        elif a == "--scale":
            g["scale"] = float(argv[i + 1]); i += 2
        elif a == "--title":
            g["title"] = argv[i + 1]; i += 2
        elif a == "--tol":
            g["tol"] = float(argv[i + 1]); i += 2
        elif a == "--no-shared-check":
            g["shared_check"] = False; i += 1
        elif a == "--input":
            cur = {"file": argv[i + 1], "blocks": [], "renames": {}, "prefix": ""}
            g["inputs"].append(cur); i += 2
        elif a == "--prefix":
            if cur is None:
                sys.exit("--prefix must follow an --input")
            cur["prefix"] = argv[i + 1]; i += 2
        elif a == "--block":
            if cur is None:
                sys.exit("--block must follow an --input")
            cur["blocks"].append(argv[i + 1]); i += 2
        elif a == "--rename":
            if cur is None:
                sys.exit("--rename must follow an --input")
            old, new = argv[i + 1].split("=", 1)
            cur["renames"][old] = new; i += 2
        elif a in ("-h", "--help"):
            print(__doc__); sys.exit(0)
        else:
            sys.exit(f"unknown argument {a!r}")
    if g["out"] is None or not g["inputs"]:
        sys.exit("need -o OUT and at least one --input file")
    return g


def load_input(spec, scale):
    ds = Dataset(spec["file"])
    nd = ds.dimensions["num_dim"].size
    X = read_coords(ds) * scale
    nb = ds.dimensions["num_el_blk"].size
    bnames = decode_names(ds.variables["eb_names"]) if "eb_names" in ds.variables else [""] * nb
    bids = _arr(ds.variables["eb_prop1"])
    # block renames: "NEW" for a single block, "OLD=NEW" otherwise
    bmap = {}
    for b in spec["blocks"]:
        if "=" in b:
            old, new = b.split("=", 1); bmap[old] = new
        elif nb == 1:
            bmap[bnames[0]] = b
        else:
            sys.exit(f"{spec['file']}: --block NEW needs a single-block file, use OLD=NEW")
    pre = spec["prefix"]
    blocks = []
    for k in range(nb):
        cv = ds.variables[f"connect{k + 1}"]
        conn = _arr(cv).astype(np.int64)
        etype = ELEM_TYPE.get((nd, conn.shape[1]), str(getattr(cv, "elem_type", "")).upper())
        blocks.append({"name": bmap.get(bnames[k], pre + bnames[k]), "src_name": bnames[k],
                       "src_id": int(bids[k]), "conn": conn, "type": etype})
    used = spec["renames"]
    sidesets = []
    if "num_side_sets" in ds.dimensions:
        n = ds.dimensions["num_side_sets"].size
        names = decode_names(ds.variables["ss_names"]) if "ss_names" in ds.variables else [""] * n
        ids = _arr(ds.variables["ss_prop1"])
        for k in range(n):
            sidesets.append({"name": used.get(names[k], pre + names[k]), "src_name": names[k],
                             "src_id": int(ids[k]),
                             "elem": _arr(ds.variables[f"elem_ss{k + 1}"]).astype(np.int64),
                             "side": _arr(ds.variables[f"side_ss{k + 1}"]).astype(np.int64)})
    nodesets = []
    if "num_node_sets" in ds.dimensions:
        n = ds.dimensions["num_node_sets"].size
        names = decode_names(ds.variables["ns_names"]) if "ns_names" in ds.variables else [""] * n
        ids = _arr(ds.variables["ns_prop1"])
        for k in range(n):
            nodesets.append({"name": used.get(names[k], pre + names[k]), "src_name": names[k],
                             "src_id": int(ids[k]),
                             "nodes": _arr(ds.variables[f"node_ns{k + 1}"]).astype(np.int64)})
    unknown = set(used) - {s["src_name"] for s in sidesets + nodesets}
    if unknown:
        sys.exit(f"{spec['file']}: --rename targets not found: {sorted(unknown)}")
    title = str(getattr(ds, "title", ""))
    ds.close()
    nelem = sum(b["conn"].shape[0] for b in blocks)
    return {"file": spec["file"], "nd": nd, "X": X, "blocks": blocks, "sidesets": sidesets,
            "nodesets": nodesets, "nelem": nelem, "title": title}


def check_shared_nodes(inputs, tol):
    """Report coincident nodes between inputs (they are NOT merged)."""
    for a in range(len(inputs)):
        tree = cKDTree(inputs[a]["X"])
        for b in range(a + 1, len(inputs)):
            d, _ = tree.query(inputs[b]["X"])
            print(f"  shared-node check {inputs[a]['file']} <-> {inputs[b]['file']}: "
                  f"min dist {d.min():.3e}, {(d < tol).sum()} nodes within tol {tol:g}")


def main():
    g = parse_args(sys.argv[1:])
    inputs = [load_input(s, g["scale"]) for s in g["inputs"]]
    nd = {i["nd"] for i in inputs}
    if len(nd) != 1:
        sys.exit(f"spatial dimensions differ: {nd}")
    nd = nd.pop()
    if g["shared_check"] and len(inputs) > 1:
        check_shared_nodes(inputs, g["tol"])

    # global offsets: nodes and elements are appended in input order
    node_off, elem_off = [], []
    nn = ne = 0
    for inp in inputs:
        node_off.append(nn); elem_off.append(ne)
        nn += inp["X"].shape[0]; ne += inp["nelem"]

    blocks, sidesets, nodesets, info = [], [], [], []
    for k, inp in enumerate(inputs):
        src = os.path.basename(inp["file"])
        for b in inp["blocks"]:
            blocks.append({**b, "conn": b["conn"] + node_off[k], "id": len(blocks) + 1})
            info.append(f"block {blocks[-1]['id']} '{b['name']}' <- {src} "
                        f"block '{b['src_name']}' (id {b['src_id']}), {b['conn'].shape[0]} {b['type']}")
        for s in inp["sidesets"]:
            sidesets.append({**s, "elem": s["elem"] + elem_off[k], "id": len(sidesets) + 1})
            info.append(f"sideset {sidesets[-1]['id']} '{s['name']}' <- {src} "
                        f"sideset '{s['src_name']}' (id {s['src_id']}), {len(s['elem'])} faces")
        for s in inp["nodesets"]:
            nodesets.append({**s, "nodes": s["nodes"] + node_off[k], "id": len(nodesets) + 1})
            info.append(f"nodeset {nodesets[-1]['id']} '{s['name']}' <- {src} "
                        f"nodeset '{s['src_name']}' (id {s['src_id']}), {len(s['nodes'])} nodes")
    for kind, items in (("block", blocks), ("side set", sidesets), ("node set", nodesets)):
        names = [x["name"] for x in items]
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            sys.exit(f"duplicate {kind} names after merge: {dup}; use --block/--rename/--prefix")

    X = np.vstack([inp["X"] for inp in inputs])
    if nn > np.iinfo(np.int32).max or ne > np.iinfo(np.int32).max:
        sys.exit("merged mesh exceeds 32-bit ids")

    title = g["title"] or ("merged: " + " + ".join(i["file"] for i in inputs))
    out = Dataset(g["out"], "w", format="NETCDF3_64BIT_OFFSET")
    out.api_version = np.float32(7.03)
    out.version = np.float32(7.03)
    out.floating_point_word_size = np.int32(8)
    out.file_size = np.int32(1)
    out.maximum_name_length = np.int32(32)
    out.int64_status = np.int32(0)
    out.title = title[:80]
    out.createDimension("len_string", 33)
    out.createDimension("len_line", 81)
    out.createDimension("four", 4)
    out.createDimension("len_name", 33)
    out.createDimension("time_step", None)
    out.createDimension("num_dim", nd)
    out.createDimension("num_nodes", nn)
    out.createDimension("num_elem", ne)
    out.createDimension("num_el_blk", len(blocks))
    out.createVariable("time_whole", "f8", ("time_step",))
    out.createVariable("eb_status", "i4", ("num_el_blk",))[:] = np.ones(len(blocks), np.int32)
    ebp = out.createVariable("eb_prop1", "i4", ("num_el_blk",))
    ebp.setncattr("name", "ID")
    ebp[:] = np.array([b["id"] for b in blocks], np.int32)
    for c, name in zip(("coordx", "coordy", "coordz"), ("x", "y", "z")):
        if c[-1] in "xyz"[:nd]:
            out.createVariable(c, "f8", ("num_nodes",))[:] = X[:, "xyz".index(name)]
    write_names(out.createVariable("eb_names", "S1", ("num_el_blk", "len_name")),
                [b["name"] for b in blocks])
    write_names(out.createVariable("coor_names", "S1", ("num_dim", "len_name")), list("xyz"[:nd]))
    for i, b in enumerate(blocks, start=1):
        out.createDimension(f"num_el_in_blk{i}", b["conn"].shape[0])
        out.createDimension(f"num_nod_per_el{i}", b["conn"].shape[1])
        cv = out.createVariable(f"connect{i}", "i4", (f"num_el_in_blk{i}", f"num_nod_per_el{i}"))
        cv.elem_type = b["type"]
        cv[:] = b["conn"].astype(np.int32)
    if sidesets:
        out.createDimension("num_side_sets", len(sidesets))
        out.createVariable("ss_status", "i4", ("num_side_sets",))[:] = np.ones(len(sidesets), np.int32)
        ssp = out.createVariable("ss_prop1", "i4", ("num_side_sets",))
        ssp.setncattr("name", "ID")
        ssp[:] = np.array([s["id"] for s in sidesets], np.int32)
        write_names(out.createVariable("ss_names", "S1", ("num_side_sets", "len_name")),
                    [s["name"] for s in sidesets])
        for i, s in enumerate(sidesets, start=1):
            out.createDimension(f"num_side_ss{i}", len(s["elem"]))
            out.createVariable(f"elem_ss{i}", "i4", (f"num_side_ss{i}",))[:] = s["elem"].astype(np.int32)
            out.createVariable(f"side_ss{i}", "i4", (f"num_side_ss{i}",))[:] = s["side"].astype(np.int32)
    if nodesets:
        out.createDimension("num_node_sets", len(nodesets))
        out.createVariable("ns_status", "i4", ("num_node_sets",))[:] = np.ones(len(nodesets), np.int32)
        nsp = out.createVariable("ns_prop1", "i4", ("num_node_sets",))
        nsp.setncattr("name", "ID")
        nsp[:] = np.array([s["id"] for s in nodesets], np.int32)
        write_names(out.createVariable("ns_names", "S1", ("num_node_sets", "len_name")),
                    [s["name"] for s in nodesets])
        for i, s in enumerate(nodesets, start=1):
            out.createDimension(f"num_nod_ns{i}", len(s["nodes"]))
            out.createVariable(f"node_ns{i}", "i4", (f"num_nod_ns{i}",))[:] = s["nodes"].astype(np.int32)
    # provenance: QA record + one info line per merged entity
    now = datetime.now()
    out.createDimension("num_qa_rec", 1)
    qa = out.createVariable("qa_records", "S1", ("num_qa_rec", "four", "len_string"))
    buf = np.zeros((1, 4, 33), dtype="S1")
    for j, s in enumerate(["merge_exodus.py", "1.0", now.strftime("%Y-%m-%d"), now.strftime("%H:%M:%S")]):
        for c, ch in enumerate(s.encode()[:32]):
            buf[0, j, c] = bytes([ch])
    qa[:] = buf
    info = [f"coordinate scale factor {g['scale']:g}"] + info
    out.createDimension("num_info", len(info))
    write_names(out.createVariable("info_records", "S1", ("num_info", "len_line")), info)
    out.close()

    print(f"wrote {g['out']}: dim {nd}, {nn} nodes, {ne} elements, "
          f"{len(blocks)} blocks, {len(sidesets)} side sets, {len(nodesets)} node sets")
    for line in info:
        print("  " + line)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"done in {time.time() - t0:.1f}s")
