#!/usr/bin/env python3
# File       : test_exodus_to_cgns.py
# Description: Test exodus_to_cgns.py with synthetic 2D and 3D Exodus meshes.
#
# Builds small Exodus files, converts them, and validates the CGNS output
# with h5py: base/zone dimensions, element section types, connectivity,
# and boundary (sideset) sections. cgnscheck (CGNS tools) validates the
# files when installed.
#
#   usage:  test_exodus_to_cgns.py [path/to/exodus_to_cgns.py]

import os
import shutil
import subprocess
import sys
import tempfile

import h5py
import netCDF4
import numpy as np


def s(key):
    """h5py group key -> str (bytes on some h5py versions)."""
    return key.decode() if isinstance(key, bytes) else key


# --- CGNS element types (cgnslib.h) ----------------------------------------
BAR_2  = 3
TRI_3  = 5
QUAD_4 = 7
HEXA_8 = 17


def write_exodus(path, coords, blocks, sidesets):
    """Write a minimal Exodus file.

    coords   : (nnode, dim) node coordinates
    blocks   : [(name, (nelem, npe) 1-based connectivity)]
    sidesets : [(name, 1-based global element ids, 1-based side ids)]
    """
    nnode, dim = coords.shape
    nelem = sum(b[1].shape[0] for b in blocks)

    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("num_dim", dim)
        d.createDimension("num_nodes", nnode)
        d.createDimension("num_elem", nelem)
        d.createDimension("num_el_blk", len(blocks))
        d.createDimension("len_name", 256)

        for i, c in enumerate("xyz"[:dim]):
            d.createVariable(f"coord{c}", "f8", ("num_nodes",))[:] = coords[:, i]

        # element blocks
        eb_names = np.zeros((len(blocks), 256), "S1")
        for i, (name, conn) in enumerate(blocks):
            d.createDimension(f"num_el_in_blk{i + 1}", conn.shape[0])
            d.createDimension(f"num_nod_per_el{i + 1}", conn.shape[1])
            d.createVariable(f"connect{i + 1}", "i8",
                             (f"num_el_in_blk{i + 1}",
                              f"num_nod_per_el{i + 1}"))[:] = conn
            eb_names[i, :len(name)] = np.frombuffer(name.encode(), "S1")
        d.createVariable("eb_names", "S1",
                         ("num_el_blk", "len_name"))[:] = eb_names

        # sidesets
        d.createDimension("num_side_sets", len(sidesets))
        ss_names = np.zeros((len(sidesets), 256), "S1")
        for k, (name, elems, sides) in enumerate(sidesets):
            d.createDimension(f"num_side_ss{k + 1}", len(elems))
            d.createVariable(f"elem_ss{k + 1}", "i8",
                             (f"num_side_ss{k + 1}",))[:] = elems
            d.createVariable(f"side_ss{k + 1}", "i8",
                             (f"num_side_ss{k + 1}",))[:] = sides
            ss_names[k, :len(name)] = np.frombuffer(name.encode(), "S1")
        d.createVariable("ss_names", "S1",
                         ("num_side_sets", "len_name"))[:] = ss_names


def convert(script, mesh, out):
    """Run exodus_to_cgns.py on mesh and return its stdout."""
    r = subprocess.run([sys.executable, script, mesh, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"conversion failed:\n{r.stdout}\n{r.stderr}"
    return r.stdout


def read_base(h5):
    """Return (cell_dim, phys_dim, {zone name: zone group}) of a CGNS base."""
    base = h5["Base"]
    dims = base[" data"][()].tolist()
    zones = {s(k): v for k, v in base.items()
             if isinstance(v, h5py.Group) and v.attrs["label"] == b"Zone_t"}
    return dims[0], dims[1], zones


def read_sections(zone):
    """Return {section name: (element type, range, flat connectivity)}."""
    out = {}
    for k, v in zone.items():
        if not isinstance(v, h5py.Group) or v.attrs["label"] != b"Elements_t":
            continue
        etype = int(v[" data"][0])
        rng = v["ElementRange/ data"][()].tolist()
        conn = v["ElementConnectivity/ data"][()].tolist()
        out[s(k)] = (etype, rng, conn)
    return out


def read_coords(zone):
    """Return {coordinate name: values} of a CGNS zone."""
    gc = zone["GridCoordinates"]
    return {s(k): v[" data"][()].tolist() for k, v in gc.items()
            if isinstance(v, h5py.Group)}


def check_section(secs, name, etype, rng, conn):
    assert name in secs, f"missing section '{name}' in {list(secs)}"
    assert secs[name] == (etype, rng, conn), \
        f"section '{name}':\n  got      {secs[name]}\n  expected {(etype, rng, conn)}"


def cgnscheck(path):
    """Run cgnscheck on path (skipped when not installed)."""
    if shutil.which("cgnscheck") is None:
        print("  cgnscheck not found, skipped")
        return
    r = subprocess.run(["cgnscheck", path], capture_output=True, text=True)
    assert r.returncode == 0, f"cgnscheck failed:\n{r.stdout}\n{r.stderr}"


def quad_grid(n0, x0, y0):
    """2x2 quad grid, 3x3 nodes; returns (coords, 1-based connectivity)."""
    coords = np.array([[x0 + i, y0 + j] for j in range(3) for i in range(3)])
    nid = lambda i, j: n0 + i + 3 * j          # 1-based node id
    conn = np.array([[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
                     for j in range(2) for i in range(2)])
    return coords, conn


def test_2d(script, tmp):
    """Two disjoint quad blocks with one sideset each."""
    coords_a, conn_a = quad_grid(1, 0.0, 0.0)        # nodes 1..9
    coords_b, conn_b = quad_grid(10, 5.0, 0.0)       # nodes 10..18
    coords = np.vstack([coords_a, coords_b])

    # side 1 = bottom edge (n1,n2); side 3 = top edge (n3,n4)
    blocks = [("blkA-QUAD", conn_a), ("blkB-QUAD", conn_b)]
    sidesets = [("inletA", [1, 2], [1, 1]),
                ("topB", [7, 8], [3, 3])]

    exo, cgns = f"{tmp}/mesh2d.exo", f"{tmp}/mesh2d.cgns"
    write_exodus(exo, coords, blocks, sidesets)
    convert(script, exo, cgns)

    with h5py.File(cgns) as h5:
        cell_dim, phys_dim, zones = read_base(h5)
        assert (cell_dim, phys_dim) == (2, 2), "base must be 2D"

        za, zb = zones["blkA-QUAD"], zones["blkB-QUAD"]

        # zones: node/cell counts, 2D coordinates only
        assert za[" data"][()].ravel().tolist() == [9, 4, 0]
        assert zb[" data"][()].ravel().tolist() == [9, 4, 0]
        ca, cb = read_coords(za), read_coords(zb)
        assert sorted(ca) == ["CoordinateX", "CoordinateY"], "2D: no Z coord"
        assert ca["CoordinateX"] == coords_a[:, 0].tolist()
        assert ca["CoordinateY"] == coords_a[:, 1].tolist()
        assert cb["CoordinateX"] == coords_b[:, 0].tolist()
        assert cb["CoordinateY"] == coords_b[:, 1].tolist()

        # volume sections: QUAD_4 with local (1-based) connectivity
        sa, sb = read_sections(za), read_sections(zb)
        check_section(sa, "blkA-QUAD_cells", QUAD_4, [1, 4], conn_a.ravel().tolist())
        check_section(sb, "blkB-QUAD_cells", QUAD_4, [1, 4],
                      (conn_b - 9).ravel().tolist())

        # boundary sections: BAR_2 edges, numbered after the cells
        check_section(sa, "inletA", BAR_2, [5, 6], [1, 2, 2, 3])
        check_section(sb, "topB", BAR_2, [5, 6], [8, 7, 9, 8])

    cgnscheck(cgns)


def test_2d_tri(script, tmp):
    """Single tri block with one sideset."""
    coords = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], float)
    conn = np.array([[1, 2, 3], [2, 4, 3]])
    blocks = [("blkT-TRI", conn)]
    sidesets = [("right", [2], [2])]            # side 2 = edge (n2,n3)

    exo, cgns = f"{tmp}/mesh2dt.exo", f"{tmp}/mesh2dt.cgns"
    write_exodus(exo, coords, blocks, sidesets)
    convert(script, exo, cgns)

    with h5py.File(cgns) as h5:
        cell_dim, phys_dim, zones = read_base(h5)
        assert (cell_dim, phys_dim) == (2, 2), "base must be 2D"

        z = zones["blkT-TRI"]
        assert z[" data"][()].ravel().tolist() == [4, 2, 0]
        assert sorted(read_coords(z)) == ["CoordinateX", "CoordinateY"]

        secs = read_sections(z)
        check_section(secs, "blkT-TRI_cells", TRI_3, [1, 2], conn.ravel().tolist())
        check_section(secs, "right", BAR_2, [3, 3], [4, 3])

    cgnscheck(cgns)


def test_3d(script, tmp):
    """Single hex block with one sideset."""
    coords = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
                       [2, 0, 0], [2, 1, 0], [2, 0, 1], [2, 1, 1]], float)
    conn = np.array([[1, 2, 3, 4, 5, 6, 7, 8],
                     [2, 9, 10, 3, 6, 11, 12, 7]])
    blocks = [("blkC-HEX", conn)]
    sidesets = [("outlet", [2], [2])]                # face 2 = +x face

    exo, cgns = f"{tmp}/mesh3d.exo", f"{tmp}/mesh3d.cgns"
    write_exodus(exo, coords, blocks, sidesets)
    convert(script, exo, cgns)

    with h5py.File(cgns) as h5:
        cell_dim, phys_dim, zones = read_base(h5)
        assert (cell_dim, phys_dim) == (3, 3), "base must be 3D"

        z = zones["blkC-HEX"]
        assert z[" data"][()].ravel().tolist() == [12, 2, 0]
        assert sorted(read_coords(z)) == ["CoordinateX", "CoordinateY",
                                          "CoordinateZ"]

        secs = read_sections(z)
        check_section(secs, "blkC-HEX_cells", HEXA_8, [1, 2], conn.ravel().tolist())
        check_section(secs, "outlet", QUAD_4, [3, 3], [9, 10, 12, 11])

    cgnscheck(cgns)


def test_real_mesh(script, tmp):
    """Repository 2D overset mesh (skipped when not present)."""
    repo = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    exo = os.path.join(repo, "work", "overset_mesh", "overset_x12_large.exo")
    if not os.path.isfile(exo):
        print("  overset_x12_large.exo not found, skipped")
        return

    cgns = f"{tmp}/overset.cgns"
    convert(script, exo, cgns)

    nedge = 0
    with h5py.File(cgns) as h5:
        cell_dim, phys_dim, zones = read_base(h5)
        assert (cell_dim, phys_dim) == (2, 2), "base must be 2D"

        with netCDF4.Dataset(exo) as d:
            gx = np.asarray(d.variables["coordx"][:])
            gy = np.asarray(d.variables["coordy"][:])
            names = [s.strip() for s in
                     netCDF4.chartostring(d.variables["eb_names"][:])]
            blocks = {n: np.asarray(d.variables[f"connect{i + 1}"][:]) - 1
                      for i, n in enumerate(names)}

        for name, zone in zones.items():
            quad = blocks[name]

            # zone size and volume section
            assert zone[" data"][()].ravel().tolist() == \
                [np.unique(quad).size, quad.shape[0], 0]
            secs = read_sections(zone)
            vol = [s for s in secs.values() if s[0] == QUAD_4]
            assert len(vol) == 1, f"zone '{name}': one QUAD_4 section"

            # local renumbering: coordinates follow sorted global node ids
            local = np.unique(quad)
            c = read_coords(zone)
            assert c["CoordinateX"] == gx[local].tolist()
            assert c["CoordinateY"] == gy[local].tolist()

            # every BAR_2 boundary edge is an edge of some quad
            # (section connectivity holds 1-based local ids)
            edges = set()
            for q in quad:
                for a, b in zip(q, np.roll(q, -1)):
                    edges.add((min(a, b), max(a, b)))
            for sec in secs.values():
                if sec[0] != BAR_2:
                    continue
                nedge += len(sec[2]) // 2
                for a, b in zip(sec[2][::2], sec[2][1::2]):
                    ga, gb = local[a - 1], local[b - 1]
                    assert (min(ga, gb), max(ga, gb)) in edges, \
                        f"boundary edge ({a},{b}) is not an element edge"

        # sideset sizes of overset_x12_large.exo (fringe/wall/int/out)
        assert nedge == 9852 + 9852 + 1600 + 1600, "all sidesets transferred"

    cgnscheck(cgns)


def main():
    script = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "exodus_to_cgns.py")

    tests = [("2D quad mesh", test_2d),
             ("2D tri mesh", test_2d_tri),
             ("3D hex mesh", test_3d),
             ("real overset mesh", test_real_mesh)]
    failed = 0
    for name, fn in tests:
        print(f"[test] {name}")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                fn(script, tmp)
            print("[ok]")
        except AssertionError as e:
            failed += 1
            print(f"[FAILED] {e}")

    if failed:
        sys.exit(f"{failed} of {len(tests)} tests failed")
    print("all tests passed")


if __name__ == "__main__":
    main()
