"""Generate the small 2D box mesh for the settlingColumn verification case.

Exodus II layout mirrors examples/cavity/mesh.e: one QUAD4 element block named
`fluid`, plus four side sets. `per1`/`per2` are the x-normal faces paired by the
translational-periodicity interface; `bottom`/`top` are the y-normal faces.
"""
import netCDF4 as nc
import numpy as np

NX, NY = 4, 100          # elements
LX, LY = 0.004, 0.05     # domain size [m]
PATH = "mesh.e"

nnx, nny = NX + 1, NY + 1
num_nodes, num_elem = nnx * nny, NX * NY

x = np.tile(np.linspace(0.0, LX, nnx), nny)
y = np.repeat(np.linspace(0.0, LY, nny), nnx)

# QUAD4 connectivity, counter-clockwise, 1-based
nid = lambda i, j: j * nnx + i + 1
connect = np.array(
    [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
     for j in range(NY) for i in range(NX)],
    dtype=np.int32,
)

eid = lambda i, j: j * NX + i + 1
# Exodus QUAD4 side numbering: 1=(n1,n2) bottom, 2=(n2,n3) right,
#                              3=(n3,n4) top,    4=(n4,n1) left
sidesets = {
    "left":   [(eid(0, j), 4) for j in range(NY)],        # x = 0
    "right":  [(eid(NX - 1, j), 2) for j in range(NY)],   # x = LX
    "bottom": [(eid(i, 0), 1) for i in range(NX)],        # y = 0
    "top":    [(eid(i, NY - 1), 3) for i in range(NX)],   # y = LY
}

d = nc.Dataset(PATH, "w", format="NETCDF3_64BIT_OFFSET")
d.api_version = np.float32(7.03)
d.version = np.float32(7.03)
d.floating_point_word_size = np.int32(8)
d.file_size = np.int32(1)
d.maximum_name_length = np.int32(32)
d.int64_status = np.int32(0)
d.title = "slipRelaxation Euler-Euler verification case"

d.createDimension("len_string", 33)
d.createDimension("len_line", 81)
d.createDimension("four", 4)
d.createDimension("len_name", 256)
d.createDimension("time_step", None)
d.createDimension("num_dim", 2)
d.createDimension("num_nodes", num_nodes)
d.createDimension("num_elem", num_elem)
d.createDimension("num_el_blk", 1)
d.createDimension("num_el_in_blk1", num_elem)
d.createDimension("num_nod_per_el1", 4)
d.createDimension("num_side_sets", len(sidesets))

d.createVariable("time_whole", "f8", ("time_step",))

v = d.createVariable("eb_status", "i4", ("num_el_blk",)); v[:] = [1]
v = d.createVariable("eb_prop1", "i4", ("num_el_blk",)); v.setncattr("name", "ID"); v[:] = [1]

d.createVariable("coordx", "f8", ("num_nodes",))[:] = x
d.createVariable("coordy", "f8", ("num_nodes",))[:] = y

def put_names(varname, dim, names):
    v = d.createVariable(varname, "S1", (dim, "len_string"))
    arr = np.zeros((len(names), 33), dtype="S1")
    for i, nm in enumerate(names):
        for k, ch in enumerate(nm):
            arr[i, k] = ch.encode()
    v[:] = arr

put_names("eb_names", "num_el_blk", ["fluid"])
put_names("coor_names", "num_dim", ["xcoor", "ycoor"])

v = d.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
v.elem_type = "QUAD4"
v[:] = connect

v = d.createVariable("ss_status", "i4", ("num_side_sets",)); v[:] = [1] * len(sidesets)
v = d.createVariable("ss_prop1", "i4", ("num_side_sets",)); v.setncattr("name", "ID")
v[:] = list(range(1, len(sidesets) + 1))
put_names("ss_names", "num_side_sets", list(sidesets))

for i, (name, faces) in enumerate(sidesets.items(), start=1):
    d.createDimension(f"num_side_ss{i}", len(faces))
    d.createVariable(f"elem_ss{i}", "i4", (f"num_side_ss{i}",))[:] = \
        np.array([f[0] for f in faces], dtype=np.int32)
    d.createVariable(f"side_ss{i}", "i4", (f"num_side_ss{i}",))[:] = \
        np.array([f[1] for f in faces], dtype=np.int32)

d.close()
print(f"wrote {PATH}: {num_nodes} nodes, {num_elem} elements, sidesets {list(sidesets)}")
