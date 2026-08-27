#!/usr/bin/env python3
"""Add an OpenFOAM-style p_rgh nodal variable to an Exodus result file.

The OpenAccel pressure field is stored with a fixed reference-density
hydrostatic shift.  For a two-phase result, this script reconstructs the
local-density form

    p_rgh = pressure + (rho_mix - rho_ref) * (-g_y) * (y - y_ref)

for the case gravity convention ``g = (0, -g_y)``.  The output is a new
Exodus file so the original result is never overwritten unless --in-place is
explicitly requested.

Example (the bubble-column case):

    python3 add_prgh_to_exodus.py results.e \
        --output results_with_prgh.e

For results generated with the original gauge-pressure input, add the
absolute pressure reference with ``--pressure-offset 1e5``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def _decode_name(value: object) -> str:
    """Decode an Exodus fixed-width character-array entry."""
    if isinstance(value, np.ndarray):
        value = b"".join(np.asarray(value, dtype="S1").tolist())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore").rstrip("\x00 ")
    return str(value).rstrip("\x00 ")


def _nodal_variable_names(ds: Dataset) -> list[str]:
    names = ds.variables["name_nod_var"][:]
    return [_decode_name(row) for row in names]


def _copy_attributes(source: object, target: object) -> None:
    for name in source.ncattrs():
        # netCDF4 reserves _FillValue for createVariable(fill_value=...).
        # The Exodus arrays copied here do not rely on a custom fill value.
        if name == "_FillValue":
            continue
        target.setncattr(name, source.getncattr(name))


def _copy_dataset(source: Dataset, destination: Dataset) -> None:
    """Copy the complete Exodus dataset, enlarging num_nod_var by one."""
    for name, dimension in source.dimensions.items():
        size = len(dimension)
        if name == "num_nod_var":
            size += 1
        destination.createDimension(name, None if dimension.isunlimited() else size)

    for name, variable in source.variables.items():
        # The variable name table is rewritten below after adding p_rgh.
        if name == "name_nod_var":
            continue
        target = destination.createVariable(
            name,
            variable.datatype,
            variable.dimensions,
            zlib=getattr(variable, "compression", None) is not None,
            complevel=getattr(variable, "compression_opts", None)
            if getattr(variable, "compression", None) is not None
            else 0,
            shuffle=getattr(variable, "shuffle", True),
            fletcher32=getattr(variable, "fletcher32", False),
            chunksizes=getattr(variable, "chunking", lambda: None)()
            if isinstance(getattr(variable, "chunking", lambda: None)(), tuple)
            else None,
        )
        _copy_attributes(variable, target)
        target[:] = variable[:]

    source_names = _nodal_variable_names(source)
    name_table = destination.createVariable(
        "name_nod_var", source.variables["name_nod_var"].datatype, ("num_nod_var", "len_name")
    )
    _copy_attributes(source.variables["name_nod_var"], name_table)
    name_table[:] = np.full(
        (len(source_names) + 1, len(source.dimensions["len_name"])),
        b" ",
        dtype="S1",
    )
    name_table[: len(source_names), :] = source.variables["name_nod_var"][:]
    name_table[len(source_names), : len("p_rgh")] = np.array(
        list("p_rgh"), dtype="S1"
    )

    _copy_attributes(source, destination)


def _create_prgh_file(
    input_path: Path,
    output_path: Path,
    rho_air: float,
    rho_water: float,
    rho_ref: float,
    gravity_y: float,
    y_ref: float,
    pressure_offset: float,
) -> None:
    with Dataset(input_path, "r") as source:
        names = _nodal_variable_names(source)
        required = {"pressure", "volume_fraction.air", "volume_fraction.water"}
        missing = sorted(required.difference(names))
        if missing:
            raise ValueError(
                "Missing nodal variables: "
                + ", ".join(missing)
                + ". Available names: "
                + ", ".join(names)
            )

        if "p_rgh" in names:
            raise ValueError("The input already contains a p_rgh nodal variable")

        pressure_id = names.index("pressure") + 1
        alpha_air_id = names.index("volume_fraction.air") + 1
        alpha_water_id = names.index("volume_fraction.water") + 1
        n_steps = len(source.dimensions["time_step"])
        n_nodes = len(source.dimensions["num_nodes"])
        y = np.asarray(source.variables["coordy"][:], dtype=np.float64)
        if y.shape != (n_nodes,):
            raise ValueError("coordy does not match the Exodus num_nodes dimension")

        temporary = output_path.with_suffix(output_path.suffix + ".tmp")
        if temporary.exists():
            temporary.unlink()
        try:
            with Dataset(temporary, "w", format=source.file_format) as destination:
                _copy_dataset(source, destination)
                result = destination.createVariable(
                    "vals_nod_var8",
                    source.variables[f"vals_nod_var{pressure_id}"].datatype,
                    ("time_step", "num_nodes"),
                )
                _copy_attributes(source.variables[f"vals_nod_var{pressure_id}"], result)
                result.setncattr("long_name", "OpenFOAM-style p_rgh")

                pressure = source.variables[f"vals_nod_var{pressure_id}"]
                alpha_air = source.variables[f"vals_nod_var{alpha_air_id}"]
                alpha_water = source.variables[f"vals_nod_var{alpha_water_id}"]
                for step in range(n_steps):
                    rho_mix = (
                        rho_air * np.asarray(alpha_air[step, :], dtype=np.float64)
                        + rho_water * np.asarray(alpha_water[step, :], dtype=np.float64)
                    )
                    # OpenAccel stores p - rho_ref*(g dot (r-r_ref)).
                    # With g=(0, gravity_y), the local-density correction is:
                    # p_rgh = pressure + (rho_mix-rho_ref)*(-gravity_y)*(y-y_ref).
                    result[step, :] = (
                        np.asarray(pressure[step, :], dtype=np.float64)
                        + pressure_offset
                        + (rho_mix - rho_ref)
                        * (-gravity_y)
                        * (y - y_ref)
                    )
        except Exception:
            if temporary.exists():
                temporary.unlink()
            raise
        temporary.replace(output_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input Exodus result file")
    parser.add_argument(
        "--output",
        type=Path,
        help="Output Exodus file (default: <input stem>_with_prgh.e)",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Replace the input file after successful conversion",
    )
    parser.add_argument("--rho-air", type=float, default=1.2)
    parser.add_argument("--rho-water", type=float, default=1000.0)
    parser.add_argument("--rho-ref", type=float, default=1000.0)
    parser.add_argument("--gravity-y", type=float, default=-9.81)
    parser.add_argument("--y-ref", type=float, default=0.0)
    parser.add_argument(
        "--pressure-offset",
        type=float,
        default=0.0,
        help="Constant pressure offset, e.g. 1e5 for the original gauge input",
    )
    args = parser.parse_args()

    input_path = args.input.resolve()
    if not input_path.is_file():
        parser.error(f"input file does not exist: {input_path}")

    if args.in_place and args.output is not None:
        parser.error("use either --output or --in-place, not both")

    if args.in_place:
        output_path = input_path.with_suffix(input_path.suffix + ".new")
    elif args.output is None:
        output_path = input_path.with_name(input_path.stem + "_with_prgh" + input_path.suffix)
    else:
        output_path = args.output.resolve()

    if output_path == input_path:
        parser.error("output must differ from input unless --in-place is used")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _create_prgh_file(
        input_path,
        output_path,
        args.rho_air,
        args.rho_water,
        args.rho_ref,
        args.gravity_y,
        args.y_ref,
        args.pressure_offset,
    )

    if args.in_place:
        backup = input_path.with_suffix(input_path.suffix + ".bak")
        shutil.copy2(input_path, backup)
        output_path.replace(input_path)
        print(f"Added p_rgh to {input_path} (backup: {backup})")
    else:
        print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
