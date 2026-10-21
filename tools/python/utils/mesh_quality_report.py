#!/usr/bin/env python3
"""
Comprehensive 3-D mesh-quality inspection using VTK/ParaView.

The script is intended especially for Exodus II meshes, but it also reads
VTU/PVTU/VTM and legacy VTK unstructured-grid files. It uses the same VTK
quality implementation used by ParaView, so the scaled-Jacobian values should
match ParaView's Mesh Quality filter for supported linear cell types.

Recommended invocation when ParaView is installed:

    /path/to/paraview/bin/pvpython mesh_quality_report.py mesh.e

Or with a normal Python installation containing VTK:

    python3 -m pip install vtk
    python3 mesh_quality_report.py mesh.e

Outputs, under <mesh_stem>_quality_report/ by default:

    summary.txt
    all_volume_cells.csv
    problem_cells.csv
    annotated_volume_mesh.vtu

The annotated VTU contains cell arrays that can be coloured or thresholded in
ParaView, including ScaledJacobian, Jacobian, SignedVolume, MinEdgeLength,
EdgeRatio, MinFaceArea, MaxFaceWarpageDeg, and QualityClass.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import vtk  # type: ignore
except ImportError as exc:
    raise SystemExit(
        "VTK is not available in this Python environment.\n"
        "Run this script with ParaView's pvpython, for example:\n"
        "  /path/to/paraview/bin/pvpython mesh_quality_report.py mesh.e\n"
        "or install VTK with:\n"
        "  python3 -m pip install vtk"
    ) from exc


# Classification codes written to the annotated VTU.
CLASS_CODES = {
    "GOOD": 0,
    "ACCEPTABLE": 1,
    "POOR": 2,
    "BAD": 3,
    "DEGENERATE": 4,
    "INVERTED": 5,
    "UNSUPPORTED": 6,
}

SUPPORTED_3D_TYPES = {
    vtk.VTK_TETRA: "TETRA4",
    vtk.VTK_PYRAMID: "PYRAMID5",
    vtk.VTK_WEDGE: "WEDGE6",
    vtk.VTK_HEXAHEDRON: "HEX8",
}

QUALITY_FUNCTIONS = {
    vtk.VTK_TETRA: {
        "scaled_jacobian": "TetScaledJacobian",
        "jacobian": "TetJacobian",
        "volume": "TetVolume",
        "shape": "TetShape",
        "condition": "TetCondition",
        "max_aspect_frobenius": "TetAspectFrobenius",
        "equiangle_skew": "TetEquiangleSkew",
    },
    vtk.VTK_PYRAMID: {
        "scaled_jacobian": "PyramidScaledJacobian",
        "jacobian": "PyramidJacobian",
        "volume": "PyramidVolume",
        "shape": "PyramidShape",
        "condition": None,
        "max_aspect_frobenius": None,
        "equiangle_skew": "PyramidEquiangleSkew",
    },
    vtk.VTK_WEDGE: {
        "scaled_jacobian": "WedgeScaledJacobian",
        "jacobian": "WedgeJacobian",
        "volume": "WedgeVolume",
        "shape": "WedgeShape",
        "condition": "WedgeCondition",
        "max_aspect_frobenius": "WedgeMaxAspectFrobenius",
        "equiangle_skew": "WedgeEquiangleSkew",
    },
    vtk.VTK_HEXAHEDRON: {
        "scaled_jacobian": "HexScaledJacobian",
        "jacobian": "HexJacobian",
        "volume": "HexVolume",
        "shape": "HexShape",
        "condition": "HexCondition",
        "max_aspect_frobenius": "HexMaxAspectFrobenius",
        "equiangle_skew": "HexEquiangleSkew",
    },
}

CSV_FIELDS = [
    "block_flat_index",
    "block_name",
    "local_cell_id",
    "original_element_id",
    "cell_type",
    "vtk_cell_type",
    "classification",
    "scaled_jacobian",
    "jacobian",
    "signed_volume",
    "absolute_volume",
    "shape",
    "condition",
    "max_aspect_frobenius",
    "equiangle_skew",
    "min_edge_length",
    "max_edge_length",
    "edge_ratio",
    "min_face_area",
    "max_face_warpage_deg",
    "pyramid_relative_height",
    "pyramid_base_warpage_deg",
    "repeated_node_id_flag",
    "repeated_node_ids",
    "coincident_point_flag",
    "coincident_local_pairs",
    "edge_degenerate_flag",
    "face_degenerate_flag",
    "volume_degenerate_flag",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "point_ids",
    "global_point_ids",
    "point_coordinates",
]


def finite(value: float) -> bool:
    return math.isfinite(value)


def vec_sub(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_add(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vec_scale(a: Sequence[float], s: float) -> Tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a: Sequence[float]) -> float:
    return math.sqrt(dot(a, a))


def distance(a: Sequence[float], b: Sequence[float]) -> float:
    return norm(vec_sub(a, b))


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not finite(denominator):
        return math.inf
    return numerator / denominator


def percentile(sorted_values: Sequence[float], fraction: float) -> float:
    if not sorted_values:
        return math.nan
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return sorted_values[lo]
    weight = position - lo
    return sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight


def format_number(value: Any) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return f"{value:.12g}"
    return str(value)


def safe_vtk_quality(function_name: Optional[str], cell: Any) -> float:
    if not function_name:
        return math.nan
    function = getattr(vtk.vtkMeshQuality, function_name, None)
    if function is None:
        return math.nan
    try:
        return float(function(cell))
    except Exception:
        return math.nan


def read_mesh(path: Path) -> Any:
    suffix = path.suffix.lower()
    name_lower = path.name.lower()

    if suffix in {".e", ".exo", ".ex2", ".g", ".gen"} or name_lower.endswith(".e-s"):
        reader = vtk.vtkExodusIIReader()
        reader.SetFileName(str(path))
        reader.UpdateInformation()

        # Request all element blocks and, where available, original/global IDs.
        try:
            reader.SetAllArrayStatus(vtk.vtkExodusIIReader.ELEM_BLOCK, 1)
        except Exception:
            pass
        for method_name in (
            "GenerateObjectIdCellArrayOn",
            "GenerateGlobalElementIdArrayOn",
            "GenerateGlobalNodeIdArrayOn",
        ):
            method = getattr(reader, method_name, None)
            if method is not None:
                try:
                    method()
                except Exception:
                    pass

    elif suffix == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(path))
    elif suffix == ".pvtu":
        reader = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(str(path))
    elif suffix == ".vtm":
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader.SetFileName(str(path))
    elif suffix == ".vtk":
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(str(path))
    else:
        raise ValueError(
            f"Unsupported file extension '{suffix}'. Supported: Exodus (.e/.exo/.ex2/.g), "
            ".vtu, .pvtu, .vtm, and legacy .vtk."
        )

    reader.Update()
    output = reader.GetOutputDataObject(0)
    if output is None:
        raise RuntimeError(f"VTK could not read '{path}'.")
    return output


def iter_leaf_datasets(data_object: Any) -> Iterable[Tuple[int, str, Any]]:
    """Yield (flat index, block name, vtkDataSet) for each composite leaf."""
    if data_object.IsA("vtkCompositeDataSet"):
        iterator = data_object.NewIterator()
        iterator.VisitOnlyLeavesOn()
        iterator.SkipEmptyNodesOn()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            leaf = iterator.GetCurrentDataObject()
            flat_index = int(iterator.GetCurrentFlatIndex())
            block_name = f"block_{flat_index}"
            try:
                metadata = iterator.GetCurrentMetaData()
                if metadata is not None and metadata.Has(vtk.vtkCompositeDataSet.NAME()):
                    block_name = str(metadata.Get(vtk.vtkCompositeDataSet.NAME()))
            except Exception:
                pass
            if leaf is not None and leaf.IsA("vtkDataSet"):
                yield flat_index, block_name, leaf
            iterator.GoToNextItem()
    elif data_object.IsA("vtkDataSet"):
        yield 0, "mesh", data_object
    else:
        raise TypeError(f"Unsupported VTK output type: {data_object.GetClassName()}")


def find_id_array(data_attributes: Any, candidate_names: Sequence[str]) -> Any:
    if data_attributes is None:
        return None
    try:
        global_ids = data_attributes.GetGlobalIds()
        if global_ids is not None:
            return global_ids
    except Exception:
        pass
    for name in candidate_names:
        array = data_attributes.GetArray(name)
        if array is not None:
            return array
    return None


def get_integral_value(array: Any, index: int) -> Optional[int]:
    if array is None:
        return None
    try:
        return int(round(float(array.GetTuple1(index))))
    except Exception:
        return None


def get_cell_points(grid: Any, cell: Any) -> Tuple[List[int], List[Tuple[float, float, float]]]:
    point_ids: List[int] = []
    points: List[Tuple[float, float, float]] = []
    for local_id in range(cell.GetNumberOfPoints()):
        point_id = int(cell.GetPointId(local_id))
        point_ids.append(point_id)
        point = grid.GetPoint(point_id)
        points.append((float(point[0]), float(point[1]), float(point[2])))
    return point_ids, points


def repeated_ids(point_ids: Sequence[int]) -> List[int]:
    counts = Counter(point_ids)
    return sorted(point_id for point_id, count in counts.items() if count > 1)


def coincident_pairs(
    points: Sequence[Sequence[float]], tolerance: float
) -> List[Tuple[int, int]]:
    pairs: List[Tuple[int, int]] = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if distance(points[i], points[j]) <= tolerance:
                pairs.append((i, j))
    return pairs


def cell_edge_lengths(grid: Any, cell: Any) -> List[float]:
    lengths: List[float] = []
    for edge_id in range(cell.GetNumberOfEdges()):
        edge = cell.GetEdge(edge_id)
        if edge is None or edge.GetNumberOfPoints() < 2:
            continue
        p0 = grid.GetPoint(edge.GetPointId(0))
        p1 = grid.GetPoint(edge.GetPointId(1))
        lengths.append(distance(p0, p1))
    return lengths


def polygon_area(points: Sequence[Sequence[float]]) -> float:
    """Triangulated area using a fan from vertex 0."""
    if len(points) < 3:
        return 0.0
    area = 0.0
    p0 = points[0]
    for i in range(1, len(points) - 1):
        a = vec_sub(points[i], p0)
        b = vec_sub(points[i + 1], p0)
        area += 0.5 * norm(cross(a, b))
    return area


def quad_warpage_angle_deg(points: Sequence[Sequence[float]]) -> float:
    """Angle between the normals of the two triangles split by diagonal 0-2."""
    if len(points) != 4:
        return 0.0
    n1 = cross(vec_sub(points[1], points[0]), vec_sub(points[2], points[0]))
    n2 = cross(vec_sub(points[2], points[0]), vec_sub(points[3], points[0]))
    mag1 = norm(n1)
    mag2 = norm(n2)
    if mag1 <= 0.0 or mag2 <= 0.0:
        return 180.0
    cosine = max(-1.0, min(1.0, dot(n1, n2) / (mag1 * mag2)))
    return math.degrees(math.acos(cosine))


def face_metrics(grid: Any, cell: Any) -> Tuple[float, float]:
    face_areas: List[float] = []
    max_warpage = 0.0
    for face_id in range(cell.GetNumberOfFaces()):
        face = cell.GetFace(face_id)
        if face is None:
            continue
        face_points = [grid.GetPoint(face.GetPointId(i)) for i in range(face.GetNumberOfPoints())]
        face_areas.append(polygon_area(face_points))
        if len(face_points) == 4:
            max_warpage = max(max_warpage, quad_warpage_angle_deg(face_points))
    min_face_area = min(face_areas) if face_areas else math.nan
    return min_face_area, max_warpage


def pyramid_diagnostics(points: Sequence[Sequence[float]]) -> Tuple[float, float]:
    """Return apex-height/base-edge and base warpage for a VTK PYRAMID5."""
    if len(points) != 5:
        return math.nan, math.nan
    base = points[:4]
    apex = points[4]
    base_center = (0.0, 0.0, 0.0)
    for point in base:
        base_center = vec_add(base_center, point)
    base_center = vec_scale(base_center, 0.25)

    # Newell normal is robust for a quad, including mildly warped bases.
    nx = ny = nz = 0.0
    for i in range(4):
        current = base[i]
        nxt = base[(i + 1) % 4]
        nx += (current[1] - nxt[1]) * (current[2] + nxt[2])
        ny += (current[2] - nxt[2]) * (current[0] + nxt[0])
        nz += (current[0] - nxt[0]) * (current[1] + nxt[1])
    normal = (nx, ny, nz)
    normal_length = norm(normal)
    if normal_length <= 0.0:
        relative_height = 0.0
    else:
        unit_normal = vec_scale(normal, 1.0 / normal_length)
        height = abs(dot(vec_sub(apex, base_center), unit_normal))
        base_edges = [distance(base[i], base[(i + 1) % 4]) for i in range(4)]
        mean_base_edge = sum(base_edges) / 4.0 if base_edges else 0.0
        relative_height = safe_ratio(height, mean_base_edge)

    return relative_height, quad_warpage_angle_deg(base)


def cell_centroid(points: Sequence[Sequence[float]]) -> Tuple[float, float, float]:
    if not points:
        return (math.nan, math.nan, math.nan)
    total = (0.0, 0.0, 0.0)
    for point in points:
        total = vec_add(total, point)
    return vec_scale(total, 1.0 / len(points))


def classify_cell(
    metrics: Dict[str, Any],
    args: argparse.Namespace,
) -> str:
    quality = metrics["scaled_jacobian"]
    jacobian = metrics["jacobian"]
    signed_volume = metrics["signed_volume"]
    local_volume_tolerance = metrics["local_volume_tolerance"]

    if not finite(quality):
        return "UNSUPPORTED"

    # Negative scaled/raw Jacobian is the strongest inversion indicator.
    if quality < 0.0 or (finite(jacobian) and jacobian < 0.0):
        return "INVERTED"
    if finite(signed_volume) and signed_volume < -local_volume_tolerance:
        return "INVERTED"

    if (
        metrics["repeated_node_id_flag"]
        or metrics["coincident_point_flag"]
        or metrics["edge_degenerate_flag"]
        or metrics["face_degenerate_flag"]
        or metrics["volume_degenerate_flag"]
        or quality <= args.degenerate_jacobian
    ):
        return "DEGENERATE"
    if quality < args.bad_jacobian:
        return "BAD"
    if quality < args.poor_jacobian:
        return "POOR"
    if quality < args.good_jacobian:
        return "ACCEPTABLE"
    return "GOOD"


def make_array(name: str, integer: bool = False, int64: bool = False) -> Any:
    if int64:
        array = vtk.vtkTypeInt64Array()
    elif integer:
        array = vtk.vtkIntArray()
    else:
        array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(1)
    return array


def make_annotated_grid(points: Any) -> Tuple[Any, Dict[str, Any]]:
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    arrays = {
        "ScaledJacobian": make_array("ScaledJacobian"),
        "Jacobian": make_array("Jacobian"),
        "SignedVolume": make_array("SignedVolume"),
        "AbsoluteVolume": make_array("AbsoluteVolume"),
        "Shape": make_array("Shape"),
        "Condition": make_array("Condition"),
        "MaxAspectFrobenius": make_array("MaxAspectFrobenius"),
        "EquiangleSkew": make_array("EquiangleSkew"),
        "MinEdgeLength": make_array("MinEdgeLength"),
        "MaxEdgeLength": make_array("MaxEdgeLength"),
        "EdgeRatio": make_array("EdgeRatio"),
        "MinFaceArea": make_array("MinFaceArea"),
        "MaxFaceWarpageDeg": make_array("MaxFaceWarpageDeg"),
        "PyramidRelativeHeight": make_array("PyramidRelativeHeight"),
        "PyramidBaseWarpageDeg": make_array("PyramidBaseWarpageDeg"),
        "RepeatedNodeFlag": make_array("RepeatedNodeFlag", integer=True),
        "CoincidentPointFlag": make_array("CoincidentPointFlag", integer=True),
        "EdgeDegenerateFlag": make_array("EdgeDegenerateFlag", integer=True),
        "FaceDegenerateFlag": make_array("FaceDegenerateFlag", integer=True),
        "VolumeDegenerateFlag": make_array("VolumeDegenerateFlag", integer=True),
        "QualityClass": make_array("QualityClass", integer=True),
        "OriginalBlockFlatIndex": make_array("OriginalBlockFlatIndex", integer=True),
        "OriginalCellId": make_array("OriginalCellId", int64=True),
        "OriginalElementId": make_array("OriginalElementId", int64=True),
    }
    for array in arrays.values():
        grid.GetCellData().AddArray(array)
    return grid, arrays


def append_annotated_values(arrays: Dict[str, Any], row: Dict[str, Any]) -> None:
    numeric_map = {
        "ScaledJacobian": "scaled_jacobian",
        "Jacobian": "jacobian",
        "SignedVolume": "signed_volume",
        "AbsoluteVolume": "absolute_volume",
        "Shape": "shape",
        "Condition": "condition",
        "MaxAspectFrobenius": "max_aspect_frobenius",
        "EquiangleSkew": "equiangle_skew",
        "MinEdgeLength": "min_edge_length",
        "MaxEdgeLength": "max_edge_length",
        "EdgeRatio": "edge_ratio",
        "MinFaceArea": "min_face_area",
        "MaxFaceWarpageDeg": "max_face_warpage_deg",
        "PyramidRelativeHeight": "pyramid_relative_height",
        "PyramidBaseWarpageDeg": "pyramid_base_warpage_deg",
    }
    for array_name, row_name in numeric_map.items():
        arrays[array_name].InsertNextValue(float(row[row_name]))

    arrays["RepeatedNodeFlag"].InsertNextValue(int(row["repeated_node_id_flag"]))
    arrays["CoincidentPointFlag"].InsertNextValue(int(row["coincident_point_flag"]))
    arrays["EdgeDegenerateFlag"].InsertNextValue(int(row["edge_degenerate_flag"]))
    arrays["FaceDegenerateFlag"].InsertNextValue(int(row["face_degenerate_flag"]))
    arrays["VolumeDegenerateFlag"].InsertNextValue(int(row["volume_degenerate_flag"]))
    arrays["QualityClass"].InsertNextValue(CLASS_CODES[row["classification"]])
    arrays["OriginalBlockFlatIndex"].InsertNextValue(int(row["block_flat_index"]))
    arrays["OriginalCellId"].InsertNextValue(int(row["local_cell_id"]))
    original_element_id = row["original_element_id"]
    arrays["OriginalElementId"].InsertNextValue(
        int(original_element_id) if original_element_id not in (None, "") else -1
    )


def analyze_cell(
    grid: Any,
    cell: Any,
    block_flat_index: int,
    block_name: str,
    local_cell_id: int,
    global_point_ids_array: Any,
    global_cell_ids_array: Any,
    args: argparse.Namespace,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    cell_type_id = int(cell.GetCellType())
    cell_type_name = SUPPORTED_3D_TYPES.get(
        cell_type_id,
        vtk.vtkCellTypes.GetClassNameFromTypeId(cell_type_id) or f"VTK_TYPE_{cell_type_id}",
    )
    point_ids, points = get_cell_points(grid, cell)
    global_point_ids = [get_integral_value(global_point_ids_array, pid) for pid in point_ids]
    original_element_id = get_integral_value(global_cell_ids_array, local_cell_id)

    edge_lengths = cell_edge_lengths(grid, cell)
    min_edge = min(edge_lengths) if edge_lengths else math.nan
    max_edge = max(edge_lengths) if edge_lengths else math.nan
    edge_ratio = safe_ratio(max_edge, min_edge) if edge_lengths else math.nan
    characteristic_length = max_edge if finite(max_edge) and max_edge > 0.0 else 0.0

    point_tolerance = max(args.absolute_tolerance, args.coincident_relative_tolerance * characteristic_length)
    edge_tolerance = max(args.absolute_tolerance, args.edge_relative_tolerance * characteristic_length)
    area_tolerance = args.area_relative_tolerance * characteristic_length**2
    volume_tolerance = args.volume_relative_tolerance * characteristic_length**3

    duplicate_ids = repeated_ids(point_ids)
    duplicate_coordinate_pairs = coincident_pairs(points, point_tolerance)
    min_face_area, max_face_warpage = face_metrics(grid, cell)

    quality_functions = QUALITY_FUNCTIONS.get(cell_type_id, {})
    quality_metrics: Dict[str, float] = {}
    for metric_name in (
        "scaled_jacobian",
        "jacobian",
        "volume",
        "shape",
        "condition",
        "max_aspect_frobenius",
        "equiangle_skew",
    ):
        quality_metrics[metric_name] = safe_vtk_quality(quality_functions.get(metric_name), cell)

    signed_volume = quality_metrics["volume"]
    absolute_volume = abs(signed_volume) if finite(signed_volume) else math.nan
    pyramid_relative_height = math.nan
    pyramid_base_warpage = math.nan
    if cell_type_id == vtk.VTK_PYRAMID:
        pyramid_relative_height, pyramid_base_warpage = pyramid_diagnostics(points)

    geometry_flags = {
        "repeated_node_id_flag": bool(duplicate_ids),
        "coincident_point_flag": bool(duplicate_coordinate_pairs),
        "edge_degenerate_flag": finite(min_edge) and min_edge <= edge_tolerance,
        "face_degenerate_flag": finite(min_face_area) and min_face_area <= area_tolerance,
        "volume_degenerate_flag": finite(absolute_volume) and absolute_volume <= volume_tolerance,
        "local_volume_tolerance": volume_tolerance,
    }

    classification_inputs: Dict[str, Any] = {
        **quality_metrics,
        "signed_volume": signed_volume,
        **geometry_flags,
    }
    classification = classify_cell(classification_inputs, args)
    centroid = cell_centroid(points)

    row: Dict[str, Any] = {
        "block_flat_index": block_flat_index,
        "block_name": block_name,
        "local_cell_id": local_cell_id,
        "original_element_id": original_element_id if original_element_id is not None else "",
        "cell_type": cell_type_name,
        "vtk_cell_type": cell_type_id,
        "classification": classification,
        "scaled_jacobian": quality_metrics["scaled_jacobian"],
        "jacobian": quality_metrics["jacobian"],
        "signed_volume": signed_volume,
        "absolute_volume": absolute_volume,
        "shape": quality_metrics["shape"],
        "condition": quality_metrics["condition"],
        "max_aspect_frobenius": quality_metrics["max_aspect_frobenius"],
        "equiangle_skew": quality_metrics["equiangle_skew"],
        "min_edge_length": min_edge,
        "max_edge_length": max_edge,
        "edge_ratio": edge_ratio,
        "min_face_area": min_face_area,
        "max_face_warpage_deg": max_face_warpage,
        "pyramid_relative_height": pyramid_relative_height,
        "pyramid_base_warpage_deg": pyramid_base_warpage,
        "repeated_node_id_flag": int(bool(duplicate_ids)),
        "repeated_node_ids": ";".join(str(value) for value in duplicate_ids),
        "coincident_point_flag": int(bool(duplicate_coordinate_pairs)),
        "coincident_local_pairs": ";".join(f"{a}-{b}" for a, b in duplicate_coordinate_pairs),
        "edge_degenerate_flag": int(geometry_flags["edge_degenerate_flag"]),
        "face_degenerate_flag": int(geometry_flags["face_degenerate_flag"]),
        "volume_degenerate_flag": int(geometry_flags["volume_degenerate_flag"]),
        "centroid_x": centroid[0],
        "centroid_y": centroid[1],
        "centroid_z": centroid[2],
        "point_ids": ";".join(str(value) for value in point_ids),
        "global_point_ids": ";".join(
            "" if value is None else str(value) for value in global_point_ids
        ),
        "point_coordinates": "|".join(
            f"{point[0]:.12g},{point[1]:.12g},{point[2]:.12g}" for point in points
        ),
    }

    internal = {
        "local_volume_tolerance": volume_tolerance,
        "point_ids_object": cell.GetPointIds(),
    }
    return row, internal


def write_summary(
    summary_path: Path,
    input_path: Path,
    args: argparse.Namespace,
    total_volume_cells: int,
    ignored_non_volume_cells: int,
    type_counts: Counter,
    class_counts: Counter,
    type_class_counts: Dict[str, Counter],
    qualities_by_type: Dict[str, List[float]],
    issue_counts: Counter,
    worst_rows: List[Dict[str, Any]],
) -> None:
    lines: List[str] = []
    lines.append("MESH QUALITY REPORT")
    lines.append("=" * 80)
    lines.append(f"Input mesh: {input_path.resolve()}")
    lines.append(f"VTK version: {vtk.vtkVersion.GetVTKVersion()}")
    lines.append("")
    lines.append("Classification thresholds")
    lines.append("-" * 80)
    lines.append(f"INVERTED: scaled/raw Jacobian < 0")
    lines.append(
        f"DEGENERATE: scaled Jacobian <= {args.degenerate_jacobian:g}, repeated/coincident "
        "vertices, or near-zero edge/face/volume"
    )
    lines.append(f"BAD:        {args.degenerate_jacobian:g} < scaled Jacobian < {args.bad_jacobian:g}")
    lines.append(f"POOR:       {args.bad_jacobian:g} <= scaled Jacobian < {args.poor_jacobian:g}")
    lines.append(f"ACCEPTABLE: {args.poor_jacobian:g} <= scaled Jacobian < {args.good_jacobian:g}")
    lines.append(f"GOOD:       scaled Jacobian >= {args.good_jacobian:g}")
    lines.append("")
    lines.append("Mesh totals")
    lines.append("-" * 80)
    lines.append(f"3-D volume cells inspected: {total_volume_cells}")
    lines.append(f"Non-volume cells ignored:   {ignored_non_volume_cells}")
    for cell_type, count in sorted(type_counts.items()):
        lines.append(f"  {cell_type:<24} {count:>12}")
    lines.append("")
    lines.append("Classification totals")
    lines.append("-" * 80)
    for class_name in CLASS_CODES:
        lines.append(f"  {class_name:<24} {class_counts.get(class_name, 0):>12}")
    lines.append("")
    lines.append("Detected geometric issues")
    lines.append("-" * 80)
    for issue_name in (
        "repeated_node_ids",
        "coincident_points",
        "near_zero_edges",
        "near_zero_faces",
        "near_zero_volumes",
    ):
        lines.append(f"  {issue_name:<24} {issue_counts.get(issue_name, 0):>12}")

    lines.append("")
    lines.append("Scaled-Jacobian statistics by topology")
    lines.append("-" * 80)
    lines.append(
        f"{'Type':<16} {'Count':>10} {'Min':>12} {'P1':>12} {'P5':>12} "
        f"{'Median':>12} {'P95':>12} {'Max':>12}"
    )
    for cell_type in sorted(type_counts):
        values = sorted(value for value in qualities_by_type.get(cell_type, []) if finite(value))
        if not values:
            lines.append(f"{cell_type:<16} {0:>10} {'n/a':>12}")
            continue
        lines.append(
            f"{cell_type:<16} {len(values):>10d} "
            f"{format_number(values[0]):>12} "
            f"{format_number(percentile(values, 0.01)):>12} "
            f"{format_number(percentile(values, 0.05)):>12} "
            f"{format_number(percentile(values, 0.50)):>12} "
            f"{format_number(percentile(values, 0.95)):>12} "
            f"{format_number(values[-1]):>12}"
        )
        classes = type_class_counts.get(cell_type, Counter())
        class_text = ", ".join(
            f"{name}={classes.get(name, 0)}" for name in CLASS_CODES if classes.get(name, 0)
        )
        lines.append(f"  classes: {class_text or 'none'}")

    lines.append("")
    lines.append(f"Worst {len(worst_rows)} cells by scaled Jacobian")
    lines.append("-" * 80)
    lines.append(
        f"{'Block':>7} {'Cell':>10} {'Element ID':>12} {'Type':<12} "
        f"{'Class':<12} {'ScaledJ':>12} {'Volume':>12} {'Min edge':>12}"
    )
    for row in worst_rows:
        lines.append(
            f"{row['block_flat_index']:>7} {row['local_cell_id']:>10} "
            f"{str(row['original_element_id']):>12} {row['cell_type']:<12} "
            f"{row['classification']:<12} {format_number(row['scaled_jacobian']):>12} "
            f"{format_number(row['signed_volume']):>12} "
            f"{format_number(row['min_edge_length']):>12}"
        )

    lines.append("")
    lines.append("Interpretation notes")
    lines.append("-" * 80)
    lines.append("* Scaled Jacobian is dimensionless; negative means inverted and zero means singular.")
    lines.append("* Near-zero edge/face/volume tests are relative to each cell's maximum edge length.")
    lines.append("* PYRAMID5 relative height is apex distance from the base plane divided by mean base-edge length.")
    lines.append("* MaxFaceWarpageDeg is 0 for planar quad faces and rises as the face twists.")
    lines.append("* Quality limits are screening thresholds, not universal solver-acceptance guarantees.")

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect 3-D mesh quality using VTK/ParaView quality metrics."
    )
    parser.add_argument("mesh", type=Path, help="Input mesh (.e/.exo/.vtu/.pvtu/.vtm/.vtk)")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; default: <mesh_stem>_quality_report",
    )
    parser.add_argument(
        "--degenerate-jacobian",
        type=float,
        default=1.0e-3,
        help="Scaled-Jacobian threshold at or below which a cell is considered degenerate (default: 1e-3)",
    )
    parser.add_argument(
        "--bad-jacobian",
        type=float,
        default=0.1,
        help="Cells below this threshold are BAD after degeneracy checks (default: 0.1)",
    )
    parser.add_argument(
        "--poor-jacobian",
        type=float,
        default=0.3,
        help="Cells below this threshold are POOR after BAD checks (default: 0.3)",
    )
    parser.add_argument(
        "--good-jacobian",
        type=float,
        default=0.5,
        help="Cells at or above this threshold are GOOD (default: 0.5)",
    )
    parser.add_argument(
        "--coincident-relative-tolerance",
        type=float,
        default=1.0e-10,
        help="Coincident-point tolerance relative to cell maximum edge (default: 1e-10)",
    )
    parser.add_argument(
        "--edge-relative-tolerance",
        type=float,
        default=1.0e-10,
        help="Near-zero edge tolerance relative to cell maximum edge (default: 1e-10)",
    )
    parser.add_argument(
        "--area-relative-tolerance",
        type=float,
        default=1.0e-12,
        help="Near-zero face-area tolerance relative to max_edge^2 (default: 1e-12)",
    )
    parser.add_argument(
        "--volume-relative-tolerance",
        type=float,
        default=1.0e-12,
        help="Near-zero volume tolerance relative to max_edge^3 (default: 1e-12)",
    )
    parser.add_argument(
        "--absolute-tolerance",
        type=float,
        default=0.0,
        help="Absolute lower bound for point/edge coincidence tests (default: 0)",
    )
    parser.add_argument(
        "--worst-count",
        type=int,
        default=50,
        help="Number of worst cells listed in summary.txt (default: 50)",
    )
    parser.add_argument(
        "--skip-all-cells-csv",
        action="store_true",
        help="Do not write all_volume_cells.csv; problem_cells.csv is still written",
    )
    parser.add_argument(
        "--no-annotated-mesh",
        action="store_true",
        help="Do not write annotated_volume_mesh.vtu",
    )
    return parser.parse_args()


def validate_thresholds(args: argparse.Namespace) -> None:
    if not (
        args.degenerate_jacobian <= args.bad_jacobian <= args.poor_jacobian <= args.good_jacobian
    ):
        raise ValueError(
            "Thresholds must satisfy: degenerate <= bad <= poor <= good."
        )
    for name in (
        "coincident_relative_tolerance",
        "edge_relative_tolerance",
        "area_relative_tolerance",
        "volume_relative_tolerance",
        "absolute_tolerance",
    ):
        if getattr(args, name) < 0.0:
            raise ValueError(f"--{name.replace('_', '-')} must be non-negative.")


def main() -> int:
    args = parse_arguments()
    validate_thresholds(args)
    input_path = args.mesh.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Mesh file not found: {input_path}")

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = input_path.parent / f"{input_path.stem}_quality_report"
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {input_path}")
    data_object = read_mesh(input_path)
    leaves = list(iter_leaf_datasets(data_object))
    if not leaves:
        raise RuntimeError("No VTK datasets were found in the input mesh.")

    all_csv_file = None
    all_writer = None
    if not args.skip_all_cells_csv:
        all_csv_file = (output_dir / "all_volume_cells.csv").open("w", newline="", encoding="utf-8")
        all_writer = csv.DictWriter(all_csv_file, fieldnames=CSV_FIELDS)
        all_writer.writeheader()

    problem_csv_file = (output_dir / "problem_cells.csv").open("w", newline="", encoding="utf-8")
    problem_writer = csv.DictWriter(problem_csv_file, fieldnames=CSV_FIELDS)
    problem_writer.writeheader()

    type_counts: Counter = Counter()
    class_counts: Counter = Counter()
    type_class_counts: Dict[str, Counter] = defaultdict(Counter)
    qualities_by_type: Dict[str, List[float]] = defaultdict(list)
    issue_counts: Counter = Counter()
    total_volume_cells = 0
    ignored_non_volume_cells = 0
    worst_heap: List[Tuple[float, int, Dict[str, Any]]] = []
    worst_serial = 0
    annotated_grids: List[Any] = []

    try:
        for block_flat_index, block_name, grid in leaves:
            point_global_ids = find_id_array(
                grid.GetPointData(),
                ("GlobalNodeId", "GlobalNodeIds", "NodeId", "NodeIds", "vtkOriginalPointIds"),
            )
            cell_global_ids = find_id_array(
                grid.GetCellData(),
                (
                    "GlobalElementId",
                    "GlobalElementIds",
                    "ElementId",
                    "ElementIds",
                    "ObjectId",
                    "vtkOriginalCellIds",
                ),
            )

            annotated_grid = None
            annotated_arrays = None
            if not args.no_annotated_mesh:
                annotated_grid, annotated_arrays = make_annotated_grid(grid.GetPoints())
                try:
                    annotated_grid.GetPointData().ShallowCopy(grid.GetPointData())
                except Exception:
                    pass

            print(
                f"Inspecting block {block_flat_index} ({block_name}): "
                f"{grid.GetNumberOfCells()} cells"
            )
            for local_cell_id in range(grid.GetNumberOfCells()):
                cell = grid.GetCell(local_cell_id)
                if cell is None or cell.GetCellDimension() != 3:
                    ignored_non_volume_cells += 1
                    continue

                total_volume_cells += 1
                row, _internal = analyze_cell(
                    grid,
                    cell,
                    block_flat_index,
                    block_name,
                    local_cell_id,
                    point_global_ids,
                    cell_global_ids,
                    args,
                )

                cell_type = row["cell_type"]
                classification = row["classification"]
                type_counts[cell_type] += 1
                class_counts[classification] += 1
                type_class_counts[cell_type][classification] += 1
                quality = float(row["scaled_jacobian"])
                if finite(quality):
                    qualities_by_type[cell_type].append(quality)

                if row["repeated_node_id_flag"]:
                    issue_counts["repeated_node_ids"] += 1
                if row["coincident_point_flag"]:
                    issue_counts["coincident_points"] += 1
                if row["edge_degenerate_flag"]:
                    issue_counts["near_zero_edges"] += 1
                if row["face_degenerate_flag"]:
                    issue_counts["near_zero_faces"] += 1
                if row["volume_degenerate_flag"]:
                    issue_counts["near_zero_volumes"] += 1

                if all_writer is not None:
                    all_writer.writerow(row)
                if classification not in {"GOOD", "ACCEPTABLE"}:
                    problem_writer.writerow(row)

                if finite(quality) and args.worst_count > 0:
                    heap_entry = (-quality, worst_serial, dict(row))
                    worst_serial += 1
                    if len(worst_heap) < args.worst_count:
                        heapq.heappush(worst_heap, heap_entry)
                    elif quality < -worst_heap[0][0]:
                        heapq.heapreplace(worst_heap, heap_entry)

                if annotated_grid is not None and annotated_arrays is not None:
                    annotated_grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
                    append_annotated_values(annotated_arrays, row)

            if annotated_grid is not None and annotated_grid.GetNumberOfCells() > 0:
                annotated_grids.append(annotated_grid)
    finally:
        if all_csv_file is not None:
            all_csv_file.close()
        problem_csv_file.close()

    if total_volume_cells == 0:
        raise RuntimeError("No 3-D volume cells were found in the mesh.")

    worst_rows = sorted((entry[2] for entry in worst_heap), key=lambda row: row["scaled_jacobian"])

    summary_path = output_dir / "summary.txt"
    write_summary(
        summary_path,
        input_path,
        args,
        total_volume_cells,
        ignored_non_volume_cells,
        type_counts,
        class_counts,
        type_class_counts,
        qualities_by_type,
        issue_counts,
        worst_rows,
    )

    if not args.no_annotated_mesh and annotated_grids:
        append_filter = vtk.vtkAppendFilter()
        if hasattr(append_filter, "MergePointsOff"):
            append_filter.MergePointsOff()
        for annotated_grid in annotated_grids:
            append_filter.AddInputData(annotated_grid)
        append_filter.Update()

        annotated_path = output_dir / "annotated_volume_mesh.vtu"
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(str(annotated_path))
        writer.SetInputData(append_filter.GetOutput())
        if hasattr(writer, "SetDataModeToBinary"):
            writer.SetDataModeToBinary()
        if hasattr(writer, "SetCompressorTypeToZLib"):
            writer.SetCompressorTypeToZLib()
        if writer.Write() != 1:
            raise RuntimeError(f"Failed to write annotated mesh: {annotated_path}")

    print("")
    print(f"Report written to: {output_dir}")
    print(f"  {summary_path.name}")
    if not args.skip_all_cells_csv:
        print("  all_volume_cells.csv")
    print("  problem_cells.csv")
    if not args.no_annotated_mesh:
        print("  annotated_volume_mesh.vtu")
    print("")
    print("Classification totals:")
    for class_name in CLASS_CODES:
        print(f"  {class_name:<12} {class_counts.get(class_name, 0)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit("Interrupted by user.")
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
