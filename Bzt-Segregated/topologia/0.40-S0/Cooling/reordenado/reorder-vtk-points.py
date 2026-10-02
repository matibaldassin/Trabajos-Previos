#!/usr/bin/env python3
"""
Reorder VTK legacy STRUCTURED_GRID points by (x, y, z) ascending
and reorder associated POINT_DATA arrays accordingly:
- SCALARS polar_module double
- SCALARS identidad double
- VECTORS polar_vector double

Usage:
  python reorder_vtk_points.py input.vtk -o output.vtk
  python reorder_vtk_points.py Pol-*.vtk --inplace
"""

import argparse
import re
from pathlib import Path
import numpy as np
import sys
from typing import Tuple, List

def parse_structured_grid_legacy(text: str):
    lines = text.splitlines()

    # Basic checks
    if not lines or not lines[0].lstrip().startswith("# vtk"):
        raise ValueError("Not a VTK legacy file (missing '# vtk' header).")
    if "DATASET STRUCTURED_GRID" not in text:
        raise ValueError("Expected 'DATASET STRUCTURED_GRID' dataset.")

    # DIMENSIONS
    dim_line_idx = None
    for i, l in enumerate(lines):
        if l.strip().startswith("DIMENSIONS"):
            dim_line_idx = i
            break
    if dim_line_idx is None:
        raise ValueError("DIMENSIONS not found.")
    dims = list(map(int, lines[dim_line_idx].split()[1:4]))
    nx, ny, nz = dims
    n_expected = nx * ny * nz

    # Points header
    points_hdr_idx = None
    for i, l in enumerate(lines):
        if re.match(r"^\s*Points\s+\d+\s+double", l):
            points_hdr_idx = i
            break
    if points_hdr_idx is None:
        raise ValueError("Points header not found.")
    n_points = int(lines[points_hdr_idx].split()[1])
    if n_points != n_expected:
        raise ValueError(f"Point count {n_points} != DIMENSIONS product {n_expected}")

    # Gather point lines until POINT_DATA
    point_start = points_hdr_idx + 1
    point_end = None
    for i in range(point_start, len(lines)):
        if lines[i].strip().startswith("POINT_DATA"):
            point_end = i
            break
    if point_end is None:
        point_end = len(lines)

    # Parse points (triplets)
    pts = []
    for l in lines[point_start:point_end]:
        if not l.strip():
            continue
        vals = l.split()
        for j in range(0, len(vals), 3):
            trip = vals[j:j+3]
            if len(trip) == 3:
                pts.append(tuple(map(float, trip)))
    pts = np.array(pts, dtype=float)
    if pts.shape != (n_points, 3):
        raise ValueError(f"Parsed {pts.shape[0]} points but expected {n_points}")

    # POINT_DATA section
    pd_idx = None
    for i in range(point_end, len(lines)):
        if lines[i].strip().startswith("POINT_DATA"):
            pd_idx = i
            break
    if pd_idx is None:
        raise ValueError("POINT_DATA not found.")
    n_pd = int(lines[pd_idx].split()[1])
    if n_pd != n_points:
        raise ValueError("POINT_DATA count does not match number of points.")

    # SCALARS polar_module
    sc_polar_idx = None
    for i in range(pd_idx + 1, len(lines)):
        if re.match(r"^\s*SCALARS\s+polar_module\s+double", lines[i]):
            sc_polar_idx = i
            break
    if sc_polar_idx is None:
        raise ValueError("SCALARS polar_module double not found.")
    if not lines[sc_polar_idx + 1].strip().startswith("LOOKUP_TABLE"):
        raise ValueError("Missing LOOKUP_TABLE after polar_module.")
    pol_start = sc_polar_idx + 2

    # SCALARS identidad
    idn_idx = None
    for i in range(pol_start, len(lines)):
        if re.match(r"^\s*SCALARS\s+identidad\s+double", lines[i]):
            idn_idx = i
            break
    if idn_idx is None:
        raise ValueError("SCALARS identidad double not found.")
    if not lines[idn_idx + 1].strip().startswith("LOOKUP_TABLE"):
        raise ValueError("Missing LOOKUP_TABLE after identidad.")
    idn_start = idn_idx + 2

    # VECTORS polar_vector
    vec_idx = None
    for i in range(idn_start, len(lines)):
        if re.match(r"^\s*VECTORS\s+polar_vector\s+double", lines[i]):
            vec_idx = i
            break
    if vec_idx is None:
        raise ValueError("VECTORS polar_vector double not found.")
    vec_start = vec_idx + 1

    # Parse polar_module values
    pol_vals = []
    for l in lines[pol_start:idn_idx]:
        if l.strip():
            pol_vals.extend(map(float, l.split()))
    polar_module = np.array(pol_vals, dtype=float)
    if polar_module.shape[0] != n_points:
        raise ValueError("polar_module count mismatch.")

    # Parse identidad values
    idn_vals = []
    for l in lines[idn_start:vec_idx]:
        if l.strip():
            idn_vals.extend(map(float, l.split()))
    identidad = np.array(idn_vals, dtype=float)
    if identidad.shape[0] != n_points:
        raise ValueError("identidad count mismatch.")

    # Parse vectors (triplets)
    vecs = []
    for l in lines[vec_start:]:
        if not l.strip():
            continue
        vals = l.split()
        for j in range(0, len(vals), 3):
            trip = vals[j:j+3]
            if len(trip) == 3:
                vecs.append(tuple(map(float, trip)))
    vectors = np.array(vecs, dtype=float)
    if vectors.shape != (n_points, 3):
        raise ValueError("polar_vector count mismatch.")

    # Return pieces and indices to reconstruct header
    return {
        "lines": lines,
        "dims": (nx, ny, nz),
        "n_points": n_points,
        "points_hdr_idx": points_hdr_idx,
        "points": pts,
        "polar_module": polar_module,
        "identidad": identidad,
        "vectors": vectors,
        "pd_idx": pd_idx,
    }

def reorder_by_xyz(points: np.ndarray) -> np.ndarray:
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    # numpy.lexsort: last key is primary => (z,y,x) gives primary x, then y, then z
    return np.lexsort((z, y, x))

def write_structured_grid_legacy(out_path: Path, parsed, order: np.ndarray):
    lines = parsed["lines"]
    n_points = parsed["n_points"]
    points_hdr_idx = parsed["points_hdr_idx"]

    pts_sorted = parsed["points"][order]
    polar_module_sorted = parsed["polar_module"][order]
    identidad_sorted = parsed["identidad"][order]
    vectors_sorted = parsed["vectors"][order]

    with out_path.open("w", encoding="utf-8") as out:
        # copy original header up to Points line
        for i in range(points_hdr_idx):
            out.write(lines[i] + "\n")

        out.write(f"Points {n_points}  double\n")
        for p in pts_sorted:
            out.write(f"{p[0]:8.3f}    {p[1]:7.3f}    {p[2]:7.3f}\n")

        out.write(f"POINT_DATA {n_points}\n")

        out.write(" SCALARS polar_module double\n")
        out.write(" LOOKUP_TABLE default\n")
        for i in range(0, n_points, 10):
            chunk = polar_module_sorted[i:i+10]
            out.write(" " + "  ".join(f"{v:.2f}" for v in chunk) + "\n")

        out.write(" SCALARS identidad double\n")
        out.write(" LOOKUP_TABLE default\n")
        for i in range(0, n_points, 16):
            chunk = identidad_sorted[i:i+16]
            out.write(" " + " ".join(f"{v:.1f}" for v in chunk) + "\n")

        out.write(" VECTORS polar_vector double    \n")
        for v in vectors_sorted:
            out.write(f" {v[0]:6.2f}    {v[1]:6.2f}    {v[2]:6.2f}\n")

def process_file(in_path: Path, out_path: Path):
    text = in_path.read_text(encoding="utf-8", errors="ignore")
    parsed = parse_structured_grid_legacy(text)
    order = reorder_by_xyz(parsed["points"])
    write_structured_grid_legacy(out_path, parsed, order)

def main():
    ap = argparse.ArgumentParser(description="Reorder STRUCTURED_GRID VTK points by (x,y,z) and keep point-data aligned.")
    ap.add_argument("inputs", nargs="+", help="Input VTK file(s)")
    ap.add_argument("-o", "--output", help="Output file (only when processing a single input)")
    ap.add_argument("--inplace", action="store_true", help="Overwrite inputs in-place")
    ap.add_argument("--suffix", default="_sorted", help="When not in-place, append this suffix before extension (default: _sorted)")
    args = ap.parse_args()

    if args.output and (len(args.inputs) != 1):
        ap.error("When -o/--output is used, exactly one input must be provided.")

    for in_name in args.inputs:
        in_path = Path(in_name)
        if not in_path.exists():
            print(f"[skip] {in_path} not found", file=sys.stderr)
            continue

        if args.inplace:
            tmp_out = in_path.with_suffix(in_path.suffix + ".tmp")
            process_file(in_path, tmp_out)
            tmp_out.replace(in_path)
            print(f"[ok] Reordered in-place: {in_path}")
        else:
            if args.output:
                out_path = Path(args.output)
            else:
                out_path = in_path.with_name(in_path.stem + args.suffix + in_path.suffix)
            process_file(in_path, out_path)
            print(f"[ok] {in_path} -> {out_path}")

if __name__ == "__main__":
    main()
