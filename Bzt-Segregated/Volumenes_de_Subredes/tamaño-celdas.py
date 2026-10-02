#!/usr/bin/env python3
"""
Compute per-cell edge lengths (a_x, a_y, a_z) from DL_POLY data using triclinic PBC,
filter out cells with any edge > max_edge, and report stats split by B-species (Ti vs Zr).

Inputs:
  - CONFIG file (only used to identify species labels; NOT used for box)
  - tabpol.dat (for each B-site center: the 29 indices including the 8 Ba/Ca corners)
    Header lines look like:  "<center_index>  29  <Ti|Zr>"
    where <center_index> is the CONFIG index of the B-core (e.g., 2, 12, 22, ...)
  - lattice file (e.g., entrada-vectores-de-red-200.dat), 3 lines:
        Ax Ay Az
        Bx By Bz
        Cx Cy Cz
  - One or more pos-promedios files ("id  Rx Ry Rz" per line)

Outputs (per pos-promedios file):
  - cell_lengths_<stem>_kept.csv      -> cells kept (used in averages; all edges ≤ max_edge)
  - cell_lengths_<stem>_excluded.csv  -> cells excluded (any edge > max_edge)
  - cell_lengths_summary_<stem>.txt   -> stats over KEPT cells:
        * overall means ± std
        * Ti-only means ± std
        * Zr-only means ± std
"""

import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
from typing import Dict, Tuple, List

# ---------------------------
# Parsing utilities
# ---------------------------

def read_config_labels(config_path: Path) -> Dict[int, str]:
    """Read CONFIG just to get {index -> label}."""
    with open(config_path, "r") as f:
        lines = f.readlines()
    labels: Dict[int, str] = {}
    i = 5  # after header
    while i < len(lines):
        label_line = lines[i].rstrip("\n")
        i += 2  # skip coordinate line
        parts = label_line.split()
        if len(parts) < 2 or not parts[-1].isdigit():
            m = re.match(r"^\s*([A-Za-z\-]+)\s+(\d+)\s*$", label_line)
            if not m:
                raise ValueError(f"Could not parse label/index line: {label_line!r}")
            label = m.group(1); idx = int(m.group(2))
        else:
            label = " ".join(parts[:-1]); idx = int(parts[-1])
        labels[idx] = label.strip()
    return labels

def read_posavg(pos_path: Path) -> np.ndarray:
    """Read pos-promedios file: returns an array pos with shape (max_id+1, 3) indexed by atom id."""
    max_id = 0
    coords: Dict[int, Tuple[float, float, float]] = {}
    with open(pos_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            idx = int(parts[0])
            x, y, z = map(float, parts[1:4])
            coords[idx] = (x, y, z)
            if idx > max_id:
                max_id = idx
    pos = np.zeros((max_id+1, 3), dtype=float)
    for k, v in coords.items():
        pos[k] = v
    return pos

def read_tabpol_blocks(tabpol_path: Path) -> List[Tuple[int, str, List[int]]]:
    """
    Read tabpol.dat and return a list of (center_index, B_species, indices[29]).
    center_index is the CONFIG index of the B-core (e.g., 2, 12, 22, ...).
    """
    blocks = []
    with open(tabpol_path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f]
    pat = re.compile(r"^\s*(\d+)\s+29\s+(Ti|Zr)\s*$")
    i = 0
    while i < len(lines):
        m = pat.match(lines[i])
        if m:
            center_idx = int(m.group(1))
            b_spec = m.group(2)
            idxs = []
            j = i + 1
            while j < len(lines) and len(idxs) < 29 and not pat.match(lines[j]):
                if lines[j].strip():
                    idxs.extend([int(x) for x in lines[j].split()])
                j += 1
            if len(idxs) != 29:
                raise ValueError(f"Center {center_idx} does not have 29 indices (got {len(idxs)}).")
            blocks.append((center_idx, b_spec, idxs))
            i = j
        else:
            i += 1
    if not blocks:
        raise ValueError("No center blocks found in tabpol.dat.")
    return blocks

def read_lattice(lattice_path: Path) -> np.ndarray:
    """
    Read 3x3 lattice matrix (columns = A, B, C) from a file.
    Uses ONLY the first 3 numeric tokens of the first 3 non-empty lines.
    Extra tokens (e.g. stresses, strain) are ignored.
    """
    vecs = []
    with open(lattice_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(("#", "!", ";", "//")):
                continue
            tokens = line.replace(",", " ").split()
            floats = []
            for tok in tokens:
                try:
                    floats.append(float(tok.replace("D", "E").replace("d", "e")))
                except ValueError:
                    continue
            if len(floats) >= 3:
                vecs.append(floats[:3])
            if len(vecs) == 3:
                break
    if len(vecs) != 3:
        raise ValueError(f"Could not extract 3 lattice vectors from {lattice_path}")
    A = np.array(vecs[0], dtype=float)
    B = np.array(vecs[1], dtype=float)
    C = np.array(vecs[2], dtype=float)
    return np.column_stack((A, B, C))

# ---------------------------
# Geometry helpers (triclinic MIC)
# ---------------------------

def wrap_minimum_image(dr_cart: np.ndarray, H: np.ndarray, Hinv: np.ndarray) -> np.ndarray:
    """Nearest-image in triclinic box: fractional -> wrap to [-0.5,0.5) -> back to Cartesian."""
    s = Hinv @ dr_cart
    s_wrapped = s - np.round(s)
    return H @ s_wrapped

def sign_with_tol(val: float, tol: float = 1e-3) -> int:
    if val > tol: return 1
    if val < -tol: return -1
    return 0

def edge_length_from_groups(rel_vecs: np.ndarray, axis: int, tol: float = 1e-3) -> float:
    """
    Compute edge length along 'axis' by grouping on signs of the other two axes,
    picking opposite-sign pairs within each group, averaging distances.
    Fallback: 2 * median(|component|).
    """
    other = [0, 1, 2]; other.remove(axis)
    g1, g2 = other
    groups = {}
    for v in rel_vecs:
        key = (sign_with_tol(v[g1], tol), sign_with_tol(v[g2], tol))
        groups.setdefault(key, []).append(v)
    lengths = []
    for vs in groups.values():
        if len(vs) < 2:
            continue
        best_sep = -1.0; best_pair = None
        for i in range(len(vs)):
            for j in range(i+1, len(vs)):
                if sign_with_tol(vs[i][axis], tol) * sign_with_tol(vs[j][axis], tol) == -1:
                    sep = abs(vs[i][axis] - vs[j][axis])
                    if sep > best_sep:
                        best_sep = sep; best_pair = (i, j)
        if best_pair is not None:
            vi, vj = vs[best_pair[0]], vs[best_pair[1]]
            lengths.append(float(np.linalg.norm(vi - vj)))
    if lengths:
        return float(np.mean(lengths))
    comps = np.abs(rel_vecs[:, axis])
    return float(2.0 * np.median(comps))

# ---------------------------
# Core computation
# ---------------------------

def compute_lengths_for_pos(labels, H, blocks, pos, tol_sign=1e-3):
    """
    For a single pos-promedios dataset:
      - labels: {index -> label} from CONFIG
      - H: 3x3 cell matrix with columns [A, B, C]
      - blocks: [(center_index, B_species, [29 indices]), ...]
      - pos: array indexed by atom id -> (x,y,z) in Cartesian Å
    Returns (DataFrame_all, skipped_missing).
    """
    Hinv = np.linalg.inv(H)
    rows = []
    skipped = 0
    for cell_id, (center_idx, b_spec, idxs) in enumerate(blocks, start=1):
        if center_idx >= pos.shape[0]:
            skipped += 1
            continue
        center_pos = pos[center_idx]
        corner_indices = [i for i in idxs if labels.get(i, "") in {"Ba", "Ca"}]
        if len(corner_indices) != 8 or any(i >= pos.shape[0] for i in corner_indices):
            skipped += 1
            continue
        rel = []
        for corner_idx in corner_indices:
            dr = pos[corner_idx] - center_pos
            rel.append(wrap_minimum_image(dr, H, Hinv))
        rel = np.array(rel, dtype=float)
        a_x = edge_length_from_groups(rel, axis=0, tol=tol_sign)
        a_y = edge_length_from_groups(rel, axis=1, tol=tol_sign)
        a_z = edge_length_from_groups(rel, axis=2, tol=tol_sign)
        rows.append((cell_id, center_idx, b_spec, a_x, a_y, a_z))
    df_all = pd.DataFrame(rows, columns=["cell_id", "center_index", "B_species", "a_x", "a_y", "a_z"]).sort_values("cell_id")
    return df_all, skipped

# ---------------------------
# CLI
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Compute Ba/Ca-defined edge lengths around Ti/Zr B-sites using triclinic PBC; filter out edges > threshold; report split stats for Ti vs Zr.")
    ap.add_argument("--config", "-c", required=True, help="Path to CONFIG file (for labels)")
    ap.add_argument("--tabpol", "-t", required=True, help="Path to tabpol.dat")
    ap.add_argument("--lattice", "-l", required=True, help="Path to lattice vectors file (3 lines: A, B, C)")
    ap.add_argument("--pos", "-p", nargs="+", required=True, help="One or more pos-promedios files")
    ap.add_argument("--outdir", "-o", default=".", help="Output directory (default: current)")
    ap.add_argument("--tol", type=float, default=1e-3, help="Sign tolerance for grouping (Å)")
    ap.add_argument("--max_edge", type=float, default=6.0, help="Max allowed edge length (Å) to KEEP a cell")
    args = ap.parse_args()

    config_path = Path(args.config)
    tabpol_path = Path(args.tabpol)
    lattice_path = Path(args.lattice)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    labels = read_config_labels(config_path)
    blocks = read_tabpol_blocks(tabpol_path)
    H = read_lattice(lattice_path)   # columns = A,B,C
    _ = np.linalg.inv(H)             # early check for invertibility

    for pos_path_str in args.pos:
        pos_path = Path(pos_path_str)
        pos = read_posavg(pos_path)

        df_all, skipped_missing = compute_lengths_for_pos(labels, H, blocks, pos, tol_sign=args.tol)

        # Filter: keep only cells where all edges ≤ max_edge
        mask_keep = (df_all["a_x"] <= args.max_edge) & (df_all["a_y"] <= args.max_edge) & (df_all["a_z"] <= args.max_edge)
        df_kept = df_all[mask_keep].copy()
        df_excl = df_all[~mask_keep].copy()

        stem = pos_path.stem
        kept_csv = outdir / f"cell_lengths_{stem}_kept.csv"
        excl_csv = outdir / f"cell_lengths_{stem}_excluded.csv"
        df_kept.to_csv(kept_csv, index=False, float_format="%.6f")
        df_excl.to_csv(excl_csv, index=False, float_format="%.6f")

        # Summary (over KEPT only), including Ti/Zr splits
        def grp_stats(df):
            if len(df) == 0:
                return None
            return {
                "a_x_mean": float(df["a_x"].mean()),
                "a_x_std": float(df["a_x"].std(ddof=1)) if len(df) > 1 else 0.0,
                "a_y_mean": float(df["a_y"].mean()),
                "a_y_std": float(df["a_y"].std(ddof=1)) if len(df) > 1 else 0.0,
                "a_z_mean": float(df["a_z"].mean()),
                "a_z_std": float(df["a_z"].std(ddof=1)) if len(df) > 1 else 0.0,
                "n": int(len(df)),
            }

        overall = grp_stats(df_kept)
        ti_only = grp_stats(df_kept[df_kept["B_species"] == "Ti"])
        zr_only = grp_stats(df_kept[df_kept["B_species"] == "Zr"])

        summary_path = outdir / f"cell_lengths_summary_{stem}.txt"
        with open(summary_path, "w") as f:
            f.write("Cell edge lengths (Å) computed from 8 Ba/Ca corners per Ti/Zr B-site center\n")
            f.write("Triclinic PBC with minimum-image; edges from quadrant-paired Ba/Ca atoms.\n\n")
            f.write(f"Cells computed: {len(df_all)}\n")
            f.write(f"Cells skipped (missing center or Ba/Ca corner data): {skipped_missing}\n")
            f.write(f"Cells excluded (> {args.max_edge:.3f} Å in any edge): {len(df_excl)}\n")
            f.write(f"Cells kept (used in averages): {len(df_kept)}\n\n")

            if overall:
                f.write("== Overall (KEPT) ==\n")
                f.write(f"a_x mean ± std: {overall['a_x_mean']:.6f} ± {overall['a_x_std']:.6f}  (n={overall['n']})\n")
                f.write(f"a_y mean ± std: {overall['a_y_mean']:.6f} ± {overall['a_y_std']:.6f}\n")
                f.write(f"a_z mean ± std: {overall['a_z_mean']:.6f} ± {overall['a_z_std']:.6f}\n\n")
            else:
                f.write("== Overall (KEPT) ==\nNo cells kept.\n\n")

            f.write("== Split by B-species (KEPT) ==\n")
            if ti_only:
                f.write(f"[Ti] a_x mean ± std: {ti_only['a_x_mean']:.6f} ± {ti_only['a_x_std']:.6f}  (n={ti_only['n']})\n")
                f.write(f"[Ti] a_y mean ± std: {ti_only['a_y_mean']:.6f} ± {ti_only['a_y_std']:.6f}\n")
                f.write(f"[Ti] a_z mean ± std: {ti_only['a_z_mean']:.6f} ± {ti_only['a_z_std']:.6f}\n")
            else:
                f.write("[Ti] No kept cells.\n")

            if zr_only:
                f.write(f"[Zr] a_x mean ± std: {zr_only['a_x_mean']:.6f} ± {zr_only['a_x_std']:.6f}  (n={zr_only['n']})\n")
                f.write(f"[Zr] a_y mean ± std: {zr_only['a_y_mean']:.6f} ± {zr_only['a_y_std']:.6f}\n")
                f.write(f"[Zr] a_z mean ± std: {zr_only['a_z_mean']:.6f} ± {zr_only['a_z_std']:.6f}\n")
            else:
                f.write("[Zr] No kept cells.\n")

        print(f"Wrote kept cells: {kept_csv}")
        print(f"Wrote excluded cells: {excl_csv}")
        print(f"Wrote summary: {summary_path}")

if __name__ == "__main__":
    main()
