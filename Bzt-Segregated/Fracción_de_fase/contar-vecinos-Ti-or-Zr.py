#!/usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator, NullLocator

#example
# python3 contar-vecinos-Ti-or-Zr-cod-0000004.py REVCON 5   (file radio de corte)

# === User-modifiable parameters ===
CUTOFF = 4.1                 # Å (used if not provided via CLI)
BAR_WIDTH = 0.85
SPECIES = "Ti"               # default species label(s); comma-separated allowed, e.g. "Ti,Ti-sh"

# B-site definition (cores only by default). Edit if you ever want to include shells.
B_SITE_LABELS = ["Ti", "Zr"]

# Range handling for x-axis in plots:
#   "full"     -> bins from 0..(N_species-1)
#   "observed" -> only values that actually appear
#   "auto"     -> use "full" if N_species-1 <= FULL_RANGE_MAX_SPAN, else "observed"
RANGE_MODE = "auto"
FULL_RANGE_MAX_SPAN = 200

# Bar annotations
ANNOTATE_PMF = True
MIN_ANNOTATE_FRAC = 0.03     # set to 0.0 to annotate all nonzero bars
ANNOTATE_COUNTS = True
MIN_ANNOTATE_COUNT = 2

# ---------------------------------------------------------------------------

def read_config(filename):
    """
    Robust DL_POLY CONFIG reader (levcfg = 0,1,2) that does NOT trust natms.
    Scans atoms until EOF and tolerates trailing blanks.
    Returns: title, levcfg, imcon, cell(3x3), atoms=[(label, pos(3), cfg_idx)]
    """
    def _maybe_int(x):
        try:
            return int(x)
        except Exception:
            return None

    with open(filename, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f]

    if len(lines) < 8:
        raise RuntimeError("CONFIG too short.")

    title  = lines[0]
    header = lines[1].split()
    if len(header) < 2:
        raise RuntimeError(f"Cannot parse levcfg/imcon from line 2: {lines[1]!r}")

    levcfg     = _maybe_int(header[0])
    imcon      = _maybe_int(header[1])
    natms_hint = _maybe_int(header[2]) if len(header) >= 3 else None
    if levcfg is None or imcon is None:
        raise RuntimeError(f"Cannot parse levcfg/imcon from line 2: {lines[1]!r}")

    # Cell vectors (3 lines)
    try:
        cell = np.array([
            list(map(float, lines[2].split())),
            list(map(float, lines[3].split())),
            list(map(float, lines[4].split())),
        ], dtype=float)
    except Exception as e:
        raise RuntimeError(f"Error parsing cell vectors (lines 3–5): {e}")

    atoms = []
    i = 5
    cfg_idx = 1

    while i < len(lines):
        # Skip trailing empty lines
        if not lines[i].strip():
            break

        # Label line
        toks = lines[i].split()
        label = toks[0] if toks else lines[i][:8].strip()
        i += 1
        if i >= len(lines):
            break  # truncated gracefully

        # Coords line (need at least 3 floats)
        ctoks = lines[i].split()
        if len(ctoks) < 3:
            break
        try:
            pos = np.array(list(map(float, ctoks[:3])), dtype=float)
        except Exception:
            break
        i += 1

        # Optional velocity / force lines based on levcfg
        if levcfg >= 1:
            if i >= len(lines): break
            i += 1
        if levcfg >= 2:
            if i >= len(lines): break
            i += 1

        atoms.append((label, pos, cfg_idx))
        cfg_idx += 1

    if natms_hint is not None and natms_hint != len(atoms):
        print(f"[warn] Header natms={natms_hint}, parsed {len(atoms)} atoms.")

    return title, levcfg, imcon, cell, atoms

def mic_delta(delta, cell, inv_cell):
    d_frac = delta @ inv_cell
    d_frac -= np.round(d_frac)    # wrap to [-0.5, 0.5)
    return d_frac @ cell

def count_species_and_bsite_neighbors(atoms, cell, cutoff, species_labels, bsite_labels):
    """
    For each atom with label in species_labels:
      - count how many OTHER atoms from the same set lie within cutoff (PBC+MIC)
      - count how many B-site atoms (labels in bsite_labels) lie within cutoff (excluding self)
      - return absolute CONFIG indices.
    Returns:
      species_positions (Nx3),
      species_counts (list N),            # same-species (e.g., Ti–Ti)
      species_cfg_indices (list N),
      bsite_counts (list N)               # Ti or Zr
    """
    species_set = set(species_labels)
    bsite_set = set(bsite_labels)

    species_positions, species_cfg_indices = [], []
    b_positions, b_cfg_indices = [], []

    for lab, pos, cfg_idx in atoms:
        if lab in bsite_set:
            b_positions.append(pos)
            b_cfg_indices.append(cfg_idx)
        if lab in species_set:
            species_positions.append(pos)
            species_cfg_indices.append(cfg_idx)

    if not species_positions:
        return np.zeros((0, 3)), [], [], []

    species_positions = np.array(species_positions, dtype=float)
    b_positions = np.array(b_positions, dtype=float)

    N = len(species_positions)
    M = len(b_positions)

    species_counts = [0] * N
    bsite_counts   = [0] * N

    inv_cell = np.linalg.inv(cell)

    # Same-species counts (Ti–Ti, etc.)
    for i in range(N):
        pos_i = species_positions[i]
        cnt_same = 0
        for j in range(N):
            if i == j:  # skip self
                continue
            delta = species_positions[j] - pos_i
            delta = mic_delta(delta, cell, inv_cell)
            if np.linalg.norm(delta) < cutoff:   # strict "<" (keep as-is)
                cnt_same += 1
        species_counts[i] = cnt_same

    # B-site counts (Ti or Zr), excluding self via config index
    for i in range(N):
        pos_i = species_positions[i]
        cfg_i = species_cfg_indices[i]
        cnt_b = 0
        for j in range(M):
            if b_cfg_indices[j] == cfg_i:
                continue  # exclude self
            delta = b_positions[j] - pos_i
            delta = mic_delta(delta, cell, inv_cell)
            if np.linalg.norm(delta) < cutoff:
                cnt_b += 1
        bsite_counts[i] = cnt_b

    return species_positions, species_counts, species_cfg_indices, bsite_counts

def choose_major_step(span):
    if span <= 16:  return 1
    if span <= 40:  return 2
    if span <= 80:  return 4
    if span <= 160: return 5
    return None

def style_axes_for_integer_x(ax, xs):
    x_min = int(xs.min()); x_max = int(xs.max())
    span = x_max - x_min + 1
    ax.set_xlim(x_min - 0.5, x_max + 0.5)

    step = choose_major_step(span)
    if step is None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
        ax.xaxis.set_minor_locator(NullLocator())
    else:
        ax.xaxis.set_major_locator(MultipleLocator(step))
        if span <= 200:
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        else:
            ax.xaxis.set_minor_locator(NullLocator())

    ax.grid(True, which='both', axis='x', linestyle=':', alpha=0.3)
    ax.grid(True, which='major', axis='y', linestyle=':', alpha=0.3)

def pmf_and_counts(counts_arr, total_n, range_mode=RANGE_MODE):
    max_possible = max(total_n - 1, 0)
    mode = range_mode.lower()
    if mode == "auto":
        use_full = (max_possible <= FULL_RANGE_MAX_SPAN)
    elif mode == "full":
        use_full = True
    elif mode == "observed":
        use_full = False
    else:
        use_full = True

    if use_full:
        xs = np.arange(0, max_possible + 1)
        raw = np.bincount(counts_arr, minlength=max_possible + 1)
        pmf = raw / total_n
        print(f"Histogram range: full 0..{max_possible} (N_species-1)")
    else:
        uniq = np.unique(counts_arr)
        xs = uniq
        full = np.bincount(counts_arr, minlength=int(uniq.max()) + 1)
        raw = full[xs]
        pmf = raw / total_n
        print(f"Histogram range: observed {xs.min()}..{xs.max()}")

    return xs, pmf, raw

def plot_pmf(xs, pmf, cutoff, out_png, species_tag):
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.bar(xs, pmf, width=BAR_WIDTH, edgecolor='black', linewidth=1.0)
    style_axes_for_integer_x(ax, xs)
    ax.set_xlabel(f"Number of {species_tag} neighbors within {cutoff:.1f} Å", fontweight='bold')
    ax.set_ylabel("Fraction of atoms", fontweight='bold')
    ax.set_title(f"Normalized histogram of {species_tag} neighbor counts", fontweight='bold')

    if pmf.size:
        ymax = float(np.max(pmf))
        ypad = 0.05 * max(ymax, 0.1)
        ax.set_ylim(0, ymax + 1.6 * ypad)

        def fmt_pct(p):
            if p == 0:
                return "0%"
            elif p < 0.005:
                return "<1%"
            return f"{p*100:.0f}%"

        for xi, pi in zip(xs, pmf):
            if (MIN_ANNOTATE_FRAC is None) or (pi >= MIN_ANNOTATE_FRAC):
                y_text = pi + ypad if pi > 0 else 0.4 * ypad
                ax.annotate(fmt_pct(pi), xy=(xi, pi), xytext=(0, 2), textcoords='offset points',
                            ha='center', va='bottom', fontsize=8, fontweight='bold', clip_on=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def plot_counts(xs, raw_counts, cutoff, out_png, species_tag):
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.bar(xs, raw_counts, width=BAR_WIDTH, edgecolor='black', linewidth=1.0)
    style_axes_for_integer_x(ax, xs)
    ax.set_xlabel(f"Number of {species_tag} neighbors within {cutoff:.1f} Å", fontweight='bold')
    ax.set_ylabel("Number of atoms", fontweight='bold')
    ax.set_title(f"Histogram of {species_tag} neighbor counts (raw)", fontweight='bold')

    if ANNOTATE_COUNTS and raw_counts.size:
        ymax = raw_counts.max()
        pad = 0.02 * (ymax if ymax > 0 else 1)
        for xi, ci in zip(xs, raw_counts):
            if ci >= MIN_ANNOTATE_COUNT:
                ax.text(xi, ci + pad, f"{int(ci)}", ha='center', va='bottom',
                        fontsize=8, fontweight='bold')

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def main():
    # CLI:
    #   1) CONFIG path (optional)
    #   2) cutoff Å (optional, default CUTOFF)
    #   3) species labels (comma-separated) e.g. "Ti" or "Ti,Ti-sh" (optional)
    config_file = "CONFIG"
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else CUTOFF
    species_arg = sys.argv[3] if len(sys.argv) > 3 else SPECIES
    species_labels = [s.strip() for s in species_arg.split(",") if s.strip()]
    species_tag = "+".join(species_labels)

    # Read CONFIG and count neighbors
    title, levcfg, imcon, cell, atoms = read_config(config_file)
    print(f"Parsed {len(atoms)} atoms (levcfg={levcfg}, imcon={imcon}).")

    positions, counts_same, cfg_indices, counts_bsite = count_species_and_bsite_neighbors(
        atoms, cell, cutoff, species_labels, B_SITE_LABELS
    )
    total_n = len(positions)

    # Output data (with CONFIG index in col 3 and B-site count in col 4)
    out_dat = f"vecinos-{species_tag}.dat"
    with open(out_dat, 'w') as f:
        f.write(f"# Title: {title}\n")
        f.write(f"# levcfg: {levcfg}  imcon: {imcon}  (analysis uses full 3D PBC/MIC)\n")
        f.write(f"# Neighbor counts for species = {species_labels}  (cutoff = {cutoff} Å)\n")
        f.write(f"# B-site labels considered: {B_SITE_LABELS}  (cores only)\n")
        f.write(f"# Total atoms of selected species: {total_n}\n")
        f.write("# Columns: 1) species_index (1..N_species), "
                "2) neighbors_within_cutoff_same_species, "
                "3) config_atom_index (1..N_atoms), "
                "4) neighbors_within_cutoff_B_site\n")
        for i, (c_same, cfg_idx, c_b) in enumerate(zip(counts_same, cfg_indices, counts_bsite), start=1):
            f.write(f"{i:5d}   {c_same:3d}   {cfg_idx:7d}   {c_b:3d}\n")

    print(f"Found {total_n} atoms of {species_labels}. Wrote neighbor counts to '{out_dat}'.")

    if total_n == 0:
        print("No atoms of the requested species found; skipping plots.")
        return

    # Sanity checks (for same-species counts, used in plots)
    counts_arr = np.array(counts_same, dtype=int)
    min_neighbors = int(counts_arr.min())
    max_neighbors = int(counts_arr.max())
    unique, freq = np.unique(counts_arr, return_counts=True)
    dist = dict(zip(unique.tolist(), freq.tolist()))
    max_possible = max(total_n - 1, 0)

    with open(out_dat, 'a') as f:
        f.write(f"# Min neighbors (same species, observed): {min_neighbors}\n")
        f.write(f"# Max neighbors (same species, observed): {max_neighbors}\n")
        f.write(f"# Max neighbors (same species, theoretical N-1): {max_possible}\n")

    print(f"Min neighbors (same species) = {min_neighbors}, Max = {max_neighbors}")
    print("Neighbor distribution (same species, k: count):", dist)
    print(f"Max neighbors (same species, theoretical) = {max_possible}")

    # Build xs / pmf / raw with safe range selection (same-species counts)
    xs, pmf, raw_counts = pmf_and_counts(counts_arr, total_n)

    # Plots (for same-species counts)
    out_pmf = f"histograma-{species_tag}-pmf-{cutoff:.1f}A-{RANGE_MODE}.png"
    plot_pmf(xs, pmf, cutoff, out_pmf, species_tag)
    print(f"Saved normalized plot to '{out_pmf}'.")

    out_raw = f"histograma-{species_tag}-counts-{cutoff:.1f}A-{RANGE_MODE}.png"
    plot_counts(xs, raw_counts, cutoff, out_raw, species_tag)
    print(f"Saved raw-count plot to '{out_raw}'.")

if __name__ == "__main__":
    main()
