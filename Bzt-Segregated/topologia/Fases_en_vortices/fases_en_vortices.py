# ------ SCRIPT PARAMS ---------
CONC = "0.40"
SEG = "S100"

# Grouped SOM files: classify every Ti cell into R/O/T/C for each temperature
GROUPED_DIR = f"/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Fracción_de_fase/MultiSOM/Grouped_files/{CONC}-{SEG}"
GROUPED_FILE = lambda i: f"{GROUPED_DIR}/Grouped_SOM_ss_xx_{i}.dump"
N_TEMPS = 50

# Vortex cell-membership files: which cells (celda) belong to each vortex
VORTICE_DIR = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/topologia/Fases_en_vortices"
VORTICES = {
    "agregado": f"{VORTICE_DIR}/Datos_vortice_agregado_0.csv",
    "intersticial": f"{VORTICE_DIR}/Datos_vortice_intersticial_0.csv",
}

# English title shown on each plot, keyed by the internal vortex name above
VORTEX_TITLE = {
    "agregado": "Aggregate Vortex",
    "intersticial": "Interstitial Vortex",
}

# identidad == 2 -> Ti, identidad == 1 -> Zr (excluded from this analysis)
TI_IDENTIDAD = 2

# Group id -> phase label (from Fracción_de_fase/group_hist_analysis.py group_dict)
GROUP_TO_PHASE = {0: "T", 1: "R", 2: "O", 3: "C"}
PHASE_ORDER = ["R", "O", "T", "C"]
PHASE_STYLE = {
    "R": dict(marker="s", color="red", label="R"),
    "O": dict(marker="^", color="blue", label="O"),
    "T": dict(marker="o", color="darkgreen", label="T"),
    "C": dict(marker="d", color="orange", label="C"),
}

# Temperature values: files 0-39 -> 10..400 K (step 10), files 40-49 -> 410..590 K (step 20)
T_VALS = list(range(10, 410, 10)) + list(range(410, 600, 20))
RESCALE_FACTOR = 1.48

OUTPUT_DATA_TEMPLATE = f"{VORTICE_DIR}/Fraccion_fases_vortice_{{name}}.csv"
OUTPUT_PLOT_TEMPLATE = f"{VORTICE_DIR}/Fraccion_fases_vortice_{{name}}.png"

# ------------------------------

import csv
import numpy as np
import matplotlib.pyplot as plt

# Plot style, matched to Fracción_de_fase/graficas-finales.py
Marker_size = 9
line_width = 2.5
Xlabel_size = 24
Ylabel_size = 24
Xlimit = (0, 800)
Xtick_step = 200
Ylimit = (-0.05, 1.1)
Ytick_step = 0.5
legend_size = 14
Grid_style = "--"
Grid_alpha = 0.0
Grid_width = 1.5
Spine_width = 1.5
Tick_major_width = 2.0
Tick_major_length = 7.0
Tick_minor_width = 1.5
Tick_minor_length = 3.0
Tick_labelsize = 16
Subtitles_size = 20
Subtitles_type = "bold"


def load_ti_celdas(csv_path):
    """Return the Ti-cell ids ('celda') referenced by a vortex file."""
    celdas = set()
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row["identidad"]) == TI_IDENTIDAD:
                celdas.add(int(row["celda"]))
    return celdas


def count_phases_in_dump(dump_path, celdas_by_vortex):
    """
    Single pass over a Grouped SOM dump file. Returns, for each vortex name,
    a dict {phase_label: count} tallied over the cells listed in celdas_by_vortex.
    """
    counts = {name: {p: 0 for p in PHASE_ORDER} for name in celdas_by_vortex}
    with open(dump_path) as f:
        for _ in range(9):
            next(f)  # skip LAMMPS dump header
        for line in f:
            parts = line.split()
            if not parts:
                continue
            group = int(round(float(parts[-2])))
            cell_id = int(parts[-1])
            phase = GROUP_TO_PHASE.get(group)
            if phase is None:
                continue  # group 4 -> Zr, not classified into a phase
            for name, celdas in celdas_by_vortex.items():
                if cell_id in celdas:
                    counts[name][phase] += 1
    return counts


def main():
    celdas_by_vortex = {name: load_ti_celdas(path) for name, path in VORTICES.items()}
    totals = {name: len(celdas) for name, celdas in celdas_by_vortex.items()}
    for name, total in totals.items():
        print(f"Vortice '{name}': {total} celdas de Ti")

    # fractions[name][phase] = list of fractions, one per temperature index
    fractions = {
        name: {p: [] for p in PHASE_ORDER} for name in VORTICES
    }

    for i in range(N_TEMPS):
        dump_path = GROUPED_FILE(i)
        counts = count_phases_in_dump(dump_path, celdas_by_vortex)
        for name in VORTICES:
            classified = sum(counts[name].values())
            if classified != totals[name]:
                print(
                    f"[aviso] T index {i} ({dump_path}): {classified}/{totals[name]} "
                    f"celdas de Ti del vortice '{name}' encontradas en el dump"
                )
            for phase in PHASE_ORDER:
                frac = counts[name][phase] / totals[name] if totals[name] else 0.0
                fractions[name][phase].append(frac)

    T_rescaled = [t * RESCALE_FACTOR for t in T_VALS]

    for name in VORTICES:
        # Save the underlying data
        data_path = OUTPUT_DATA_TEMPLATE.format(name=name)
        with open(data_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["T_K", "T_rescaled_K"] + [f"f_{p}" for p in PHASE_ORDER])
            for j, t in enumerate(T_VALS):
                writer.writerow(
                    [t, T_rescaled[j]] + [fractions[name][p][j] for p in PHASE_ORDER]
                )
        print(f"Datos guardados en {data_path}")

        # Plot (same style as Fracción_de_fase/graficas-finales.py)
        fig, ax = plt.subplots(figsize=(8, 6))
        for phase in PHASE_ORDER:
            style = PHASE_STYLE[phase]
            ax.plot(
                T_rescaled,
                fractions[name][phase],
                marker=style["marker"],
                linestyle="-",
                color=style["color"],
                markersize=Marker_size,
                linewidth=line_width,
                label=style["label"],
            )
        ax.set_title(VORTEX_TITLE[name], fontsize=Subtitles_size, fontweight=Subtitles_type)
        ax.set_xlim(Xlimit)
        ax.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
        ax.set_ylim(Ylimit)
        ax.set_yticks(np.arange(0, 1.1, Ytick_step))
        ax.set_xlabel(r"Temperature $[K]$", fontsize=Xlabel_size)
        ax.set_ylabel(r"Ti Phase Fraction", fontsize=Ylabel_size)
        ax.tick_params(axis="both", which="major", width=Tick_major_width,
                        length=Tick_major_length, labelsize=Tick_labelsize)
        ax.tick_params(axis="both", which="minor", width=Tick_minor_width,
                        length=Tick_minor_length)
        ax.minorticks_on()
        ax.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
        for spine in ax.spines.values():
            spine.set_linewidth(Spine_width)
        ax.legend(fontsize=legend_size, loc="upper left")
        fig.tight_layout()

        plot_path = OUTPUT_PLOT_TEMPLATE.format(name=name)
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Grafico guardado en {plot_path}")


if __name__ == "__main__":
    main()
