#!/usr/bin/env python3
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import os

# Etiquetas
CONC = "0.40"
AGREGACION = ["S0", "S50"]   # se procesa cada agregación de la lista

# Valores de campo eléctrico [kV/mm] a usar en el ajuste (campo aplicado en la dirección [111])
FIELDS_KV_MM = [0.0, 3.0, 5.0, 7.0, 10.0]

# Factor de escala de temperatura (solo para graficar)
T_SCALE = 1.48

# Carpetas (relativas a la carpeta de este script)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Unit conversions
U_P_UC_CM2_TO_C_M2 = 1e-2   # μC/cm^2 -> C/m^2
U_E_KV_MM_TO_V_M   = 1e6    # kV/mm   -> V/m
EPS0 = 8.854187817e-12      # F/m

# Regex to parse "5kV_mm"
field_re = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*kV_mm\s*$")


def read_polarization(input_file):
    """Devuelve {T: [(E [V/m], P_[111] [C/m^2])]} usando solo los campos de FIELDS_KV_MM."""
    by_T_111 = defaultdict(list)
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            # Expect: T, "EkV_mm", Px, Py, Pz, |P|
            if len(parts) < 6:
                continue

            try:
                T = int(float(parts[0]))

                m = field_re.match(parts[1])
                if not m:
                    continue
                E_kV_mm = float(m.group(1))
                if not any(np.isclose(E_kV_mm, fld) for fld in FIELDS_KV_MM):
                    continue
                E_V_m = E_kV_mm * U_E_KV_MM_TO_V_M  # campo aplicado a lo largo de [111]

                Px = float(parts[2]) * U_P_UC_CM2_TO_C_M2
                Py = float(parts[3]) * U_P_UC_CM2_TO_C_M2
                Pz = float(parts[4]) * U_P_UC_CM2_TO_C_M2

                # Projection onto [111]
                by_T_111[T].append((E_V_m, (Px + Py + Pz) / np.sqrt(3.0)))
            except ValueError:
                continue
    return by_T_111


def fit_eps_vs_E(data_pairs):
    """
    data_pairs: list of (E, P_component) in SI.
    Returns (eps_r, slope, intercept, r2) where slope = dP/dE [F/m], eps_r = 1 + slope/eps0
    """
    data_pairs = sorted(data_pairs, key=lambda x: x[0])
    E = np.array([d[0] for d in data_pairs], dtype=float)
    P = np.array([d[1] for d in data_pairs], dtype=float)
    slope, intercept = np.polyfit(E, P, 1)
    eps_r = 1.0 + slope / EPS0
    ss_res = np.sum((P - (slope * E + intercept)) ** 2)
    ss_tot = np.sum((P - P.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return eps_r, slope, intercept, r2


def process(agregacion):
    input_file = os.path.join(BASE_DIR, "Datos_polarizaciones", f"polmacro-{CONC}-{agregacion}")
    table_111  = os.path.join(BASE_DIR, "Datos_permitividad", f"permittivity-table-{CONC}-{agregacion}.dat")
    ref_111    = os.path.join(BASE_DIR, "Datos_permitividad", f"permittivity-e111-table-{CONC}-{agregacion}.dat")
    fits_dir   = os.path.join(BASE_DIR, "Graficas_EvsP", f"{CONC}-{agregacion}")

    if not os.path.exists(input_file):
        print(f"[{agregacion}] No se encontró el archivo de datos: {input_file}")
        return

    by_T_111 = read_polarization(input_file)

    # ---- [111] results ----
    fit_details_111 = []  # (T, eps_r_111, slope_111, r2)
    for T in sorted(by_T_111.keys()):
        if len(by_T_111[T]) < 2:
            continue
        eps_r_111, slope_111, _, r2_111 = fit_eps_vs_E(by_T_111[T])
        fit_details_111.append((T, eps_r_111, slope_111, r2_111))

    if not fit_details_111:
        print(f"[{agregacion}] No valid [111] results to process.")
        return

    # Save [111] table
    os.makedirs(os.path.dirname(table_111), exist_ok=True)
    with open(table_111, "w") as fout:
        fout.write("# T(K)    eps_111(rel)        slope_dP_dE_SI(C_per_m2_per_V_per_m)    R2\n")
        for T, epsr, slope, r2 in sorted(fit_details_111):
            fout.write(f"{T:6d}  {epsr:12.1f}    {slope:.6e}    {r2:.4f}\n")

    # Gráficas E vs P con el ajuste lineal, una por temperatura (carpeta aparte)
    os.makedirs(fits_dir, exist_ok=True)
    for T in sorted(by_T_111.keys()):
        if len(by_T_111[T]) < 2:
            continue
        pairs = sorted(by_T_111[T], key=lambda x: x[0])
        E = np.array([p[0] for p in pairs])
        P = np.array([p[1] for p in pairs])
        slope, intercept = np.polyfit(E, P, 1)
        E_fit = np.linspace(E.min(), E.max(), 100)
        ss_res = np.sum((P - (slope * E + intercept)) ** 2)
        ss_tot = np.sum((P - P.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(E / 1e6, P, "o", label="Datos")
        ax.plot(E_fit / 1e6, slope * E_fit + intercept, "-",
                label=f"Ajuste lineal\n" + r"$\varepsilon_r$" + f" = {1.0 + slope / EPS0:.1f}\n" + r"$R^2$" + f" = {r2:.4f}")
        ax.set_xlabel("E [kV/mm]", fontsize=14, fontweight='bold')
        ax.set_ylabel(r"$P_{111}$ [C/m$^2$]", fontsize=14, fontweight='bold')
        ax.set_title(f"E vs P, T = {T_SCALE*T:.1f} K, x = {CONC}, {agregacion}", fontsize=14, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(fits_dir, f"EvsP-{CONC}-{agregacion}-T{T}.png"), dpi=200)
        plt.close(fig)

    # Plot [111] (se muestra al final, junto con los de las otras agregaciones)
    fit_details_111.sort(key=lambda x: x[0])
    Ts_111   = T_SCALE*np.array([t for t, _, _, _ in fit_details_111])
    epsr_111 = [e for _, e, _, _ in fit_details_111]

    plt.figure(figsize=(10, 6))
    plt.plot(Ts_111, epsr_111, marker='o', linewidth=1.8, label="Campos seleccionados")

    # Comparación con la tabla de referencia permittivity-e111-table
    if os.path.exists(ref_111):
        ref = np.loadtxt(ref_111, usecols=(0, 1), comments="#")
        plt.plot(T_SCALE*ref[:, 0], ref[:, 1], marker='s', linewidth=1.8, linestyle='--',
                 label="permittivity-e111-table")
    else:
        print(f"[{agregacion}] No se encontró la tabla de referencia: {ref_111}")
    plt.legend(fontsize=12)
    plt.xlabel("Temperatura [K]", fontsize=14, fontweight='bold')
    plt.ylabel("Permitividad Relativa "+r"$(\varepsilon_{r,111})$", fontsize=14, fontweight='bold')
    plt.title(f"Permitividad Relativa vs Temperatura, x = {CONC}, {agregacion} ", fontsize=16, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.minorticks_on()
    plt.tight_layout()


for agregacion in AGREGACION:
    process(agregacion)

plt.show()
