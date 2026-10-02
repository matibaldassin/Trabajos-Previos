#!/usr/bin/env python3
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import os

# Etiquetas
CONC = "0.20"
AGREGACION = "S100"

# File names (relativos a la carpeta de este script)
BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE   = os.path.join(BASE_DIR, "Datos_polarizaciones", f"polmacro-{CONC}-{AGREGACION}")
TABLE_111    = os.path.join(BASE_DIR, "Datos_permitividad", f"permittivity-e111-table-{CONC}-{AGREGACION}.dat")
PLOT_111     = f"permittivity-e111-{CONC}-{AGREGACION}.png"

# Unit conversions
U_P_UC_CM2_TO_C_M2 = 1e-2   # μC/cm^2 -> C/m^2
U_E_KV_MM_TO_V_M   = 1e6    # kV/mm   -> V/m
EPS0 = 8.854187817e-12      # F/m

# Regex to parse "5kV_mm"
field_re = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*kV_mm\s*$")

# Corrección de etiquetas: en los archivos figura "10kV_mm" pero el campo real es 100 kV/mm
FIELD_FIX_KV_MM = {10.0: 10.0}

# Data holders by temperature
by_T_111 = defaultdict(list)   # (E, P_parallel)
by_T_Px  = defaultdict(list)   # (E, Px)
by_T_Py  = defaultdict(list)   # (E, Py)
by_T_Pz  = defaultdict(list)   # (E, Pz)

with open(INPUT_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        # Expect: T, "EkV_mm", Px, Py, Pz, |P|
        if len(parts) < 6:
            continue

        try:
            T_str, E_str = parts[0], parts[1]
            Px_str, Py_str, Pz_str = parts[2], parts[3], parts[4]

            T = int(float(T_str))

            m = field_re.match(E_str)
            if not m:
                continue
            E_kV_mm = float(m.group(1))
            E_kV_mm = FIELD_FIX_KV_MM.get(E_kV_mm, E_kV_mm)
            E_V_m = E_kV_mm * U_E_KV_MM_TO_V_M

            Px = float(Px_str) * U_P_UC_CM2_TO_C_M2
            Py = float(Py_str) * U_P_UC_CM2_TO_C_M2
            Pz = float(Pz_str) * U_P_UC_CM2_TO_C_M2

            # Projection onto [111]
            P_parallel = (Px + Py + Pz) / np.sqrt(3.0)

            by_T_111[T].append((E_V_m, P_parallel))
            by_T_Px[T].append((E_V_m, Px))
            by_T_Py[T].append((E_V_m, Py))
            by_T_Pz[T].append((E_V_m, Pz))

        except ValueError:
            continue

def fit_eps_vs_E(data_pairs):
    """
    data_pairs: list of (E, P_component) in SI.
    Returns (eps_r, slope, intercept) where slope = dP/dE [F/m], eps_r = 1 + slope/eps0
    """
    data_pairs = sorted(data_pairs, key=lambda x: x[0])
    E = np.array([d[0] for d in data_pairs], dtype=float)
    P = np.array([d[1] for d in data_pairs], dtype=float)
    slope, intercept = np.polyfit(E, P, 1)
    eps_r = 1.0 + slope / EPS0
    return eps_r, slope, intercept

# ---- [111] results ----
fit_details_111 = []  # (T, eps_r_111, slope_111)
for T in sorted(by_T_111.keys()):
    if len(by_T_111[T]) < 2:
        continue
    eps_r_111, slope_111, _ = fit_eps_vs_E(by_T_111[T])
    fit_details_111.append((T, eps_r_111, slope_111))

# Save [111] table (as before)
with open(TABLE_111, "w") as fout:
    fout.write("# T(K)    exx(rel)        slope_dP_dE_SI(C_per_m2_per_V_per_m)\n")
    for T, epsr, slope in sorted(fit_details_111):
        fout.write(f"{T:6d}  {epsr:12.1f}    {slope:.6e}\n")

# Plot [111]
if fit_details_111:
    fit_details_111.sort(key=lambda x: x[0])
    Ts_111   = 1.45*np.array([t for t, _, _ in fit_details_111])
    epsr_111 = [e for _, e, _ in fit_details_111]

    plt.figure(figsize=(10,6))
    plt.plot(Ts_111, epsr_111, marker='o', linewidth=1.8)
    plt.xlabel("Temperatura [K]", fontsize=14, fontweight='bold')
    plt.ylabel("Permitividad Relativa "+r"$(\varepsilon_{r,111})$", fontsize=14, fontweight='bold')
    plt.title(f"Permitividad Relativa vs Temperatura, x = {CONC}, {AGREGACION} ", fontsize=16, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()
else:
    print("No valid [111] results to process.")

