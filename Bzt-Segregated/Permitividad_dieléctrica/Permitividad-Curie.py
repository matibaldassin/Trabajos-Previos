#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

# Etiquetas
CONC = "0.40"
AGREGACION = ["S0", "S50", "S100"]   # una curva por agregación

# Factor de escala de temperatura (se aplica antes de graficar y de ajustar)
T_SCALE = 1.48

# Ventana del promedio móvil para suavizar 1/eps_r vs T (1 = sin suavizado)
SMOOTH_N = 3

# Cantidad de puntos (los de mayor temperatura) que entran en el ajuste de Curie-Weiss
N_FIT = 5

# Carpetas (relativas a la carpeta de este script)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "Datos_permitividad")


def smooth(y, n):
    """Promedio móvil centrado de n puntos; en los bordes promedia los puntos disponibles."""
    y = np.asarray(y, dtype=float)
    if n <= 1 or len(y) < 2:
        return y
    kernel = np.ones(n)
    suma = np.convolve(y, kernel, mode="same")
    cuenta = np.convolve(np.ones_like(y), kernel, mode="same")
    return suma / cuenta


def fit_curie_weiss(T, inv_eps):
    """
    Ajuste lineal 1/eps_r = (T - T0)/C.
    Devuelve (C, T0, slope, intercept, r2).
    """
    slope, intercept = np.polyfit(T, inv_eps, 1)
    C = 1.0 / slope
    T0 = -intercept / slope
    ss_res = np.sum((inv_eps - (slope * T + intercept)) ** 2)
    ss_tot = np.sum((inv_eps - inv_eps.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return C, T0, slope, intercept, r2


fig, (ax_eps, ax_inv) = plt.subplots(1, 2, figsize=(16, 6))
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

for i, agregacion in enumerate(AGREGACION):
    data_file = os.path.join(DATA_DIR, f"permittivity-table-{CONC}-{agregacion}.dat")
    if not os.path.exists(data_file):
        print(f"[{agregacion}] No se encontró el archivo de datos: {data_file}")
        continue

    data = np.loadtxt(data_file, usecols=(0, 1), comments="#")
    data = data[np.argsort(data[:, 0])]
    T = T_SCALE * data[:, 0]
    eps = data[:, 1]
    inv_eps = 1.0 / eps
    inv_eps_s = smooth(inv_eps, SMOOTH_N)
    c = colors[i % len(colors)]

    # ---- Permitividad relativa ----
    ax_eps.plot(T, eps, marker='o', linewidth=1.8, color=c, label=agregacion)

    # ---- Inversa de la permitividad + ajuste de Curie-Weiss ----
    ax_inv.plot(T, inv_eps, marker='o', linestyle='none', alpha=0.35, color=c)
    ax_inv.plot(T, inv_eps_s, linewidth=2.0, color=c, label=f"{agregacion} (suavizado n={SMOOTH_N})")

    n_fit = min(N_FIT, len(T))
    if n_fit < 2:
        print(f"[{agregacion}] No hay suficientes puntos para el ajuste.")
        continue
    T_fit, inv_fit = T[-n_fit:], inv_eps_s[-n_fit:]
    C, T0, slope, intercept, r2 = fit_curie_weiss(T_fit, inv_fit)
    print(f"[{agregacion}] C = {C:.4g} K, T0 = {T0:.1f} K, R2 = {r2:.4f} ({n_fit} puntos)")

    ax_inv.plot(T_fit, inv_fit, marker='s', markersize=7, linestyle='none',
                markerfacecolor='none', markeredgecolor=c)
    T_line = np.linspace(max(T0, T.min()), T.max(), 200)
    ax_inv.plot(T_line, slope * T_line + intercept, linestyle='--', linewidth=1.6, color=c,
                label=f"Curie-Weiss {agregacion}: C = {C:.3g} K, " + r"$T_0$" + f" = {T0:.0f} K, "
                      + r"$R^2$" + f" = {r2:.3f}")

ax_eps.set_xlabel("Temperatura [K]", fontsize=14, fontweight='bold')
ax_eps.set_ylabel("Permitividad Relativa " + r"$(\varepsilon_{r,111})$", fontsize=14, fontweight='bold')
ax_eps.set_title(f"Permitividad Relativa vs Temperatura, x = {CONC}", fontsize=15, fontweight='bold')

ax_inv.set_xlabel("Temperatura [K]", fontsize=14, fontweight='bold')
ax_inv.set_ylabel(r"$1/\varepsilon_{r,111}$", fontsize=14, fontweight='bold')
ax_inv.set_title(f"Inversa de la Permitividad, x = {CONC} (ajuste: {N_FIT} puntos)", fontsize=15, fontweight='bold')

for ax in (ax_eps, ax_inv):
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.minorticks_on()

fig.tight_layout()
plt.show()
