"""
Calcula la permitividad dielectrica relativa del sistema a partir de las
fluctuaciones temporales de la polarizacion macroscopica (metodo de
fluctuacion-disipacion, sin campo externo aplicado):

    eps_ab(T) = delta_ab + V / (eps0 * kB * T) * [ <Pa Pb>_t - <Pa>_t <Pb>_t ]

donde P(t) es la polarizacion macroscopica (promedio espacial de las
N_SITIOS polarizaciones locales) en cada paso temporal, V es el volumen de
la supercelda simulada y T la temperatura en Kelvin de la corrida.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from leer_P_local import N_SITIOS, datos

DATA_CAMPO = Path(__file__).parent.parent / "Datos_permitividad" / "permittivity-e111-table-0.40-S0.dat"

EPS0 = 8.854187817e-12   # F/m
KB = 1.380649e-23        # J/K
A_RED = 4.0e-10          # m, constante de red de la celda unidad
V = N_SITIOS * A_RED**3  # m^3, volumen de la supercelda

U_P_UC_CM2_TO_C_M2 = 1e-2  # uC/cm^2 -> C/m^2


def permitividad_por_fluctuaciones(P, T, descartar=0):
    """Tensor de permitividad relativa (3,3) a partir de P(t).

    P: array (n_pasos, N_SITIOS, 3) en uC/cm^2
    T: temperatura en Kelvin
    descartar: cantidad de pasos iniciales a excluir (equilibracion)
    """
    P_macro = P[descartar:].mean(axis=1) * U_P_UC_CM2_TO_C_M2  # (n_pasos, 3), C/m^2
    P_prom = P_macro.mean(axis=0)
    cov = (P_macro.T @ P_macro) / P_macro.shape[0] - np.outer(P_prom, P_prom)
    eps = np.eye(3) + (V / (EPS0 * KB * T)) * cov
    return eps


# ──────────────────────────────────────────────
# Permitividad para cada temperatura disponible
# ──────────────────────────────────────────────
permitividad = {}
for T, d in sorted(datos.items()):
    eps = permitividad_por_fluctuaciones(d["P"], T)
    permitividad[T] = eps
    eps_prom = np.trace(eps) / 3
    print(
        f"T={T:4d} K   eps_xx={eps[0,0]:10.2f}  eps_yy={eps[1,1]:10.2f}  "
        f"eps_zz={eps[2,2]:10.2f}  eps_prom={eps_prom:10.2f}"
    )

# A partir de aca `permitividad` tiene, por cada T, el tensor eps (3,3):
#   permitividad[70] -> array (3,3) con eps_xx, eps_yy, eps_zz en la diagonal

# ──────────────────────────────────────────────
# Graficos: promedio vs T (izq.) y su inversa (der.)
# ──────────────────────────────────────────────
T_arr = np.array(sorted(permitividad.keys()))
eps_prom = np.array([np.trace(permitividad[T]) / 3 for T in T_arr])


def smooth3(y):
    """Promedio movil de 3 puntos (centrado); los extremos quedan sin suavizar."""
    smoothed = y.copy()
    smoothed[1:-1] = (y[:-2] + y[1:-1] + y[2:]) / 3
    return smoothed


eps_prom = smooth3(eps_prom)
inv_eps = 1 / eps_prom

# ──────────────────────────────────────────────
# Permitividad calculada con campo electrico aplicado (dato de referencia)
# ──────────────────────────────────────────────
T_campo, eps_campo, _ = np.loadtxt(DATA_CAMPO, skiprows=1, unpack=True)
inv_eps_campo = 1 / eps_campo

# ──────────────────────────────────────────────
# Ajuste lineal de Curie-Weiss (1/eps = (T - Tc)/C) a alta T
# ──────────────────────────────────────────────
N_PUNTOS_AJUSTE = 6      # cantidad de puntos de mayor T usados para el ajuste

pendiente, ordenada = np.polyfit(T_arr[-N_PUNTOS_AJUSTE:], inv_eps[-N_PUNTOS_AJUSTE:], 1)
ajuste = pendiente * T_arr + ordenada
T_curie = -ordenada / pendiente  # temperatura de Curie-Weiss (interseccion con eps->inf)

pendiente_campo, ordenada_campo = np.polyfit(T_campo[-N_PUNTOS_AJUSTE:], inv_eps_campo[-N_PUNTOS_AJUSTE:], 1)
ajuste_campo = pendiente_campo * T_campo + ordenada_campo
T_curie_campo = -ordenada_campo / pendiente_campo

print(f"\nAjuste Curie-Weiss (alta T, ultimos {N_PUNTOS_AJUSTE} puntos):")
print(f"  covarianza:    1/eps = {pendiente:.6e} * T + {ordenada:.6e}   T_Curie-Weiss = {T_curie:.2f} K")
print(f"  campo elec.:   1/eps = {pendiente_campo:.6e} * T + {ordenada_campo:.6e}   T_Curie-Weiss = {T_curie_campo:.2f} K")

COLOR_COV = 'black'
COLOR_CAMPO = 'tab:green'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
ax1b = ax1.twinx()
ax2b = ax2.twinx()

ax1.plot(T_arr, eps_prom, label=r"$\varepsilon_{prom}$ (covarianza)", marker='o', markersize=6, linewidth=2.5, color=COLOR_COV)
ax1b.plot(T_campo, eps_campo, label=r"$\varepsilon_{prom}$ (campo elec., S0)", marker='s', markersize=6, linewidth=2.5, color=COLOR_CAMPO)
ax2.plot(T_arr, inv_eps, label=r"$\varepsilon_{prom}$ (covarianza)", marker='o', markersize=6, linewidth=2.5, color=COLOR_COV)
ax2.plot(T_arr, ajuste, label="Ajuste Curie-Weiss (covarianza)", linestyle='--', linewidth=2.0, color=COLOR_COV)
ax2b.plot(T_campo, inv_eps_campo, label=r"$\varepsilon_{prom}$ (campo elec., S0)", marker='s', markersize=6, linewidth=2.5, color=COLOR_CAMPO)
ax2b.plot(T_campo, ajuste_campo, label="Ajuste Curie-Weiss (campo elec.)", linestyle='--', linewidth=2.0, color=COLOR_CAMPO)

ax1.set_ylabel(r"$\varepsilon_r$ (covarianza)", size=20, color=COLOR_COV)
ax1b.set_ylabel(r"$\varepsilon_r$ (campo elec.)", size=20, color=COLOR_CAMPO)
ax2.set_ylabel(r"$1/\varepsilon_r$ (covarianza)", size=20, color=COLOR_COV)
ax2b.set_ylabel(r"$1/\varepsilon_r$ (campo elec.)", size=20, color=COLOR_CAMPO)

for ax in (ax1, ax2):
    ax.set_xlabel(r"Temperature $[K]$", size=18)
    ax.tick_params(axis='y', colors=COLOR_COV)

for axb in (ax1b, ax2b):
    axb.tick_params(axis='y', colors=COLOR_CAMPO)

for ax in (ax1, ax2, ax1b, ax2b):
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=14, width=2.0, length=7.0)
    ax.tick_params(axis='both', which='minor', width=1.5, length=3.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

for ax, axb in ((ax1, ax1b), (ax2, ax2b)):
    ax.grid(True, linestyle='--', alpha=0.6)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = axb.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=13)

fig.tight_layout()
fig.savefig("Permitividad-covarianza.png", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
