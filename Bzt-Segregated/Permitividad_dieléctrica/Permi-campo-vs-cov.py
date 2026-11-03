#!/usr/bin/env python3
"""
Permitividad de BZT x=0.40 (S0, S50): campo electrico vs covarianza de P.

  - Campo: Datos_permitividad/permittivity-table-0.40-<S>.dat (salida de Permitividad.py,
    ajuste lineal P_111 vs E con los campos de FIELDS_KV_MM).
  - Covarianza (sin campo, V(T) de Parám_de_Red, T de simulacion): eps_ab = delta_ab + V/(eps0 kB T) * (<Pa Pb> - <Pa><Pb>),
    proyectada en [111]: eps_111 = n^T eps n, n = (1,1,1)/sqrt(3).

Los P-locales se leen de P-local-temperatura/Zr-0.40-<S>/P-local-<T> (T de simulacion,
sin escalar). Se cachea P_macro(t) de cada archivo en "Permitividad con covarianza/cache_Pmacro".
"""
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CONC = "0.40"
CORRIDAS = ["S0", "S50"]

BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
CARPETA_PLOCAL = Path("/home/matias/Escritorio/P-local-temperatura")
DATOS_CAMPO = BASE_DIR / "Datos_permitividad"
CACHE = BASE_DIR / "Permitividad con covarianza" / "cache_Pmacro"
CACHE.mkdir(parents=True, exist_ok=True)

N_SITIOS = 21952
EPS0 = 8.854187817e-12      # F/m
KB = 1.380649e-23           # J/K
A3_TO_M3 = 1e-30            # A^3 -> m^3
DIR_PARAM = BASE_DIR.parent / "Parám_de_Red" / "Parámetros-de-red-DL-POLY"
U_P_UC_CM2_TO_C_M2 = 1e-2   # uC/cm^2 -> C/m^2
T_SCALE = 1.48              # factor de escala de T (solo para graficar)
T_MAX = 590                 # T de simulacion maxima a usar [K] (rango de Parám_de_Red)
N_SUAVIZADO = 7             # puntos de la media movil (impar), solo para graficar
DESCARTAR = 0               # pasos iniciales a excluir (equilibracion)
N111 = np.ones(3) / np.sqrt(3)


def p_macro(path):
    """P macroscopico (n_pasos, 3) en uC/cm^2 de un archivo P-local."""
    path = Path(path)
    cache = CACHE / f"{path.parent.name}-{path.name}.npy"
    if cache.exists() and cache.stat().st_mtime > path.stat().st_mtime:
        return np.load(cache)
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["i", "x", "y", "z"],
                     dtype=str, engine="c")
    stride = N_SITIOS + 1
    if len(df) % stride:
        raise ValueError(f"{path}: {len(df)} filas no es multiplo de {stride}")
    # las filas encabezado (paso temporal) tienen solo la primera columna
    mask_enc = df["x"].isna().to_numpy()
    assert mask_enc.sum() == len(df) // stride and mask_enc[::stride].all(), path
    datos = df.loc[~mask_enc, ["x", "y", "z"]].to_numpy(dtype=np.float64)
    P = datos.reshape(-1, N_SITIOS, 3).mean(axis=1)
    np.save(cache, P)
    return P


def volumen_vs_T(corrida):
    """Devuelve f(T) -> volumen de la caja [m^3] (N_SITIOS celdas), interpolado de
    parametros-<CONC>-<corrida> (col. 0 = T, col. 7 = volumen de celda [A^3]).
    Fuera del rango de la tabla se usa el valor del extremo."""
    T, vcelda = np.loadtxt(DIR_PARAM / f"parametros-{CONC}-{corrida}", usecols=(0, 7), unpack=True)
    return lambda t: N_SITIOS * np.interp(t, T, vcelda) * A3_TO_M3


def eps_111_cov(P_macro_uc, T, V):
    P = P_macro_uc[DESCARTAR:] * U_P_UC_CM2_TO_C_M2
    Pm = P.mean(axis=0)
    cov = (P.T @ P) / len(P) - np.outer(Pm, Pm)
    eps = np.eye(3) + V / (EPS0 * KB * T) * cov
    return N111 @ eps @ N111


def smooth(y, n=N_SUAVIZADO):
    """Media movil centrada de n puntos; los bordes usan una ventana mas chica."""
    y = np.asarray(y, dtype=float)
    h = n // 2
    return np.array([y[max(0, i - h):i + h + 1].mean() for i in range(len(y))])


def temperaturas(corrida):
    out = []
    for p in (CARPETA_PLOCAL / f"Zr-{CONC}-{corrida}").glob("P-local-*"):
        sufijo = p.name.removeprefix("P-local-")
        if sufijo.isdigit() and int(sufijo) <= T_MAX:   # descarta p.ej. "590-parte2"
            out.append((int(sufijo), p))
    return sorted(out)


if __name__ == "__main__":
    res = {}
    with ProcessPoolExecutor(max_workers=2) as ex:
        for corrida in CORRIDAS:
            lista = temperaturas(corrida)
            Ps = list(ex.map(p_macro, [p for _, p in lista]))
            T = np.array([t for t, _ in lista], dtype=float)
            vol = volumen_vs_T(corrida)
            eps_cov = np.array([eps_111_cov(P, t, vol(t)) for P, (t, _) in zip(Ps, lista)])
            Tc, eps_campo = np.loadtxt(DATOS_CAMPO / f"permittivity-table-{CONC}-{corrida}.dat",
                                       usecols=(0, 1), unpack=True)
            res[corrida] = dict(T=T, cov=eps_cov, Tc=Tc, campo=eps_campo)

            salida = BASE_DIR / f"Permitividad-cov-{CONC}-{corrida}.dat"
            with open(salida, "w") as f:
                f.write("# T_sim(K)  eps_111_cov  eps_111_campo(NaN si no hay)\n")
                for t, e in zip(T, eps_cov):
                    k = np.where(Tc == t)[0]
                    ec = eps_campo[k[0]] if len(k) else np.nan
                    f.write(f"{t:6.0f}  {e:12.1f}  {ec:12.1f}\n")
            print(f"{corrida}: {len(lista)} temperaturas, {len(Ps[0])} pasos por archivo -> {salida.name}")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex="col")
    for j, corrida in enumerate(CORRIDAS):
        r = res[corrida]
        for i, (f, ylabel) in enumerate(((lambda y: y, r"$\varepsilon_{r,111}$"),
                                         (lambda y: 1 / y, r"$1/\varepsilon_{r,111}$"))):
            ax = axes[i, j]
            ax.plot(T_SCALE * r["T"], f(smooth(r["cov"])), "o-", color="black", ms=5, lw=1.8,
                    label=f"covarianza de P (suav. n={N_SUAVIZADO})")
            ax.plot(T_SCALE * r["Tc"], f(smooth(r["campo"])), "s-", color="tab:green", ms=5, lw=1.8,
                    label=r"campo eléctrico ($E\parallel$[111])" + f" (suav. n={N_SUAVIZADO})")
            ax.set_ylabel(ylabel if j == 0 else "", fontsize=18)
            ax.grid(True, linestyle="--", alpha=0.5)
            ax.minorticks_on()
            if i == 0:
                ax.set_title(f"x = {CONC}  {corrida}", fontsize=16, fontweight="bold")
                ax.legend(fontsize=12)
            else:
                ax.set_xlabel("Temperatura [K]", fontsize=14, fontweight="bold")
                ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(BASE_DIR / f"Permi-campo-vs-cov-{CONC}.png", dpi=200, bbox_inches="tight")
    plt.show()
