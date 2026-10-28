"""
Permitividad de BZT x=0.40 (S0, S50, S100) por fluctuaciones de la polarizacion
macroscopica (sin campo) y comparacion con el metodo de campo electrico (E || [111]):

    eps_ab(T) = delta_ab + V / (eps0 kB T) * ( <Pa Pb> - <Pa><Pb> )

Para comparar con el metodo de campo se usa la proyeccion sobre [111]:
    eps_111 = n^T eps n,  n = (1,1,1)/sqrt(3)

Los P-local se leen de P-local-temperatura/Zr-0.40-Sxx/P-local-<T> (T = temperatura
de simulacion, sin escalar). Se cachea P_macro(t) de cada archivo en cache_Pmacro/.
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AQUI = Path(__file__).parent
CARPETA_PLOCAL = Path("/home/matias/Escritorio/P-local-temperatura")
DATOS_CAMPO = AQUI.parent / "Datos_permitividad"
CACHE = AQUI / "cache_Pmacro"
CACHE.mkdir(exist_ok=True)

N_SITIOS = 21952
EPS0 = 8.854187817e-12      # F/m
KB = 1.380649e-23           # J/K
A_RED = 4.0e-10             # m
V = N_SITIOS * A_RED**3     # m^3
U_P_UC_CM2_TO_C_M2 = 1e-2   # uC/cm^2 -> C/m^2
T_FACTOR = 1.48             # reescalado de T usado en el resto de las figuras
DESCARTAR = 0               # pasos iniciales a excluir (equilibracion)
CORRIDAS = ["S0", "S50", "S100"]
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


def eps_cov(P_macro_uc, T):
    P = P_macro_uc[DESCARTAR:] * U_P_UC_CM2_TO_C_M2
    Pm = P.mean(axis=0)
    cov = (P.T @ P) / len(P) - np.outer(Pm, Pm)
    return np.eye(3) + V / (EPS0 * KB * T) * cov


N_SUAVIZADO = 5             # puntos de la media movil (impar)


def smooth(y, n=N_SUAVIZADO):
    """Media movil centrada de n puntos; los bordes usan una ventana mas chica."""
    y = np.asarray(y, dtype=float)
    h = n // 2
    return np.array([y[max(0, i - h):i + h + 1].mean() for i in range(len(y))])


def temperaturas(corrida):
    out = []
    for p in (CARPETA_PLOCAL / f"Zr-0.40-{corrida}").glob("P-local-*"):
        sufijo = p.name.removeprefix("P-local-")
        if sufijo.isdigit():          # descarta p.ej. "590-parte2"
            out.append((int(sufijo), p))
    return sorted(out)


if __name__ == "__main__":
    resultados = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        for corrida in CORRIDAS:
            lista = temperaturas(corrida)
            Ps = list(ex.map(p_macro, [p for _, p in lista]))
            T = np.array([t for t, _ in lista], dtype=float)
            eps = [eps_cov(P, t) for P, (t, _) in zip(Ps, lista)]
            resultados[corrida] = dict(
                T=T,
                eps_111=np.array([N111 @ e @ N111 for e in eps]),
                eps_prom=np.array([np.trace(e) / 3 for e in eps]),
                eps_diag=np.array([np.diag(e) for e in eps]),
                n_pasos=len(Ps[0]),
            )
            print(f"{corrida}: {len(lista)} temperaturas, {len(Ps[0])} pasos por archivo")

    # ── tabla de comparacion ──
    for corrida in CORRIDAS:
        r = resultados[corrida]
        Tc, ec, _ = np.loadtxt(DATOS_CAMPO / f"permittivity-e111-table-0.40-{corrida}.dat",
                               skiprows=1, unpack=True)
        r["T_campo"], r["eps_campo"] = Tc, ec
        print(f"\n{corrida}: T_sim  eps_111(cov)  eps_prom(cov)  eps_campo")
        for t, e111, ep in zip(r["T"], r["eps_111"], r["eps_prom"]):
            k = np.where(Tc == t)[0]
            print(f"  {t:5.0f}  {e111:10.1f}  {ep:10.1f}  {ec[k[0]]:8.1f}" if len(k)
                  else f"  {t:5.0f}  {e111:10.1f}  {ep:10.1f}        -")

    # ── figura: eps y 1/eps, cov vs campo, una columna por corrida ──
    colores = {"cov": "black", "campo": "tab:green"}
    fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharex=True)
    for j, corrida in enumerate(CORRIDAS):
        r = resultados[corrida]
        T = T_FACTOR * r["T"]
        Tc = T_FACTOR * r["T_campo"]
        for i, (f, ylabel) in enumerate(((lambda y: y, r"$\varepsilon_r$"),
                                         (lambda y: 1 / y, r"$1/\varepsilon_r$"))):
            ax = axes[i, j]
            ax.plot(T, f(smooth(r["eps_111"])), "o-", color=colores["cov"], ms=5, lw=2,
                    label=r"covarianza ($\varepsilon_{111}$)")
            ax.plot(Tc, f(smooth(r["eps_campo"])), "s-", color=colores["campo"], ms=5, lw=2,
                    label="campo eléctrico (E$\\parallel$111)")
            ax.set_ylabel(ylabel if j == 0 else "", size=22)
            ax.minorticks_on()
            ax.tick_params(axis="both", which="major", labelsize=14, width=2, length=7)
            ax.tick_params(axis="both", which="minor", width=1.5, length=3)
            for s in ax.spines.values():
                s.set_linewidth(1.5)
            if i == 0:
                ax.set_title(f"x = 0.40  {corrida}", size=18, fontweight="bold")
                ax.legend(fontsize=12)
            else:
                ax.set_xlabel(r"Temperature $[K]$", size=16)
                ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(AQUI / "Permitividad-covarianza-x040.png", dpi=200, bbox_inches="tight")
    print("\nFigura: Permitividad-covarianza-x040.png")
