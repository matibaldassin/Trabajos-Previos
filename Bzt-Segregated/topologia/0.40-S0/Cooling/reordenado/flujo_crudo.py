#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uso
---
    python3 flujo_crudo.py Pol-10.vtk --sigma 1.5 --surrogates 12 

Solo requiere numpy y matplotlib.
"""

import argparse
import os

import numpy as np

PLANES = {"xy": (0, 1, 2), "xz": (0, 2, 1), "yz": (1, 2, 0)}  # (eje_a, eje_b, normal)
RADII = (1, 2, 3)  # -> ventanas 3x3, 5x5, 7x7
SCALE = 1e8  # se grafica en unidades de 10^8 para evitar el offset "1e8" de matplotlib

# (a,b) esta fijado por el orden natural de ejes que deja np al hacer el
# corte (siempre ejes crecientes, necesario para que las flechas P se
# grafiquen en el eje correcto), pero para "xz" ese orden (x,z) es el
# opuesto del par ciclico que le corresponde a normal=y (que es (z,x)).
# Sin este signo, positivo significa "antihorario" en los graficos de xy
# e yz pero "horario" en xz: mismo vortice, colores con sentido invertido.
_CYCLIC = {(0, 1, 2), (1, 2, 0), (2, 0, 1)}
PLANE_SIGN = {pl: (1.0 if (n, a, b) in _CYCLIC else -1.0)
              for pl, (a, b, n) in PLANES.items()}


# ----------------------------------------------------------------------
# lectura
# ----------------------------------------------------------------------
def read_legacy_vtk(path):
    """Lector robusto de legacy VTK ASCII.

    NO asume que los puntos vengan en orden de barrido: construye el mapeo
    indice->celda a partir de las coordenadas reales de cada punto. Un
    reshape ciego sobre un fichero desordenado revuelve el campo
    espacialmente sin dar ningun aviso, que es el error mas insidioso
    posible en este analisis.
    """
    L = open(path).read().splitlines()
    dims = None
    ip = None
    for i, l in enumerate(L[:64]):
        u = l.strip().upper()
        if u.startswith("DIMENSIONS"):
            dims = tuple(int(v) for v in l.split()[1:4])
        elif u.startswith("POINTS"):
            ip = i + 1
            break
    if dims is None or ip is None:
        raise ValueError("Cabecera incompleta en %s" % path)
    n = int(np.prod(dims))

    def grab(s, c):
        v = []
        i = s
        while len(v) < c:
            v.extend(float(x) for x in L[i].split())
            i += 1
        return np.array(v[:c]), i

    pts, cur = grab(ip, 3 * n)
    pts = pts.reshape(n, 3)
    axes = [np.unique(pts[:, c]) for c in range(3)]
    if tuple(len(a) for a in axes) != dims:
        raise ValueError("valores unicos %s != DIMENSIONS %s"
                         % (tuple(len(a) for a in axes), dims))
    idx = tuple(np.searchsorted(axes[c], pts[:, c]) for c in range(3))

    lin = (idx[0] * dims[1] + idx[1]) * dims[2] + idx[2]
    if len(np.unique(lin)) != n:
        raise ValueError("Los puntos no forman una malla completa sin repetidos.")
    if not np.array_equal(lin, np.arange(n)):
        print("  AVISO: los puntos NO vienen en orden de barrido; se reordenan "
              "usando las coordenadas.")

    scal, vec = {}, {}
    i = cur
    while i < len(L):
        u = L[i].strip().upper()
        if u.startswith("SCALARS"):
            nm = L[i].split()[1]
            j = i + 1 + (1 if L[i + 1].strip().upper().startswith("LOOKUP") else 0)
            a, i = grab(j, n)
            A = np.empty(dims)
            A[idx] = a
            scal[nm] = A
            continue
        if u.startswith("VECTORS"):
            nm = L[i].split()[1]
            a, i = grab(i + 1, 3 * n)
            a = a.reshape(n, 3)
            A = np.empty(dims + (3,))
            for c in range(3):
                A[idx + (c,)] = a[:, c]
            vec[nm] = A
            continue
        i += 1
    return dims, pts, scal, vec


# ----------------------------------------------------------------------
# operadores
# ----------------------------------------------------------------------
DX = 4.0e-8  # espaciado fisico de la grilla (celda), igual que en 01_parse_and_compute.py


def curl(F, dx=DX):
    """Rotor con diferencias centradas y PBC. F:(nx,ny,nz,3) -> (nx,ny,nz,3)"""
    d = lambda A, ax: (np.roll(A, -1, axis=ax) - np.roll(A, 1, axis=ax)) / (2 * dx)
    Fx, Fy, Fz = F[..., 0], F[..., 1], F[..., 2]
    return np.stack(
        [d(Fz, 1) - d(Fy, 2), d(Fx, 2) - d(Fz, 0), d(Fy, 0) - d(Fx, 1)], axis=-1
    )


def boxsum(A, r, axa, axb):
    """Suma sobre el parche (2r+1)x(2r+1) con PBC. Separable: O(4r) rolls."""
    out = A
    for ax in (axa, axb):
        acc = out.copy()
        for s in range(1, r + 1):
            acc = acc + np.roll(out, s, axis=ax) + np.roll(out, -s, axis=ax)
        out = acc
    return out


def boxmean(A, r, axa, axb):
    """Media aritmetica sobre el parche (2r+1)x(2r+1)."""
    return boxsum(A, r, axa, axb) / (2 * r + 1) ** 2


def gaussian_smooth(F, sigma):
    """Suavizado gaussiano periodico por FFT, componente a componente."""
    if sigma <= 0:
        return F.copy()
    ks = [np.fft.fftfreq(n) for n in F.shape[:3]]
    G = np.exp(-2 * np.pi ** 2 * sigma ** 2 *
               (ks[0][:, None, None] ** 2 + ks[1][None, :, None] ** 2
                + ks[2][None, None, :] ** 2))
    return np.stack([np.real(np.fft.ifftn(np.fft.fftn(F[..., c]) * G))
                     for c in range(3)], axis=-1)


def phase_randomize(F, rng):
    """Modelo nulo: conserva |FFT| y correlaciones cruzadas, aleatoriza fases."""
    W = np.fft.fftn(rng.normal(size=F.shape[:3]))
    m = np.abs(W)
    ph = W / np.where(m < 1e-30, 1e-30, m)
    return np.stack(
        [np.real(np.fft.ifftn(np.fft.fftn(F[..., c]) * ph)) for c in range(3)], axis=-1
    )


# ----------------------------------------------------------------------
# calculo
# ----------------------------------------------------------------------
def compute_flux(P, sigma=0.0):
    """Flujo del rotor multiescala en las tres familias de planos.

    A diferencia de flujo.py, el rotor se calcula sobre el campo de
    polarizacion CRUDO (P, opcionalmente suavizado con --sigma), sin
    normalizar por |P|. Es decir C = rot(P), no rot(P/|P|). El parche se
    promedia sin ponderar (ver docstring del modulo).

    sigma : si > 0, suaviza P antes de derivar (reduce el ruido de alta
            frecuencia del rotor sin tocar su normalizacion, ya que aqui
            no hay paso de normalizado).
    """
    mod_raw = np.linalg.norm(P, axis=-1)
    Pw = gaussian_smooth(P, sigma) if sigma > 0 else P
    C = curl(Pw)

    out = {"_curl": C, "_mod": mod_raw, "_Pw": Pw, "_P": P}
    for pl, (a, b, nrm) in PLANES.items():
        for r in RADII:
            ww = 2 * r + 1
            out["flujo_%s_%d" % (pl, ww)] = (
                PLANE_SIGN[pl] * boxmean(C[..., nrm], r, a, b))
    return out


# ----------------------------------------------------------------------
# principal
# ----------------------------------------------------------------------
def run(args):
    dims, coords, scal, vec = read_legacy_vtk(args.vtk)
    name = args.vector or sorted(vec)[0]
    P = vec[name]
    print("Malla %s | campo '%s' | escalares: %s" % (dims, name, list(scal)))
    print("sigma=%.2g  (rotor sobre P crudo, sin normalizar, parche sin ponderar)"
          % args.sigma)

    F = compute_flux(P, sigma=args.sigma)
    os.makedirs(args.out, exist_ok=True)

    null = {}
    if args.surrogates:
        rng = np.random.default_rng(2024)
        acc = {}
        for _ in range(args.surrogates):
            Fs = compute_flux(phase_randomize(P, rng), sigma=args.sigma)
            for pl in PLANES:
                for r in RADII:
                    k = "flujo_%s_%d" % (pl, 2 * r + 1)
                    acc.setdefault(k, []).append(Fs[k].std())
        null = {k: (np.mean(v), np.std(v)) for k, v in acc.items()}

    Ws = [2 * r + 1 for r in RADII]
    print("\n%-5s | %-24s | %-9s | %s"
          % ("plano", "rms  3x3 / 5x5 / 7x7", "pendiente", "z vs nulo (3,5,7)"))
    print("-" * 78)
    for pl in PLANES:
        v = [F["flujo_%s_%d" % (pl, w)].std() for w in Ws]
        slope = np.polyfit(np.log(Ws), np.log(v), 1)[0]
        zs = ""
        if null:
            z = []
            for w in Ws:
                m, s = null["flujo_%s_%d" % (pl, w)]
                z.append((F["flujo_%s_%d" % (pl, w)].std() - m) / s if s > 0 else 0.0)
            zs = "%+6.1f %+6.1f %+6.1f" % tuple(z)
        print("%-5s | %7.4f %7.4f %7.4f  | %+6.2f   | %s"
              % (pl, v[0], v[1], v[2], slope, zs))
    print("\nreferencia de pendiente: ruido puro ~ -1.5 ; estructura real ~ 0")

    np.savez_compressed(
        os.path.join(args.out, "flujo_crudo.npz"),
        **{k: v for k, v in F.items() if not k.startswith("_")},
        curl=F["_curl"], modulo=F["_mod"],
    )
    for ww in (2 * r + 1 for r in RADII):
        make_cortes(F, scal, args, ww)
    if args.vtk_out:
        write_vtk(os.path.join(args.out, "flujo_crudo.vtk"), dims, coords, F, P)
    print("\nSalidas en %s" % os.path.abspath(args.out))


def _axis_triad(ax, lab_h, lab_v, lab_n, x0=-0.05, y0=-0.05, length=0.16):
    """Triedro de orientacion, dibujado por FUERA del cuadro de datos (esquina
    inferior izquierda): flechas para los dos ejes del plano (horizontal y
    vertical en pantalla) y un circulo con punto central para el eje que sale
    hacia el/la observador/a (convencion ⊙ = saliente, ⊗ = entrante)."""
    arrow_kw = dict(xycoords="axes fraction", textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="-|>", color="black", lw=1.8,
                                     shrinkA=0, shrinkB=0),
                     annotation_clip=False, zorder=10)
    ax.annotate("", xy=(x0 + length, y0), xytext=(x0, y0), **arrow_kw)
    ax.annotate("", xy=(x0, y0 + length), xytext=(x0, y0), **arrow_kw)
    ax.annotate(lab_h, xy=(x0 + length + 0.025, y0), xycoords="axes fraction",
                ha="left", va="center", fontsize=11, fontweight="bold",
                annotation_clip=False, zorder=10)
    ax.annotate(lab_v, xy=(x0, y0 + length + 0.025), xycoords="axes fraction",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                annotation_clip=False, zorder=10)
    ax.plot(x0, y0, marker="o", markersize=11, markerfacecolor="white",
            markeredgecolor="black", markeredgewidth=1.4,
            transform=ax.transAxes, clip_on=False, zorder=10)
    ax.plot(x0, y0, marker=".", markersize=4, color="black",
            transform=ax.transAxes, clip_on=False, zorder=11)
    ax.annotate(lab_n, xy=(x0 - 0.03, y0 - 0.03), xycoords="axes fraction",
                ha="right", va="top", fontsize=11, fontweight="bold",
                annotation_clip=False, zorder=10)


# ----------------------------------------------------------------------
# figura: cortes
# ----------------------------------------------------------------------
def make_cortes(F, scal, args, ww=5):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Praw = F["_P"]
    idn = scal.get("identidad")
    dims = Praw.shape[:3]
    ncor = args.cortes

    # quiver: 'scale' es INVERSO a la longitud de la flecha.
    scale_all = args.arrow_scale * np.linalg.norm(Praw, axis=-1).max()

    # escala de color COMUN a los 3 planos, para poder comparar intensidad
    # entre filas (si cada subplot se autoescala a su propio maximo, un
    # vortice 10x mas debil se ve igual de "intenso" que uno fuerte).
    # --vmax permite fijarla a mano (p.ej. para comparar contra otra figura).
    vmax = args.vmax if args.vmax is not None else max(
        np.abs(F["flujo_%s_%d" % (pl, ww)]).max() for pl in PLANES)

    fig, ax = plt.subplots(3, ncor, figsize=(4.2 * ncor, 5.0 * 3), squeeze=False)
    im = None
    for i, pl in enumerate(PLANES):
        a, b, nrm = PLANES[pl]
        n_normal = dims[nrm]
        ks = [k * n_normal // ncor for k in range(ncor)]
        for j, k in enumerate(ks):
            sl = [slice(None), slice(None), slice(None)]
            sl[nrm] = k
            sl = tuple(sl)
            axij = ax[i][j]

            gm = F["flujo_%s_%d" % (pl, ww)][sl] / SCALE
            im = axij.imshow(gm.T, origin="lower", cmap="RdBu_r",
                             vmin=-vmax / SCALE, vmax=vmax / SCALE)

            V = Praw[sl]
            n0, n1 = V.shape[0], V.shape[1]
            xx, yy = np.meshgrid(np.arange(n0), np.arange(n1), indexing="ij")
            if idn is not None:
                m = idn[sl]
                zr, ti = np.isclose(m, 1), ~np.isclose(m, 1)
                axij.scatter(xx[zr], yy[zr], s=46, facecolors="none",
                             edgecolors="0.45", linewidths=0.6, zorder=2)
                axij.quiver(xx[ti], yy[ti], V[..., a][ti], V[..., b][ti],
                            pivot="mid", scale=scale_all, width=0.004,
                            color="k", zorder=3, label="Ti")
                axij.quiver(xx[zr], yy[zr], V[..., a][zr], V[..., b][zr],
                            pivot="mid", scale=scale_all / args.zr_boost,
                            width=0.004, color="darkgreen", zorder=4, label="Zr")
                if i == 0 and j == 0:
                    axij.legend(loc="upper right", fontsize=14, framealpha=0.9)
            else:
                axij.quiver(xx, yy, V[..., a], V[..., b], pivot="mid",
                            scale=scale_all, width=0.004, color="k")
            axij.set_xticks([])
            axij.set_yticks([])
            if j == 0:
                lab_n = ("" if PLANE_SIGN[pl] > 0 else "-") + "xyz"[nrm]
                _axis_triad(axij, "xyz"[a], "xyz"[b], lab_n)
            lab = "%s = %d" % ("xyz"[nrm], k)
            axij.set_title("%s plane, %s" % (pl, lab), fontsize=14, fontweight="bold")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle("S0 - T = %s" % args.temp, fontsize=20, y=0.98,fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02,
                        label=r"$\nabla\times\mathbf{P}\cdot\hat{n}$ $(\cdot 10^8~\mu C/cm^3)$")
    cbar.ax.ticklabel_format(style="plain", useOffset=False)
    out = os.path.join(args.out, "cortes_crudo_%dx%d.png" % (ww, ww))
    fig.savefig(out, dpi=600, bbox_inches="tight")
    print("Figura: %s" % out)


# ----------------------------------------------------------------------
# export a ParaView
# ----------------------------------------------------------------------
def write_vtk(path, dims, coords, F, P):
    o = coords.min(axis=0)
    sp = np.array([np.diff(np.unique(coords[:, c]))[0]
                   if len(np.unique(coords[:, c])) > 1 else 1.0 for c in range(3)])
    n = int(np.prod(dims))
    flat = lambda A: np.transpose(A, (2, 1, 0)).ravel()

    def dump(fh, arr, per=6):
        a = np.asarray(arr).ravel()
        k = a.size - a.size % per
        if k:
            np.savetxt(fh, a[:k].reshape(-1, per), fmt="%.5g", delimiter=" ")
        if a.size > k:
            fh.write(" ".join("%.5g" % v for v in a[k:]) + "\n")

    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nflujo del rotor (P crudo)\nASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS %d %d %d\n" % dims)
        f.write("ORIGIN %g %g %g\nSPACING %g %g %g\n" % (*o, *sp))
        f.write("POINT_DATA %d\n" % n)
        for k in sorted(k for k in F if not k.startswith("_")):
            f.write("SCALARS %s float\nLOOKUP_TABLE default\n" % k)
            dump(f, flat(F[k]))
        f.write("SCALARS modulo float\nLOOKUP_TABLE default\n")
        dump(f, flat(F["_mod"]))
        for nm, A in (("P_suavizado", F["_Pw"]), ("rotor", F["_curl"]), ("P", P)):
            f.write("VECTORS %s float\n" % nm)
            np.savetxt(f, np.stack([flat(A[..., c]) for c in range(3)], axis=1),
                       fmt="%.5g", delimiter=" ")
    print("VTK para ParaView: %s" % path)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vtk", help="fichero .vtk de entrada")
    p.add_argument("--vector", default=None, help="nombre del campo VECTORS")
    p.add_argument("--sigma", type=float, default=0.0,
                   help="suavizado de P antes de derivar el rotor (celdas)")
    p.add_argument("--vmax", type=float, default=None,
                   help="escala de color fija: rango [-vmax, vmax]; si se "
                        "omite, autoescala al maximo de cada archivo (por "
                        "tamano de ventana)")
    p.add_argument("--arrow-scale", type=float, default=16.0,
                   help="mayor = flechas mas cortas")
    p.add_argument("--zr-boost", type=float, default=3.0,
                   help="alarga las flechas de Zr por este factor (Ti no "
                        "cambia); |P|_Zr ~ 4-5x mas chico que |P|_Ti, asi "
                        "que sin esto casi no se ven")
    p.add_argument("--cortes", type=int, default=4,
                   help="cantidad de cortes por familia de planos, "
                        "equiespaciados barriendo la supercelda")
    p.add_argument("--surrogates", type=int, default=0, metavar="N",
                   help="realizaciones del modelo nulo")
    p.add_argument("--vtk-out", action="store_true", help="exportar .vtk")
    p.add_argument("--temp", required=True,
                   help="temperatura a mostrar en el titulo comun de las figuras "
                        "('Temperature = <temp>'), p.ej. --temp '300 K'")
    p.add_argument("--out", default="resultados_flujo_crudo")
    run(p.parse_args())
