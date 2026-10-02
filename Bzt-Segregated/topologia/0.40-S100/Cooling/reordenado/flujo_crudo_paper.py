#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage
-----
    python3 flujo_crudo_paper.py Pol-10.vtk --window 5 --vmax 6 --sigma 1.5
"""

import argparse
import os

import numpy as np

import flujo_crudo as fc

CORTES_S0 = [("xy", 0), ("xy", 14), ("xz", 14)]
CORTES_S100 = [("xy", 21), ("xz", 14), ("xz", 21)]
SCALE = 1e8  # se grafica en unidades de 10^8 para evitar el offset "1e8" de matplotlib


def load_dataset(path, vector, sigma):
    dims, coords, scal, vec = fc.read_legacy_vtk(path)
    name = vector or sorted(vec)[0]
    P = vec[name]
    print("Grid %s | field '%s' | file %s" % (dims, name, path))
    return fc.compute_flux(P, sigma=sigma), scal


def make_figure(panels, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ww = args.window
    vmax = args.vmax
    ncols = len(panels) // 2

    fig, ax = plt.subplots(2, ncols, figsize=(4.6 * ncols, 4.6 * 2),
                           squeeze=False, gridspec_kw={"hspace": 0.04})
    ax = ax.ravel()
    im = None
    for j, (F, scal, pl, k) in enumerate(panels):
        Praw = F["_P"]
        idn = scal.get("identidad")
        # quiver: 'scale' is INVERSE to arrow length.
        scale_all = args.arrow_scale * np.linalg.norm(Praw, axis=-1).max()

        a, b, nrm = fc.PLANES[pl]
        sl = [slice(None), slice(None), slice(None)]
        sl[nrm] = k
        sl = tuple(sl)
        axj = ax[j]

        gm = F["flujo_%s_%d" % (pl, ww)][sl] / SCALE
        im = axj.imshow(gm.T, origin="lower", cmap="RdBu_r",
                        vmin=-vmax / SCALE, vmax=vmax / SCALE)

        V = Praw[sl]
        n0, n1 = V.shape[0], V.shape[1]
        xx, yy = np.meshgrid(np.arange(n0), np.arange(n1), indexing="ij")
        if idn is not None:
            m = idn[sl]
            zr, ti = np.isclose(m, 1), ~np.isclose(m, 1)
            axj.scatter(xx[zr], yy[zr], s=46, facecolors="none",
                        edgecolors="0.45", linewidths=0.6, zorder=2)
            axj.quiver(xx[ti], yy[ti], V[..., a][ti], V[..., b][ti],
                       pivot="mid", scale=scale_all, width=0.004,
                       color="k", zorder=3, label="Ti")
            axj.quiver(xx[zr], yy[zr], V[..., a][zr], V[..., b][zr],
                       pivot="mid", scale=scale_all / args.zr_boost,
                       width=0.004, color="darkgreen", zorder=4, label="Zr")
            if j == 0:
                axj.legend(loc="upper right", fontsize=8, framealpha=0.9)
        else:
            axj.quiver(xx, yy, V[..., a], V[..., b], pivot="mid",
                       scale=scale_all, width=0.004, color="k")
        axj.set_xticks([])
        axj.set_yticks([])

    fig.subplots_adjust(top=0.92)

    cbar = fig.colorbar(im, ax=list(ax), fraction=0.03, pad=0.02,
                label=r"$(\nabla \times \mathbf{P}) \cdot \hat{n}$ ($\cdot10^8~\mu$C/cm$^3$)")
    cbar.ax.ticklabel_format(style="plain", useOffset=False)

    fig.canvas.draw()
    row0, row1 = ax[:ncols], ax[ncols:]
    x_center = (min(a.get_position().x0 for a in row0)
                + max(a.get_position().x1 for a in row0)) / 2
    top_y1 = max(a.get_position().y1 for a in row0)
    gap_center = (min(a.get_position().y0 for a in row0)
                  + max(a.get_position().y1 for a in row1)) / 2
    fig.text(x_center, top_y1 + 0.02, "a) S0", ha="center", va="bottom",
              fontsize=18, fontweight="bold")
    fig.text(x_center, gap_center, "b) S100", ha="center", va="center",
              fontsize=18, fontweight="bold")

    os.makedirs(args.out, exist_ok=True)
    out = os.path.join(args.out, "paper_%dx%d.png" % (ww, ww))
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print("Figure: %s" % out)


def run(args):
    F_s0, scal_s0 = load_dataset(args.other_vtk, args.vector, args.sigma)
    F_s100, scal_s100 = load_dataset(args.vtk, args.vector, args.sigma)
    print("window=%dx%d  sigma=%.2g  vmax=%.3g (fixed)"
          % (args.window, args.window, args.sigma, args.vmax))

    panels = ([(F_s0, scal_s0, pl, k) for pl, k in CORTES_S0]
              + [(F_s100, scal_s100, pl, k) for pl, k in CORTES_S100])
    make_figure(panels, args)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vtk", help="S100 input .vtk file (plotted last, panels d-f)")
    p.add_argument("--other-vtk",
                   default="../../../0.40-S0/Cooling/reordenado/Pol-10.vtk",
                   help="S0 input .vtk file (plotted first, panels a-c)")
    p.add_argument("--vector", default=None, help="name of the VECTORS field")
    p.add_argument("--sigma", type=float, default=0.0,
                   help="smoothing of P before differentiating the curl (cells)")
    p.add_argument("--window", type=int, default=5, choices=(3, 5, 7),
                   help="flux window size to plot")
    p.add_argument("--vmax", type=float, required=True,
                   help="fixed color scale: range [-vmax, vmax]; pick it "
                        "once and reuse it so both datasets are comparable")
    p.add_argument("--arrow-scale", type=float, default=16.0,
                   help="larger = shorter arrows")
    p.add_argument("--zr-boost", type=float, default=3.0,
                   help="lengthens Zr arrows by this factor (Ti unchanged); "
                        "|P|_Zr is ~4-5x smaller than |P|_Ti, so without "
                        "this they're barely visible")
    p.add_argument("--out", default="resultados_flujo_crudo")
    run(p.parse_args())
