import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator


# T-Burns: genera dos figuras (CONC=0.20 y 0.40) con 3x3 subplots cada una.
# Filas: 1) Curie fit (inverse permittivity + linear fit)
#        2) Volumen vs T (con ajuste lineal)
#        3) Tangentes / Polarización local vs T
# Columnas: S0, S50, S100 (usando los mismos colores que en los scripts originales)

COLOR_MAP = {'S0': '#1f77b4', 'S50': '#ff7f0e', 'S100': '#2ca02c'}
ADJUST_COLOR = '#4E4D4D'  # color gris para ajustes (como en tus scripts)

# Estilos y tamaños pensados para figuras de tesis (A4 con márgenes)
FIGSIZE = (24, 18)
ROW_TITLE_FONT = 28
LABEL_FONT = 24
TICK_FONT = 20
LEGEND_FONT = 18
SCALE_INV = 1e3  # escala para las inversas de la permitividad (reduce decimales en el plot)

# Controles de estilo (ajustables)
MARKER_SIZE = 6
LINE_WIDTH = 2
# Tick sizes (longitud de las barras en puntos)
MAJOR_TICK_LENGTH = 8
MINOR_TICK_LENGTH = 4
# Tick spacing (K)
LARGE_TICK = 200  # para los plots grandes
ZOOM_TICK = 100   # para los insets/zooms

def read_permittivity(conc, seg_label):
    # seg_label: 'S0'|'S50'|'S100' -> file name uses S-0, S-50, S-100
    seg_num = seg_label.replace('S','')
    fname = f'permittivity-e111-table-{conc}-S-{seg_num}.dat'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T = []
    e = []
    with open(fname, 'r') as f:
        for row in csv.reader(f, delimiter=' '):
            # filter out empty entries
            aux = [x for x in row if x]
            if not aux:
                continue
            # skip comment header lines starting with '#'
            if aux[0].startswith('#'):
                continue
            try:
                T.append(float(aux[0]))
                e.append(float(aux[1]))
            except Exception:
                continue
    return np.array(T) * 1.45, np.array(e)

def read_volume(conc, suffix):
    # suffix: 'sol'|'esfsol'|'esf'
    fname = f'../parametros-de-red-vs-T/parametros-{conc}-{suffix}'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T = []
    V = []
    with open(fname, 'r') as f:
        for row in csv.reader(f, delimiter=' '):
            aux = [x for x in row if x]
            if not aux:
                continue
            try:
                T.append(float(aux[0]))
                V.append(float(aux[7]))
            except Exception:
                continue
    return np.array(T) * 1.45, np.array(V)

def read_tangentes(conc, suffix):
    # suffix: 'sol'|'esfsol'|'esf'
    fname = f'../polarizaciones-rms-vs-T/p-local-datos/heating/pol-rms-locales-{conc}-{suffix}'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T = []
    p = []
    with open(fname, 'r') as f:
        for row in csv.reader(f, delimiter=' '):
            aux = [x for x in row if x]
            if not aux:
                continue
            try:
                T.append(float(aux[0]))
                # index 4 as in tu script tangentes-P-locales
                p.append(float(aux[4]))
            except Exception:
                continue
    return np.array(T) * 1.45, np.array(p)

def linear_func(x, a, b):
    return a * x + b

def make_figure_for_conc(conc):
    # Segregation mapping: columns
    cols = [('S0','sol'), ('S50','esfsol'), ('S100','esf')]

    fig, axes = plt.subplots(3, 3, figsize=FIGSIZE)

    # Row titles as requested (se colocarán por encima de cada fila al final)
    row_title_a = 'a) Ajuste de la Ley de Curie a la permitividad dieléctrica'
    row_title_b = 'b) Variación del volumen de la celda unidad con la temperatura'
    row_title_c = 'c) Variación de la polarización local de las celdas con la temperatura'

    # Storage for fit results to write to files
    curie_results = []
    volumen_results = []

    # Zoom ranges solicitados para conc 0.40 (valores en K)
    zoom_ranges_040 = {
        'S0': (200.0, 500.0),
        'S50': (270.0, 571.0),
        'S100': (400.0, 631.0)
    }

    # Row 1: Curie fits
    for j, (seg_label, suffix) in enumerate(cols):
        ax = axes[0, j]
        try:
            T_vals, e_vals = read_permittivity(conc, seg_label)
        except FileNotFoundError as e:
            ax.text(0.5,0.5,str(e), ha='center')
            continue

        # inverse
        inv_unscaled = 1.0 / e_vals
        inv_e = inv_unscaled * SCALE_INV

        # choose tail sizes consistent with Curie.py (last 6 for S0,S50; last 9 for S100)
        if seg_label == 'S100':
            nfit = 9
        else:
            nfit = 6

        if len(T_vals) < nfit:
            # fallback
            nfit = max(2, len(T_vals))

        T_fit = T_vals[-nfit:]
        inv_fit = inv_e[-nfit:]

        # fit on the SCALED inverse values, but convert back to unscaled for Curie constants
        try:
            popt, pcov = curve_fit(linear_func, T_fit, inv_fit, p0=[1.0, 1.0])
            slope_scaled, intercept_scaled = popt
        except Exception as e:
            slope_scaled, intercept_scaled = np.nan, np.nan

        # convert slope/intercept back to unscaled values for physical constants
        if not np.isnan(slope_scaled):
            slope = slope_scaled / SCALE_INV
        else:
            slope = np.nan
        if not np.isnan(intercept_scaled):
            intercept = intercept_scaled / SCALE_INV
        else:
            intercept = np.nan

        # Curie constants (computed using unscaled slope)
        curie_const = 1.0 / slope if slope and not np.isnan(slope) else np.nan
        T0 = -intercept * curie_const if not np.isnan(curie_const) else np.nan
        curie_results.append((seg_label, slope, intercept, curie_const, T0))

        # plot (use scaled values so decimals are nicer) — unir puntos con línea
        ax.plot(T_vals, inv_e, marker='o', linestyle='-', color=COLOR_MAP[seg_label], label=seg_label,
            markersize=MARKER_SIZE, linewidth=LINE_WIDTH)
        x = np.linspace(np.min(T_vals), np.max(T_vals), 200)
        # plot fit line using the SCALED fit parameters
        ax.plot(x, linear_func(x, slope_scaled, intercept_scaled), linestyle='--', color=ADJUST_COLOR,
            label='ajuste', linewidth=LINE_WIDTH)

        # Labels: sólo etiqueta Y en la primera columna, sólo X en la última fila
        if j == 0:
            ax.set_ylabel(r'$(\varepsilon_{r,111})^{-1}\ (\times 10^{-3})$', fontsize=LABEL_FONT)
        # no xlabel en filas superiores (solo en la tercera fila)

        ax.grid(True, linestyle='--', alpha=0.4)
        # fijar x-limits para que cada subplot use su propio rango de temperaturas
        ax.set_xlim(-1, T_vals[-1] + 10)  # un poco de margen a la derecha
        # asegurar ticks de temperatura cada LARGE_TICK K
        try:
            ax.xaxis.set_major_locator(MultipleLocator(LARGE_TICK))
        except Exception:
            pass
        try:
            ax.minorticks_on()
        except Exception:
            pass
        ax.tick_params(labelsize=TICK_FONT, labelbottom=False, labelleft=True, length=MAJOR_TICK_LENGTH)
        ax.tick_params(which='minor', length=MINOR_TICK_LENGTH)
        ax.legend(fontsize=LEGEND_FONT)
    # (Los títulos de fila se colocan más abajo, tras ajustar el layout)

    # Row 2: Volume
    for j, (seg_label, suffix) in enumerate(cols):
        ax = axes[1, j]
        try:
            T_vals, V_vals = read_volume(conc, suffix)
        except FileNotFoundError as e:
            ax.text(0.5,0.5,str(e), ha='center')
            continue

        # select tail sizes like volumen script: last6,last6,last9
        if seg_label == 'S100':
            nfit = 9
        else:
            nfit = 6
        if len(T_vals) < nfit:
            nfit = max(2, len(T_vals))

        T_fit = T_vals[-nfit:]
        V_fit = V_vals[-nfit:]

        try:
            popt, pcov = curve_fit(linear_func, T_fit, V_fit, p0=[1.0, 1.0])
            slope_v, intercept_v = popt
        except Exception:
            slope_v, intercept_v = np.nan, np.nan

        volumen_results.append((seg_label, slope_v, intercept_v))

        # plot and fit line (unir puntos con línea)
        ax.plot(T_vals, V_vals, marker='o', linestyle='-', color=COLOR_MAP[seg_label], label=seg_label
                , markersize=MARKER_SIZE, linewidth=LINE_WIDTH)
        x = np.linspace(np.min(T_vals), np.max(T_vals), 200)
        ax.plot(x, linear_func(x, slope_v, intercept_v), linestyle='--', color=ADJUST_COLOR, label='ajuste'
                , linewidth=LINE_WIDTH)
        if j == 0:
            ax.set_ylabel(r'Volumen [$\mathrm{\AA}^3$]', fontsize=LABEL_FONT)
        # no xlabel en esta fila (sólo en la tercera fila)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(-1, T_vals[-1] + 50)  # un poco de margen a la derecha
        # asegurar ticks de temperatura cada LARGE_TICK K
        try:
            ax.xaxis.set_major_locator(MultipleLocator(LARGE_TICK))
        except Exception:
            pass
        try:
            ax.minorticks_on()
        except Exception:
            pass
        ax.tick_params(labelsize=TICK_FONT, labelbottom=False, labelleft=True, length=MAJOR_TICK_LENGTH)
        ax.tick_params(which='minor', length=MINOR_TICK_LENGTH)
        ax.legend(fontsize=LEGEND_FONT)

        # --- Si estamos generando la figura para la concentración 0.40,
        # añadir un inset (zoom) en el mismo subplot en la parte inferior derecha
        if conc == "0.40":
            try:
                trange = zoom_ranges_040.get(seg_label, None)
                if trange is not None:
                    tmin, tmax = trange
                    # seleccionar los índices dentro del rango
                    mask = (T_vals >= tmin) & (T_vals <= tmax)
                    if np.any(mask):
                        V_zoom = V_vals[mask]
                        # calcular límites Y automáticos con margen
                        vmin = np.min(V_zoom)
                        vmax = np.max(V_zoom)
                        vrange = vmax - vmin if vmax > vmin else max(1e-6, abs(vmin)*0.01)
                        ypad = 0.06 * vrange
                        y0 = vmin - ypad
                        y1 = vmax + ypad

                        # crear inset dentro del mismo eje (lower right) pero un poco más arriba
                        # usamos bbox_to_anchor en coordenadas del eje para desplazar hacia arriba
                        axins = inset_axes(ax, width="100%", height="100%",
                                           bbox_to_anchor=(0.55, 0.12, 0.40, 0.35),
                                           bbox_transform=ax.transAxes)
                        # plot datos y ajuste en el inset con misma estética
                        axins.plot(T_vals, V_vals, marker='o', linestyle='-', color=COLOR_MAP[seg_label], 
                                   markersize=MARKER_SIZE, linewidth=LINE_WIDTH)
                        # plot ajuste en inset (usar parámetros slope_v/intercept_v)
                        xzoom = np.linspace(tmin, tmax, 100)
                        axins.plot(xzoom, linear_func(xzoom, slope_v, intercept_v), linestyle='--', color=ADJUST_COLOR,
                                    linewidth=LINE_WIDTH)
                        axins.set_xlim(tmin, tmax)
                        axins.set_ylim(y0, y1)
                        axins.minorticks_on()
                        axins.grid(True, linestyle='--', alpha=0.4)
                        axins.tick_params(labelsize=int(TICK_FONT*0.9), length=MAJOR_TICK_LENGTH)
                        axins.tick_params(which='minor', length=MINOR_TICK_LENGTH)
                        ax.grid(True, linestyle='--', alpha=0.4)
                        # mostrar escala de temperaturas (x labels) solo aquí
                        try:
                            axins.xaxis.set_major_locator(MultipleLocator(ZOOM_TICK))
                        except Exception:
                            pass
                        axins.xaxis.set_visible(True)
                        # no mostrar etiqueta Y textual en inset, solo números
                        # activar minor ticks y ajustar longitudes del eje padre también (seguridad)
                        try:
                            ax.minorticks_on()
                        except Exception:
                            pass
                        ax.tick_params(labelsize=TICK_FONT, labelbottom=False, labelleft=True, length=MAJOR_TICK_LENGTH)
                        ax.tick_params(which='minor', length=MINOR_TICK_LENGTH)

            except Exception:
                # si algo falla en el zoom, seguir sin interrumpir
                pass

    # (Row title colocado más arriba tras tight_layout)

    # Row 3: Tangentes / Polarización local
    for j, (seg_label, suffix) in enumerate(cols):
        ax = axes[2, j]
        try:
            T_vals, p_vals = read_tangentes(conc, suffix)
        except FileNotFoundError as e:
            ax.text(0.5,0.5,str(e), ha='center')
            continue

        # unir puntos con línea para visualizar curva
        ax.plot(T_vals, p_vals, marker='o', linestyle='-', color=COLOR_MAP[seg_label], label=seg_label,
                 markersize=MARKER_SIZE, linewidth=LINE_WIDTH)
        if j == 0:
            ax.set_ylabel(r'Polarización local [$\mu\mathrm{C}/\mathrm{cm}^2$]', fontsize=LABEL_FONT)
        # Sólo la tercera fila tendrá etiqueta X
        ax.set_xlabel('Temperatura [K]', fontsize=LABEL_FONT)
        ax.grid(True, linestyle='--', alpha=0.4)
        # # fijar x-limits para la fila 3 también
        ax.set_xlim(-1, T_vals[-1] + 10)  # un poco de margen a la derecha
        # asegurar ticks de temperatura cada LARGE_TICK K
        try:
            ax.xaxis.set_major_locator(MultipleLocator(LARGE_TICK))
        except Exception:
            pass
        # para fila 3 dejamos visible el labelbottom
        try:
            ax.minorticks_on()
        except Exception:
            pass
        ax.tick_params(labelsize=TICK_FONT, labelbottom=True, labelleft=True, length=MAJOR_TICK_LENGTH)
        ax.tick_params(which='minor', length=MINOR_TICK_LENGTH)
        ax.legend(fontsize=LEGEND_FONT)

    # Ajustar separación entre filas: queremos exactamente 1 cm entre filas
    CM_TO_INCH = 0.393700787
    desired_gap_inch = 3.0 * CM_TO_INCH
    # hspace es la fracción del alto medio del subplot; aproximamos:
    hspace = (desired_gap_inch * 3) / FIGSIZE[1]
    # Ajustar subplots: reducir márgenes laterales y separación horizontal
    fig.subplots_adjust(top=0.96, hspace=hspace, left=0.06, right=0.98, wspace=0.12)

    # Place the row titles centered in the added gap above each row
    try:
        for row_idx, title in [(0, row_title_a), (1, row_title_b), (2, row_title_c)]:
            bbox = axes[row_idx, 0].get_position()
            # compute offset equal to half the desired gap in figure coords
            offset = (desired_gap_inch / FIGSIZE[1]) / 4.0
            y = min(bbox.y1 + offset, 0.98)
            fig.text(0.5, y, title, ha='center', fontsize=ROW_TITLE_FONT)
    except Exception:
        # fallback: colocar títulos en posiciones fijas
        fig.text(0.5, 0.95, row_title_a, ha='center', fontsize=ROW_TITLE_FONT)
        fig.text(0.5, 0.62, row_title_b, ha='center', fontsize=ROW_TITLE_FONT)
        fig.text(0.5, 0.29, row_title_c, ha='center', fontsize=ROW_TITLE_FONT)

    # Save figure
    out_png = f'T_Burns-{conc}.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f'Saved figure: {out_png}')

    # Save fit results to text files (Curie and Volumen)
    out_curie = f'T_Burns-{conc}-Curie.txt'
    with open(out_curie, 'w') as f:
        f.write(f'Concentracion: {conc}\n')
        f.write('Curie fits (seg, slope, intercept, Curie_const, T0)\n')
        for r in curie_results:
            f.write(f'{r[0]}\t{slope_format(r[1])}\t{slope_format(r[2])}\t{slope_format(r[3])}\t{slope_format(r[4])}\n')
    print(f'Saved Curie results: {out_curie}')

    out_vol = f'T_Burns-{conc}-Volumen.txt'
    with open(out_vol, 'w') as f:
        f.write(f'Concentracion: {conc}\n')
        f.write('Volumen linear fits (seg, slope, intercept)\n')
        for r in volumen_results:
            f.write(f'{r[0]}\t{slope_format(r[1])}\t{slope_format(r[2])}\n')
    print(f'Saved Volumen results: {out_vol}')

def slope_format(x):
    try:
        return f'{float(x):.6e}'
    except Exception:
        return 'nan'

def main():
    for conc in ["0.20", "0.40"]:
        try:
            make_figure_for_conc(conc)
        except Exception as e:
            print(f'Error generating figure for {conc}: {e}')

if __name__ == '__main__':
    main()
