import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------------------------------------------------
# READING FUNCTIONS

def read_permittivity(conc, seg_label):
    fname = f'Datos_Permitividades/permittivity-e111-table-{conc}-{seg_label}.dat'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T, e = [], []
    with open(fname, 'r') as f:
        for row in csv.reader(f, delimiter=' '):
            aux = [x for x in row if x]
            if not aux or aux[0].startswith('#'):
                continue
            try:
                T.append(float(aux[0]))
                e.append(float(aux[1]))
            except Exception:
                continue
    return np.array(T) * 1.45, np.array(e)

def read_volume(conc, seg_label):
    fname = f'Datos_Parámetros_de_Red/parametros-{conc}-{seg_label}'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T, V = [], []
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

def read_tangentes(conc, seg_label):
    fname = f'Datos_P_local/P-local-{conc}-{seg_label}'
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    T, p = [], []
    with open(fname, 'r') as f:
        for row in csv.reader(f, delimiter=' '):
            aux = [x for x in row if x]
            if not aux:
                continue
            try:
                T.append(float(aux[0]))
                p.append(float(aux[4]))
            except Exception:
                continue
    return np.array(T) * 1.45, np.array(p)

def linear_func(x, a, b):
    return a * x + b

# -----------------------------------------------------------------------
# READING DATA

T_curie_020_S0,   e_020_S0   = read_permittivity("0.20", "S0")
T_curie_020_S50,  e_020_S50  = read_permittivity("0.20", "S50")
T_curie_020_S100, e_020_S100 = read_permittivity("0.20", "S100")

T_curie_040_S0,   e_040_S0   = read_permittivity("0.40", "S0")
T_curie_040_S50,  e_040_S50  = read_permittivity("0.40", "S50")
T_curie_040_S100, e_040_S100 = read_permittivity("0.40", "S100")

T_vol_020_S0,   V_020_S0   = read_volume("0.20", "S0")
T_vol_020_S50,  V_020_S50  = read_volume("0.20", "S50")
T_vol_020_S100, V_020_S100 = read_volume("0.20", "S100")

T_vol_040_S0,   V_040_S0   = read_volume("0.40", "S0")
T_vol_040_S50,  V_040_S50  = read_volume("0.40", "S50")
T_vol_040_S100, V_040_S100 = read_volume("0.40", "S100")

T_pol_020_S0,   p_020_S0   = read_tangentes("0.20", "S0")
T_pol_020_S50,  p_020_S50  = read_tangentes("0.20", "S50")
T_pol_020_S100, p_020_S100 = read_tangentes("0.20", "S100")

T_pol_040_S0,   p_040_S0   = read_tangentes("0.40", "S0")
T_pol_040_S50,  p_040_S50  = read_tangentes("0.40", "S50")
T_pol_040_S100, p_040_S100 = read_tangentes("0.40", "S100")

# -----------------------------------------------------------------------
# DATA PROCESSING

SCALE_INV = 1e3

inv_e_020_S0   = (1.0 / e_020_S0)   * SCALE_INV
inv_e_020_S50  = (1.0 / e_020_S50)  * SCALE_INV
inv_e_020_S100 = (1.0 / e_020_S100) * SCALE_INV

inv_e_040_S0   = (1.0 / e_040_S0)   * SCALE_INV
inv_e_040_S50  = (1.0 / e_040_S50)  * SCALE_INV
inv_e_040_S100 = (1.0 / e_040_S100) * SCALE_INV

nfit_S0 = 6
nfit_S50 = 6
nfit_S100 = 9

# Curie fits (scaled, then converted back)
popt, _ = curve_fit(linear_func, T_curie_020_S0[-nfit_S0:],     inv_e_020_S0[-nfit_S0:],     p0=[1., 1.])
s_curie_020_S0_sc,   i_curie_020_S0_sc   = popt
s_curie_020_S0,   i_curie_020_S0   = s_curie_020_S0_sc / SCALE_INV,   i_curie_020_S0_sc / SCALE_INV
C_020_S0   = 1.0 / s_curie_020_S0
T0_020_S0  = -i_curie_020_S0 * C_020_S0

popt, _ = curve_fit(linear_func, T_curie_020_S50[-nfit_S50:],   inv_e_020_S50[-nfit_S50:],   p0=[1., 1.])
s_curie_020_S50_sc,  i_curie_020_S50_sc  = popt
s_curie_020_S50,  i_curie_020_S50  = s_curie_020_S50_sc / SCALE_INV,  i_curie_020_S50_sc / SCALE_INV
C_020_S50  = 1.0 / s_curie_020_S50
T0_020_S50 = -i_curie_020_S50 * C_020_S50

popt, _ = curve_fit(linear_func, T_curie_020_S100[-nfit_S100:], inv_e_020_S100[-nfit_S100:], p0=[1., 1.])
s_curie_020_S100_sc, i_curie_020_S100_sc = popt
s_curie_020_S100, i_curie_020_S100 = s_curie_020_S100_sc / SCALE_INV, i_curie_020_S100_sc / SCALE_INV
C_020_S100 = 1.0 / s_curie_020_S100
T0_020_S100 = -i_curie_020_S100 * C_020_S100

popt, _ = curve_fit(linear_func, T_curie_040_S0[-nfit_S0:],     inv_e_040_S0[-nfit_S0:],     p0=[1., 1.])
s_curie_040_S0_sc,   i_curie_040_S0_sc   = popt
s_curie_040_S0,   i_curie_040_S0   = s_curie_040_S0_sc / SCALE_INV,   i_curie_040_S0_sc / SCALE_INV
C_040_S0   = 1.0 / s_curie_040_S0
T0_040_S0  = -i_curie_040_S0 * C_040_S0

popt, _ = curve_fit(linear_func, T_curie_040_S50[-nfit_S50:],   inv_e_040_S50[-nfit_S50:],   p0=[1., 1.])
s_curie_040_S50_sc,  i_curie_040_S50_sc  = popt
s_curie_040_S50,  i_curie_040_S50  = s_curie_040_S50_sc / SCALE_INV,  i_curie_040_S50_sc / SCALE_INV
C_040_S50  = 1.0 / s_curie_040_S50
T0_040_S50 = -i_curie_040_S50 * C_040_S50

popt, _ = curve_fit(linear_func, T_curie_040_S100[-nfit_S100:], inv_e_040_S100[-nfit_S100:], p0=[1., 1.])
s_curie_040_S100_sc, i_curie_040_S100_sc = popt
s_curie_040_S100, i_curie_040_S100 = s_curie_040_S100_sc / SCALE_INV, i_curie_040_S100_sc / SCALE_INV
C_040_S100 = 1.0 / s_curie_040_S100
T0_040_S100 = -i_curie_040_S100 * C_040_S100

# Volume fits
popt, _ = curve_fit(linear_func, T_vol_020_S0[-nfit_S0:],     V_020_S0[-nfit_S0:],     p0=[1., 1.])
s_vol_020_S0,   i_vol_020_S0   = popt

popt, _ = curve_fit(linear_func, T_vol_020_S50[-nfit_S50:],   V_020_S50[-nfit_S50:],   p0=[1., 1.])
s_vol_020_S50,  i_vol_020_S50  = popt

popt, _ = curve_fit(linear_func, T_vol_020_S100[-nfit_S100:], V_020_S100[-nfit_S100:], p0=[1., 1.])
s_vol_020_S100, i_vol_020_S100 = popt

popt, _ = curve_fit(linear_func, T_vol_040_S0[-nfit_S0:],     V_040_S0[-nfit_S0:],     p0=[1., 1.])
s_vol_040_S0,   i_vol_040_S0   = popt

popt, _ = curve_fit(linear_func, T_vol_040_S50[-nfit_S50:],   V_040_S50[-nfit_S50:],   p0=[1., 1.])
s_vol_040_S50,  i_vol_040_S50  = popt

popt, _ = curve_fit(linear_func, T_vol_040_S100[-nfit_S100:], V_040_S100[-nfit_S100:], p0=[1., 1.])
s_vol_040_S100, i_vol_040_S100 = popt

# -----------------------------------------------------------------------
# PLOTTING SETTINGS

fig_size = (24, 15)

Row_labels = ['(a)', '(b)',  '(c)']
Row_labels_size = 25
Row_labels_type = 'bold'

Xlabel = r'Temperature $[K]$'
Xlabel_size = 24
Xlimit = [(-5, 805), (-5, 805), (-5, 1100)]
Xtick_step = 200

Ylabel = [r'$\varepsilon^{-1}\ (\times 10^{-3})$',
    r'Volume [$\mathrm{\AA}^3$]',
    r'$P_L$ [$\mu\mathrm{C}/\mathrm{cm}^2$]']
Ylabel_size = 24

Col_colors = {'S0': '#1f77b4', 'S50': '#ff7f0e', 'S100': '#2ca02c'}
Adj_color = '#4E4D4D'
Adj_linestyle = '--'

Marker_size = 7
line_width = 2.5
Adj_linewidth = 2.0

legend_size = 20
legend_loc = 'upper center'

Grid_style = '--'
Grid_alpha = 0.1
Grid_width = 1.0

Spine_width = 1.5
Tick_major_width = 2.0
Tick_major_length = 8.0
Tick_minor_width = 1.5
Tick_minor_length = 4.0
Tick_labelsize = 20

Hspace = 0.30
Wspace = 0.20

# y-axis ticks — (ymin, ymax, step)  →  None = escala automática
# Figure x = 0.20
Yticks_020_eps_S0   = (0, 2.4, 0.8)   # ε⁻¹  S0
Yticks_020_eps_S50  = (0, 3, 1)   # ε⁻¹  S50
Yticks_020_eps_S100 = (0, 3.2, 1)   # ε⁻¹  S100

Yticks_020_vol_S0   = (65.6,66.8,0.4)   # Vol  S0
Yticks_020_vol_S50  = (65.8,67,0.4)   # Vol  S50
Yticks_020_vol_S100 = (66,67.3,0.4)   # Vol  S100

Yticks_020_pol_S0   = (5,18,3)   # Pol  S0
Yticks_020_pol_S50  = (5,20,5)   # Pol  S50
Yticks_020_pol_S100 = (5,23,4)   # Pol  S100

# Figure x = 0.40
Yticks_040_eps_S0   = (0.5,1.8,0.4)   # ε⁻¹  S0
Yticks_040_eps_S50  = (0.5,1.8,0.4)   # ε⁻¹  S50
Yticks_040_eps_S100 = (0.5,2.2,0.5)   # ε⁻¹  S100

Yticks_040_vol_S0   = (67.5,68.7,0.4)   # Vol  S0
Yticks_040_vol_S50  = (67.7,68.9,0.4)   # Vol  S50
Yticks_040_vol_S100 = (68,69.4,0.4)   # Vol  S100

Yticks_040_pol_S0   = None   # Pol  S0
Yticks_040_pol_S50  = None   # Pol  S50
Yticks_040_pol_S100 = None   # Pol  S100

# Zoom ranges for x=0.40 volume insets
Zoom_040_S0   = (200.0, 501.0)
Zoom_040_S50  = (300.0, 601.0)
Zoom_040_S100 = (400.0, 701.0)

Inset_xtick_step1 = 100
Inset_xtick_step2 = 100
Inset_xtick_step3 = 100

Inset_spine_width = 1.5
Inset_tick_major_width = 1.5
Inset_tick_major_length = 6.0
Inset_tick_minor_width = 1.0
Inset_tick_minor_length = 3.0
Inset_tick_labelsize = 14

# -----------------------------------------------------------------------
# FIGURE 1 — x = 0.20

fig1, axes1 = plt.subplots(3, 3, figsize=fig_size)

# Y-axis labels per row placed after tight_layout — see bottom of this figure block

# -----------------------------------------------------------------------
# Row 0, Col 0 — Curie S0, x = 0.20
ax00 = axes1[0, 0]
ax00.plot(T_curie_020_S0, inv_e_020_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width, label='S0')
x_fit = np.linspace(T_curie_020_S0[0], T_curie_020_S0[-1], 300)
ax00.plot(x_fit, linear_func(x_fit, s_curie_020_S0_sc, i_curie_020_S0_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax00.set_xlim(Xlimit[0])
ax00.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax00.tick_params(labelbottom=False)
ax00.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax00.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax00.minorticks_on()
ax00.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax00.spines.values():
    spine.set_linewidth(Spine_width)
ax00.legend(fontsize=legend_size, loc=legend_loc)

# Row 0, Col 1 — Curie S50, x = 0.20
ax01 = axes1[0, 1]
ax01.plot(T_curie_020_S50, inv_e_020_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width, label='S50')
x_fit = np.linspace(T_curie_020_S50[0], T_curie_020_S50[-1], 300)
ax01.plot(x_fit, linear_func(x_fit, s_curie_020_S50_sc, i_curie_020_S50_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax01.set_xlim(Xlimit[1])
ax01.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax01.tick_params(labelbottom=False)
ax01.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax01.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax01.minorticks_on()
ax01.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax01.spines.values():
    spine.set_linewidth(Spine_width)
ax01.legend(fontsize=legend_size, loc=legend_loc)

# Row 0, Col 2 — Curie S100, x = 0.20
ax02 = axes1[0, 2]
ax02.plot(T_curie_020_S100, inv_e_020_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width, label='S100')
x_fit = np.linspace(T_curie_020_S100[0], T_curie_020_S100[-1], 300)
ax02.plot(x_fit, linear_func(x_fit, s_curie_020_S100_sc, i_curie_020_S100_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax02.set_xlim(Xlimit[2])
ax02.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax02.tick_params(labelbottom=False)
ax02.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax02.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax02.minorticks_on()
ax02.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax02.spines.values():
    spine.set_linewidth(Spine_width)
ax02.legend(fontsize=legend_size, loc=legend_loc)

# -----------------------------------------------------------------------
# Row 1, Col 0 — Volume S0, x = 0.20
ax10 = axes1[1, 0]
ax10.plot(T_vol_020_S0, V_020_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_020_S0[0], T_vol_020_S0[-1], 300)
ax10.plot(x_fit, linear_func(x_fit, s_vol_020_S0, i_vol_020_S0),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax10.set_xlim(Xlimit[0])
ax10.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax10.tick_params(labelbottom=False)
ax10.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax10.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax10.minorticks_on()
ax10.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax10.spines.values():
    spine.set_linewidth(Spine_width)

# Row 1, Col 1 — Volume S50, x = 0.20
ax11 = axes1[1, 1]
ax11.plot(T_vol_020_S50, V_020_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_020_S50[0], T_vol_020_S50[-1], 300)
ax11.plot(x_fit, linear_func(x_fit, s_vol_020_S50, i_vol_020_S50),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax11.set_xlim(Xlimit[1])
ax11.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax11.tick_params(labelbottom=False)
ax11.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax11.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax11.minorticks_on()
ax11.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax11.spines.values():
    spine.set_linewidth(Spine_width)

# Row 1, Col 2 — Volume S100, x = 0.20
ax12 = axes1[1, 2]
ax12.plot(T_vol_020_S100, V_020_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_020_S100[0], T_vol_020_S100[-1], 300)
ax12.plot(x_fit, linear_func(x_fit, s_vol_020_S100, i_vol_020_S100),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax12.set_xlim(Xlimit[2])
ax12.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax12.tick_params(labelbottom=False)
ax12.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax12.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax12.minorticks_on()
ax12.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax12.spines.values():
    spine.set_linewidth(Spine_width)

# -----------------------------------------------------------------------
# Row 2, Col 0 — Polarization S0, x = 0.20  (bottom row — x labels shown)
ax20 = axes1[2, 0]
ax20.plot(T_pol_020_S0, p_020_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width)
ax20.set_xlim(Xlimit[0])
ax20.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax20.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax20.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax20.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax20.minorticks_on()
ax20.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax20.spines.values():
    spine.set_linewidth(Spine_width)


# Row 2, Col 1 — Polarization S50, x = 0.20
ax21 = axes1[2, 1]
ax21.plot(T_pol_020_S50, p_020_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width)
ax21.set_xlim(Xlimit[1])
ax21.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax21.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax21.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax21.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax21.minorticks_on()
ax21.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax21.spines.values():
    spine.set_linewidth(Spine_width)

# Row 2, Col 2 — Polarization S100, x = 0.20
ax22 = axes1[2, 2]
ax22.plot(T_pol_020_S100, p_020_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width)
ax22.set_xlim(Xlimit[2])
ax22.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax22.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax22.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax22.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax22.minorticks_on()
ax22.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax22.spines.values():
    spine.set_linewidth(Spine_width)

for ax, yticks in zip(axes1.flat, [
    Yticks_020_eps_S0,   Yticks_020_eps_S50,   Yticks_020_eps_S100,
    Yticks_020_vol_S0,   Yticks_020_vol_S50,   Yticks_020_vol_S100,
    Yticks_020_pol_S0,   Yticks_020_pol_S50,   Yticks_020_pol_S100,
]):
    if yticks is not None:
        ax.set_yticks(np.arange(yticks[0], yticks[1] + yticks[2], yticks[2]))
        ax.set_ylim(yticks[0], yticks[1])

plt.tight_layout(rect=[0.06, 0.02, 0.98, 0.97])
plt.subplots_adjust(hspace=Hspace, wspace=Wspace)

# Single Y-axis labels per row (left side), placed after layout
for ax_row, ylabel in zip([ax00, ax10, ax20], Ylabel):
    pos = ax_row.get_position()
    fig1.text(0.05, (pos.y0 + pos.y1) / 2, ylabel, va='center', ha='center',
              rotation='vertical', fontsize=Ylabel_size)

# Row labels centered above each row
for ax_l, ax_r, label in zip([ax00, ax10, ax20], [ax02, ax12, ax22], Row_labels):
    pos_l = ax_l.get_position()
    pos_r = ax_r.get_position()
    x_center = (pos_l.x0 + pos_r.x1) / 2
    fig1.text(x_center, pos_l.y1 + 0.01, label, ha='center', va='bottom',
              fontsize=Row_labels_size, fontweight=Row_labels_type)

fig1.savefig('Burns-0.20.png', dpi=300, bbox_inches='tight')
plt.show()

# -----------------------------------------------------------------------
# FIGURE 2 — x = 0.40

fig2, axes2 = plt.subplots(3, 3, figsize=fig_size)

# -----------------------------------------------------------------------
# Row 0, Col 0 — Curie S0, x = 0.40
ax00 = axes2[0, 0]
ax00.plot(T_curie_040_S0, inv_e_040_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width, label='S0')
x_fit = np.linspace(T_curie_040_S0[0], T_curie_040_S0[-1], 300)
ax00.plot(x_fit, linear_func(x_fit, s_curie_040_S0_sc, i_curie_040_S0_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax00.set_xlim(Xlimit[0])
ax00.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax00.tick_params(labelbottom=False)
ax00.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax00.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax00.minorticks_on()
ax00.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax00.spines.values():
    spine.set_linewidth(Spine_width)
ax00.legend(fontsize=legend_size, loc=legend_loc)

# Row 0, Col 1 — Curie S50, x = 0.40
ax01 = axes2[0, 1]
ax01.plot(T_curie_040_S50, inv_e_040_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width, label='S50')
x_fit = np.linspace(T_curie_040_S50[0], T_curie_040_S50[-1], 300)
ax01.plot(x_fit, linear_func(x_fit, s_curie_040_S50_sc, i_curie_040_S50_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax01.set_xlim(Xlimit[1])
ax01.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax01.tick_params(labelbottom=False)
ax01.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax01.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax01.minorticks_on()
ax01.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax01.spines.values():
    spine.set_linewidth(Spine_width)
ax01.legend(fontsize=legend_size, loc=legend_loc)

# Row 0, Col 2 — Curie S100, x = 0.40
ax02 = axes2[0, 2]
ax02.plot(T_curie_040_S100, inv_e_040_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width, label='S100')
x_fit = np.linspace(T_curie_040_S100[0], T_curie_040_S100[-1], 300)
ax02.plot(x_fit, linear_func(x_fit, s_curie_040_S100_sc, i_curie_040_S100_sc),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax02.set_xlim(Xlimit[2])
ax02.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax02.tick_params(labelbottom=False)
ax02.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax02.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax02.minorticks_on()
ax02.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax02.spines.values():
    spine.set_linewidth(Spine_width)
ax02.legend(fontsize=legend_size, loc=legend_loc)

# -----------------------------------------------------------------------
# Row 1, Col 0 — Volume S0, x = 0.40  (with inset zoom)
ax10 = axes2[1, 0]
ax10.plot(T_vol_040_S0, V_040_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_040_S0[0], T_vol_040_S0[-1], 300)
ax10.plot(x_fit, linear_func(x_fit, s_vol_040_S0, i_vol_040_S0),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax10.set_xlim(Xlimit[0])
ax10.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax10.tick_params(labelbottom=False)
ax10.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax10.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax10.minorticks_on()
ax10.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax10.spines.values():
    spine.set_linewidth(Spine_width)

# inset
tmin, tmax = Zoom_040_S0
mask = (T_vol_040_S0 >= tmin) & (T_vol_040_S0 <= tmax)
axins10 = ax10.inset_axes([0.55, 0.12, 0.40, 0.35])
axins10.plot(T_vol_040_S0, V_040_S0, marker='o', linestyle='-',
             color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width)
xzoom = np.linspace(tmin, tmax, 100)
axins10.plot(xzoom, linear_func(xzoom, s_vol_040_S0, i_vol_040_S0),
             linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
V_zoom = V_040_S0[mask]
vpad = 0.06 * (V_zoom.max() - V_zoom.min()) if V_zoom.size > 1 else 0.001
axins10.set_xlim(tmin, tmax)
axins10.set_ylim(V_zoom.min() - vpad, V_zoom.max() + vpad)
axins10.set_xticks(np.arange(tmin, tmax + 1, Inset_xtick_step1))
axins10.minorticks_on()
axins10.grid(True, linestyle=Grid_style, alpha=Grid_alpha)
axins10.tick_params(labelsize=Inset_tick_labelsize, length=Inset_tick_major_length, width=Inset_tick_major_width, axis='x')
axins10.tick_params(axis='y', which='both', left=False, labelleft=False)
for spine in axins10.spines.values():
    spine.set_linewidth(Inset_spine_width)

# Row 1, Col 1 — Volume S50, x = 0.40  (with inset zoom)
ax11 = axes2[1, 1]
ax11.plot(T_vol_040_S50, V_040_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_040_S50[0], T_vol_040_S50[-1], 300)
ax11.plot(x_fit, linear_func(x_fit, s_vol_040_S50, i_vol_040_S50),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax11.set_xlim(Xlimit[1])
ax11.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax11.tick_params(labelbottom=False)
ax11.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax11.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax11.minorticks_on()
ax11.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax11.spines.values():
    spine.set_linewidth(Spine_width)

# inset
tmin, tmax = Zoom_040_S50
mask = (T_vol_040_S50 >= tmin) & (T_vol_040_S50 <= tmax)
axins11 = ax11.inset_axes([0.55, 0.12, 0.40, 0.35])
axins11.plot(T_vol_040_S50, V_040_S50, marker='o', linestyle='-',
             color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width)
xzoom = np.linspace(tmin, tmax, 100)
axins11.plot(xzoom, linear_func(xzoom, s_vol_040_S50, i_vol_040_S50),
             linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
V_zoom = V_040_S50[mask]
vpad = 0.06 * (V_zoom.max() - V_zoom.min()) if V_zoom.size > 1 else 0.001
axins11.set_xlim(tmin, tmax)
axins11.set_ylim(V_zoom.min() - vpad, V_zoom.max() + vpad)
axins11.set_xticks(np.arange(tmin, tmax + 1, Inset_xtick_step2))
axins11.minorticks_on()
axins11.grid(True, linestyle=Grid_style, alpha=Grid_alpha)
axins11.tick_params(labelsize=Inset_tick_labelsize, length=Inset_tick_major_length, width=Inset_tick_major_width, axis='x')
axins11.tick_params(axis='y', which='both', left=False, labelleft=False)
for spine in axins11.spines.values():
    spine.set_linewidth(Inset_spine_width)

# Row 1, Col 2 — Volume S100, x = 0.40  (with inset zoom)
ax12 = axes2[1, 2]
ax12.plot(T_vol_040_S100, V_040_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width)
x_fit = np.linspace(T_vol_040_S100[0], T_vol_040_S100[-1], 300)
ax12.plot(x_fit, linear_func(x_fit, s_vol_040_S100, i_vol_040_S100),
          linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
ax12.set_xlim(Xlimit[2])
ax12.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax12.tick_params(labelbottom=False)
ax12.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax12.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax12.minorticks_on()
ax12.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax12.spines.values():
    spine.set_linewidth(Spine_width)
# inset
tmin, tmax = Zoom_040_S100
mask = (T_vol_040_S100 >= tmin) & (T_vol_040_S100 <= tmax)
axins12 = ax12.inset_axes([0.55, 0.12, 0.40, 0.35])
axins12.plot(T_vol_040_S100, V_040_S100, marker='o', linestyle='-',
             color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width)
xzoom = np.linspace(tmin, tmax, 100)
axins12.plot(xzoom, linear_func(xzoom, s_vol_040_S100, i_vol_040_S100),
             linestyle=Adj_linestyle, color=Adj_color, linewidth=Adj_linewidth)
V_zoom = V_040_S100[mask]
vpad = 0.06 * (V_zoom.max() - V_zoom.min()) if V_zoom.size > 1 else 0.001
axins12.set_xlim(tmin, tmax)
axins12.set_ylim(V_zoom.min() - vpad, V_zoom.max() + vpad)
axins12.set_xticks(np.arange(tmin, tmax + 1, Inset_xtick_step3))
axins12.minorticks_on()
axins12.grid(True, linestyle=Grid_style, alpha=Grid_alpha)
axins12.tick_params(labelsize=Inset_tick_labelsize, length=Inset_tick_major_length, width=Inset_tick_major_width, axis='x')
axins12.tick_params(axis='y', which='both', left=False, labelleft=False)
for spine in axins12.spines.values():
    spine.set_linewidth(Inset_spine_width)

# -----------------------------------------------------------------------
# Row 2, Col 0 — Polarization S0, x = 0.40  (bottom row — x labels shown)
ax20 = axes2[2, 0]
ax20.plot(T_pol_040_S0, p_040_S0, marker='o', linestyle='-',
          color=Col_colors['S0'], markersize=Marker_size, linewidth=line_width)
ax20.set_xlim(Xlimit[0])
ax20.set_xticks(np.arange(0, Xlimit[0][1], Xtick_step))
ax20.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax20.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax20.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax20.minorticks_on()
ax20.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax20.spines.values():
    spine.set_linewidth(Spine_width)


# Row 2, Col 1 — Polarization S50, x = 0.40
ax21 = axes2[2, 1]
ax21.plot(T_pol_040_S50, p_040_S50, marker='o', linestyle='-',
          color=Col_colors['S50'], markersize=Marker_size, linewidth=line_width)
ax21.set_xlim(Xlimit[1])
ax21.set_xticks(np.arange(0, Xlimit[1][1], Xtick_step))
ax21.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax21.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax21.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax21.minorticks_on()
ax21.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax21.spines.values():
    spine.set_linewidth(Spine_width)


# Row 2, Col 2 — Polarization S100, x = 0.40
ax22 = axes2[2, 2]
ax22.plot(T_pol_040_S100, p_040_S100, marker='o', linestyle='-',
          color=Col_colors['S100'], markersize=Marker_size, linewidth=line_width)
ax22.set_xlim(Xlimit[2])
ax22.set_xticks(np.arange(0, Xlimit[2][1], Xtick_step))
ax22.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax22.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax22.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax22.minorticks_on()
ax22.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax22.spines.values():
    spine.set_linewidth(Spine_width)

for ax, yticks in zip(axes2.flat, [
    Yticks_040_eps_S0,   Yticks_040_eps_S50,   Yticks_040_eps_S100,
    Yticks_040_vol_S0,   Yticks_040_vol_S50,   Yticks_040_vol_S100,
    Yticks_040_pol_S0,   Yticks_040_pol_S50,   Yticks_040_pol_S100,
]):
    if yticks is not None:
        ax.set_yticks(np.arange(yticks[0], yticks[1] + yticks[2], yticks[2]))
        ax.set_ylim(yticks[0], yticks[1])

plt.tight_layout(rect=[0.06, 0.02, 0.98, 0.97])
plt.subplots_adjust(hspace=Hspace, wspace=Wspace)

# Single Y-axis labels per row (left side), placed after layout
for ax_row, ylabel in zip([ax00, ax10, ax20], Ylabel):
    pos = ax_row.get_position()
    fig2.text(0.05, (pos.y0 + pos.y1) / 2, ylabel, va='center', ha='center',
              rotation='vertical', fontsize=Ylabel_size)

# Row labels centered above each row
for ax_l, ax_r, label in zip([ax00, ax10, ax20], [ax02, ax12, ax22], Row_labels):
    pos_l = ax_l.get_position()
    pos_r = ax_r.get_position()
    x_center = (pos_l.x0 + pos_r.x1) / 2
    fig2.text(x_center, pos_l.y1 + 0.01, label, ha='center', va='bottom',
              fontsize=Row_labels_size, fontweight=Row_labels_type)

fig2.savefig('Burns-0.40.png', dpi=300, bbox_inches='tight')
plt.show()

# -----------------------------------------------------------------------
# SAVE FIT RESULTS

def fmt(x):
    try:
        return f'{float(x):.6e}'
    except Exception:
        return 'nan'

for conc, results in [
    ("0.20", [
        ('S0',   s_curie_020_S0,   i_curie_020_S0,   C_020_S0,   T0_020_S0),
        ('S50',  s_curie_020_S50,  i_curie_020_S50,  C_020_S50,  T0_020_S50),
        ('S100', s_curie_020_S100, i_curie_020_S100, C_020_S100, T0_020_S100),
    ]),
    ("0.40", [
        ('S0',   s_curie_040_S0,   i_curie_040_S0,   C_040_S0,   T0_040_S0),
        ('S50',  s_curie_040_S50,  i_curie_040_S50,  C_040_S50,  T0_040_S50),
        ('S100', s_curie_040_S100, i_curie_040_S100, C_040_S100, T0_040_S100),
    ]),
]:
    with open(f'Burns-{conc}-Curie.txt', 'w') as f:
        f.write(f'Concentracion: {conc}\n')
        f.write('seg\tslope\tintercept\tCurie_const\tT0\n')
        for r in results:
            f.write(f'{r[0]}\t{fmt(r[1])}\t{fmt(r[2])}\t{fmt(r[3])}\t{fmt(r[4])}\n')

for conc, results in [
    ("0.20", [
        ('S0',   s_vol_020_S0,   i_vol_020_S0),
        ('S50',  s_vol_020_S50,  i_vol_020_S50),
        ('S100', s_vol_020_S100, i_vol_020_S100),
    ]),
    ("0.40", [
        ('S0',   s_vol_040_S0,   i_vol_040_S0),
        ('S50',  s_vol_040_S50,  i_vol_040_S50),
        ('S100', s_vol_040_S100, i_vol_040_S100),
    ]),
]:
    with open(f'Burns-{conc}-Volumen.txt', 'w') as f:
        f.write(f'Concentracion: {conc}\n')
        f.write('seg\tslope\tintercept\n')
        for r in results:
            f.write(f'{r[0]}\t{fmt(r[1])}\t{fmt(r[2])}\n')
