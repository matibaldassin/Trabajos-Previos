import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# INPUT FILES
COMPOSICION = '0.40'
SISTEMAS = ['S0', 'S50', 'S100']
FILE_DIV    = {s: f'../{COMPOSICION}-{s}/Cooling/div.npy' for s in SISTEMAS}
FILE_CURL_X = {s: f'../{COMPOSICION}-{s}/Cooling/curl_x.npy' for s in SISTEMAS}
FILE_CURL_Y = {s: f'../{COMPOSICION}-{s}/Cooling/curl_y.npy' for s in SISTEMAS}
FILE_CURL_Z = {s: f'../{COMPOSICION}-{s}/Cooling/curl_z.npy' for s in SISTEMAS}
FILE_ID     = {s: f'../{COMPOSICION}-{s}/Cooling/id_grid.npy' for s in SISTEMAS}

COLOR_ZR = '#FC0D04'
COLOR_TI = '#6B68FC'
BINS = 60

# PLOT SETTINGS
Fig_size=(24, 14)
DPI=300
Outfile=f'Hist-Div-RotAbs-{COMPOSICION}.png'

Suptitle=f''
Suptitle_size=13
Suptitle_type='bold'

Row_labels=['S0', 'S50', 'S100']
Titles_div='∇·P'
Titles_curl_x='|(∇×P)ₓ|'
Titles_curl_y='|(∇×P)ᵧ|'
Titles_curl_z='|(∇×P)_z|'
Titles_size=20
Titles_type='bold'

Xlabel_div='∇·P'
Xlabel_curl_x='|(∇×P)ₓ|'
Xlabel_curl_y='|(∇×P)ᵧ|'
Xlabel_curl_z='|(∇×P)_z|'
Xlabel_size=20
Xlabel_type='normal'

Ylabel_size=20
Ylabel_type='normal'

Hist_alpha=0.2
Mean_line_width=2
Zero_line_width=1
Zero_line_style='--'

legend_size=16
legend_type='normal'
legend_loc='best'

Grid_style='--'
Grid_alpha=0
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=20


# ── Estilo (coherente con 02_overview.py) ─────────

def _style(ax, title):
    ax.set_title(title, color='black', fontsize=Titles_size, fontweight=Titles_type)
    ax.set_facecolor('white')
    ax.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width, color='black')
    ax.tick_params(which='major', colors='black', labelsize=Tick_label_size,
                    width=Tick_major_width, length=Tick_major_length)
    ax.tick_params(which='minor', colors='black',
                    width=Tick_minor_width, length=Tick_minor_length)
    for sp in ax.spines.values():
        sp.set_color('black')
        sp.set_linewidth(Spine_width)


# ── Carga y separación por especie ────────────────
datos_div, medias_div = {}, {}
datos_curlx, medias_curlx = {}, {}
datos_curly, medias_curly = {}, {}
datos_curlz, medias_curlz = {}, {}
for sistema in SISTEMAS:
    id_grid = np.load(FILE_ID[sistema])
    zr, ti = id_grid == 1, id_grid == 2

    div = np.load(FILE_DIV[sistema])
    datos_div[sistema] = (div[zr], div[ti])
    medias_div[sistema] = (div[zr].mean(), div[ti].mean())

    curl_x = np.abs(np.load(FILE_CURL_X[sistema]))
    curl_y = np.abs(np.load(FILE_CURL_Y[sistema]))
    curl_z = np.abs(np.load(FILE_CURL_Z[sistema]))

    datos_curlx[sistema] = (curl_x[zr], curl_x[ti])
    medias_curlx[sistema] = (curl_x[zr].mean(), curl_x[ti].mean())
    datos_curly[sistema] = (curl_y[zr], curl_y[ti])
    medias_curly[sistema] = (curl_y[zr].mean(), curl_y[ti].mean())
    datos_curlz[sistema] = (curl_z[zr], curl_z[ti])
    medias_curlz[sistema] = (curl_z[zr].mean(), curl_z[ti].mean())

print(f'\n=== Media de la divergencia ∇·P — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Zr":>14}{"Ti":>14}')
for sistema in SISTEMAS:
    m_zr, m_ti = medias_div[sistema]
    print(f'{sistema:<8}{m_zr:>14.4e}{m_ti:>14.4e}')

# Carga de polarización: div.npy ya esta en uC/cm^3 (rho = div P), asi que
# integrar sobre el volumen es sumar cada celda * V_celda. V_celda = 64 A^3
# = 64*(1e-8 cm)^3, el volumen promedio de celda unidad usado para dx=4e-8.
V_CELL_A3 = 64.0
V_CELL = V_CELL_A3 * (1e-8) ** 3  # cm^3

print(f'\n=== Carga de polarización Q = Σ(∇·P)·V_celda, V_celda={V_CELL_A3:g} Å³ — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Q_Zr (uC)":>16}{"Q_Ti (uC)":>16}{"Q_total (uC)":>16}')
for sistema in SISTEMAS:
    div_zr, div_ti = datos_div[sistema]
    q_zr = div_zr.sum() * V_CELL
    q_ti = div_ti.sum() * V_CELL
    print(f'{sistema:<8}{q_zr:>16.4e}{q_ti:>16.4e}{q_zr + q_ti:>16.4e}')

print(f'\n=== Media de |componentes de ∇×P| — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"":>14}{"Zr":>14}{"Ti":>14}')
for sistema in SISTEMAS:
    mx_zr, mx_ti = medias_curlx[sistema]
    my_zr, my_ti = medias_curly[sistema]
    mz_zr, mz_ti = medias_curlz[sistema]
    print(f'{sistema:<8}{"|(∇×P)x|":>14}{mx_zr:>14.4e}{mx_ti:>14.4e}')
    print(f'{"":<8}{"|(∇×P)y|":>14}{my_zr:>14.4e}{my_ti:>14.4e}')
    print(f'{"":<8}{"|(∇×P)z|":>14}{mz_zr:>14.4e}{mz_ti:>14.4e}')

# Rango común de bins por columna (divergencia puede ser negativa; |rotor| no)
todos_div = np.concatenate([v for par in datos_div.values() for v in par])
vlim_div = np.percentile(np.abs(todos_div), 99)
bins_div = np.linspace(-vlim_div, vlim_div, BINS + 1)

todos_curlx = np.concatenate([v for par in datos_curlx.values() for v in par])
vlim_curlx = np.percentile(todos_curlx, 99)
bins_curlx = np.linspace(0, vlim_curlx, BINS + 1)

todos_curly = np.concatenate([v for par in datos_curly.values() for v in par])
vlim_curly = np.percentile(todos_curly, 99)
bins_curly = np.linspace(0, vlim_curly, BINS + 1)

todos_curlz = np.concatenate([v for par in datos_curlz.values() for v in par])
vlim_curlz = np.percentile(todos_curlz, 99)
bins_curlz = np.linspace(0, vlim_curlz, BINS + 1)

# ── Figura: 3 filas (sistemas) x 4 columnas (div, |curl_x|, |curl_y|, |curl_z|) ──
fig, axes = plt.subplots(3, 4, figsize=Fig_size)
fig.patch.set_facecolor('white')

columnas = [
    (datos_div, medias_div, bins_div, Titles_div, Xlabel_div, True),
    (datos_curlx, medias_curlx, bins_curlx, Titles_curl_x, Xlabel_curl_x, False),
    (datos_curly, medias_curly, bins_curly, Titles_curl_y, Xlabel_curl_y, False),
    (datos_curlz, medias_curlz, bins_curlz, Titles_curl_z, Xlabel_curl_z, False),
]

for row, sistema in enumerate(SISTEMAS):
    for col, (datos, medias, bins_c, titulo, xlabel, con_signo) in enumerate(columnas):
        ax = axes[row, col]
        c_zr, c_ti = datos[sistema]
        m_zr, m_ti = medias[sistema]
        ax.hist(c_zr, bins=bins_c, color=COLOR_ZR, alpha=Hist_alpha,
                label=f'Zr (⟨⟩={m_zr:.2e})')
        ax.hist(c_ti, bins=bins_c, color=COLOR_TI, alpha=Hist_alpha,
                label=f'Ti (⟨⟩={m_ti:.2e})')
        if con_signo:
            ax.axvline(0, color='black', linewidth=Zero_line_width, linestyle=Zero_line_style)
        ax.axvline(m_zr, color=COLOR_ZR, linewidth=Mean_line_width, linestyle='-')
        ax.axvline(m_ti, color=COLOR_TI, linewidth=Mean_line_width, linestyle='-')
        _style(ax, titulo if row == 0 else '')
        if row == len(SISTEMAS) - 1:
            ax.set_xlabel(xlabel, color='black', fontsize=Xlabel_size, fontweight=Xlabel_type)
        if col == 0:
            ax.set_ylabel(Row_labels[row], color='black', fontsize=Ylabel_size, fontweight=Ylabel_type)
        leg = ax.legend(facecolor='white', edgecolor='black', fontsize=legend_size, loc=legend_loc)
        for text in leg.get_texts():
            text.set_color('black')
            text.set_fontweight(legend_type)

fig.suptitle(Suptitle, color='black', fontsize=Suptitle_size, fontweight=Suptitle_type)
fig.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig(Outfile, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"Guardado: {Outfile}")
