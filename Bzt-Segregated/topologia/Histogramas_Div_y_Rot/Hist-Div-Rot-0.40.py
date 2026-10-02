import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

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
SCALE = 1e8  # los datos se grafican en unidades de 10^8 para evitar el offset "1e8" de matplotlib

# PLOT SETTINGS
Fig_size=(18, 21)
DPI=300
Outfile=f'Hist-Div-Rot-{COMPOSICION}.png'

Suptitle=f''
Suptitle_size=13
Suptitle_type='bold'

Row_labels=['a) S0', 'b) S50', 'c) S100']
Titles_div=''
Titles_curl_mod=''
Titles_size=20
Titles_type='bold'

Xlabel_div=r'∇·P ($\cdot10^8~\mu C/cm^3)$'
Xlabel_curl_mod=r'$|\nabla\times P|$ ($\cdot10^8~\mu C/cm^3)$'
Xlabel_size=24
Xlabel_type='normal'

Ylabel_size=24
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
Tick_label_size=22


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
    ax.ticklabel_format(style='plain', axis='x', useOffset=False)


# ── Carga y separación por especie ────────────────
datos_div, medias_div = {}, {}
datos_curl_mod, medias_curl_mod = {}, {}
for sistema in SISTEMAS:
    id_grid = np.load(FILE_ID[sistema])
    zr, ti = id_grid == 1, id_grid == 2

    div = np.load(FILE_DIV[sistema])
    datos_div[sistema] = (div[zr], div[ti])
    medias_div[sistema] = (div[zr].mean(), div[ti].mean())

    curl_x = np.load(FILE_CURL_X[sistema])
    curl_y = np.load(FILE_CURL_Y[sistema])
    curl_z = np.load(FILE_CURL_Z[sistema])
    curl_mod = np.sqrt(curl_x**2 + curl_y**2 + curl_z**2)

    datos_curl_mod[sistema] = (curl_mod[zr], curl_mod[ti])
    medias_curl_mod[sistema] = (curl_mod[zr].mean(), curl_mod[ti].mean())

print(f'\n=== Media de la divergencia ∇·P — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Zr":>14}{"Ti":>14}')
for sistema in SISTEMAS:
    m_zr, m_ti = medias_div[sistema]
    print(f'{sistema:<8}{m_zr:>14.4e}{m_ti:>14.4e}')

print(f'\n=== Media del módulo de ∇×P — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Zr":>14}{"Ti":>14}')
for sistema in SISTEMAS:
    m_zr, m_ti = medias_curl_mod[sistema]
    print(f'{sistema:<8}{m_zr:>14.4e}{m_ti:>14.4e}')

# Rango común de bins por columna (la divergencia puede ser negativa, el módulo del rotor no)
todos_div = np.concatenate([v for par in datos_div.values() for v in par])
vlim_div = np.percentile(np.abs(todos_div), 99)
bins_div = np.linspace(-vlim_div, vlim_div, BINS + 1) / SCALE

todos_curl_mod = np.concatenate([v for par in datos_curl_mod.values() for v in par])
vlim_curl_mod = np.percentile(todos_curl_mod, 99)
bins_curl_mod = np.linspace(0, vlim_curl_mod, BINS + 1) / SCALE

# ── Figura: 3 filas (sistemas) x 2 columnas (div, |curl|) ──
fig, axes = plt.subplots(3, 2, figsize=Fig_size, gridspec_kw={'hspace': 0.21})
fig.patch.set_facecolor('white')

Ylim_div = (0, 900)
Ylim_curl_mod = (0, 650)

columnas = [
    (datos_div, medias_div, bins_div, Titles_div, Xlabel_div, True, Ylim_div),
    (datos_curl_mod, medias_curl_mod, bins_curl_mod, Titles_curl_mod, Xlabel_curl_mod, False, Ylim_curl_mod),
]

for row, sistema in enumerate(SISTEMAS):
    for col, (datos, medias, bins_c, titulo, xlabel, zero_line, ylim) in enumerate(columnas):
        ax = axes[row, col]
        c_zr, c_ti = datos[sistema]
        m_zr, m_ti = medias[sistema]
        c_zr_s, c_ti_s = c_zr / SCALE, c_ti / SCALE
        m_zr_s, m_ti_s = m_zr / SCALE, m_ti / SCALE
        ax.hist(c_zr_s, bins=bins_c, color=COLOR_ZR, alpha=Hist_alpha,
                label=fr'Zr (⟨⟩={m_zr_s:.2f} $\cdot10^{{8}}$)')
        ax.hist(c_ti_s, bins=bins_c, color=COLOR_TI, alpha=Hist_alpha,
                label=fr'Ti (⟨⟩={m_ti_s:.2f} $\cdot10^{{8}}$)')
        ax.set_ylim(*ylim)
        ax.yaxis.set_major_locator(MultipleLocator(200))
        if zero_line:
            ax.axvline(0, color='black', linewidth=Zero_line_width, linestyle=Zero_line_style)
        ax.axvline(m_zr_s, color=COLOR_ZR, linewidth=Mean_line_width, linestyle='-')
        ax.axvline(m_ti_s, color=COLOR_TI, linewidth=Mean_line_width, linestyle='-')
        _style(ax, titulo if row == 0 else '')
        if row == len(SISTEMAS) - 1:
            ax.set_xlabel(xlabel, color='black', fontsize=Xlabel_size, fontweight=Xlabel_type)
        else:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        leg = ax.legend(facecolor='white', edgecolor='black', fontsize=legend_size, loc=legend_loc)
        for text in leg.get_texts():
            text.set_color('black')
            text.set_fontweight(legend_type)

fig.suptitle(Suptitle, color='black', fontsize=Suptitle_size, fontweight=Suptitle_type)
fig.tight_layout(rect=[0.005, 0, 1, 0.96])

fig.canvas.draw()
fig.supylabel('Frequency', x=0.02, color='black', fontsize=Ylabel_size, fontweight=Ylabel_type)
gap_height = (axes[0, 0].get_position().y0
              - max(axes[1, 0].get_position().y1, axes[1, 1].get_position().y1))
for row, label in enumerate(Row_labels):
    bbox_left = axes[row, 0].get_position()
    bbox_right = axes[row, 1].get_position()
    x_center = (bbox_left.x0 + bbox_right.x1) / 2
    row_top = max(bbox_left.y1, bbox_right.y1)
    if row == 0:
        y_center = row_top + gap_height / 2
    else:
        space_top = min(axes[row - 1, 0].get_position().y0, axes[row - 1, 1].get_position().y0)
        y_center = (row_top + space_top) / 2
    fig.text(x_center, y_center, label, ha='center', va='center',
              color='black', fontsize=Titles_size + 4, fontweight=Titles_type)

plt.savefig(Outfile, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"Guardado: {Outfile}")
