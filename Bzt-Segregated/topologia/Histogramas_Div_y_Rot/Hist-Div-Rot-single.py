import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# INPUT FILES
COMPOSICION = '0.40'
SISTEMA = 'S0'
ETAPA = 'Cooling'
FILE_DIV  = f'div.npy'
FILE_CURL_X = f'curl_x.npy'
FILE_CURL_Y = f'curl_y.npy'
FILE_CURL_Z = f'curl_z.npy'
FILE_ID   = f'id_grid.npy'

COLOR_ZR = '#FC0D04'
COLOR_TI = '#6B68FC'
BINS = 60

# PLOT SETTINGS
Fig_size=(20, 5)
DPI=300
Outfile=f'Hist-Div-Rot-{COMPOSICION}-{SISTEMA}.png'

Suptitle=f''
Suptitle_size=13
Suptitle_type='bold'

Title_div='a) ∇·P'
Title_curl_x='b) (∇×P)ₓ'
Title_curl_y='c) (∇×P)ᵧ'
Title_curl_z='d) (∇×P)_z'
Titles_size=20
Titles_type='bold'

Xlabel_div='∇·P'
Xlabel_curl_x='(∇×P)ₓ'
Xlabel_curl_y='(∇×P)ᵧ'
Xlabel_curl_z='(∇×P)_z'
Xlabel_size=20
Xlabel_type='normal'

Ylabel_div='Frecuencia'
Ylabel_curl='Frecuencia'
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
id_grid = np.load(FILE_ID)
zr, ti = id_grid == 1, id_grid == 2

div = np.load(FILE_DIV)
div_zr, div_ti = div[zr], div[ti]
m_div_zr, m_div_ti = div_zr.mean(), div_ti.mean()

curl_x = np.load(FILE_CURL_X)
curl_y = np.load(FILE_CURL_Y)
curl_z = np.load(FILE_CURL_Z)

curlx_zr, curlx_ti = curl_x[zr], curl_x[ti]
curly_zr, curly_ti = curl_y[zr], curl_y[ti]
curlz_zr, curlz_ti = curl_z[zr], curl_z[ti]

m_curlx_zr, m_curlx_ti = curlx_zr.mean(), curlx_ti.mean()
m_curly_zr, m_curly_ti = curly_zr.mean(), curly_ti.mean()
m_curlz_zr, m_curlz_ti = curlz_zr.mean(), curlz_ti.mean()

print(f'\n=== Media de la divergencia ∇·P — x={COMPOSICION}, {SISTEMA}, {ETAPA} ===')
print(f'{"Zr":>14}{"Ti":>14}')
print(f'{m_div_zr:>14.4e}{m_div_ti:>14.4e}')

# Carga de polarización: div.npy ya esta en uC/cm^3 (rho = div P), asi que
# integrar sobre el volumen es sumar cada celda * V_celda. V_celda = 64 A^3
# = 64*(1e-8 cm)^3, el volumen promedio de celda unidad usado para dx=4e-8.
V_CELL_A3 = 64.0
V_CELL = V_CELL_A3 * (1e-8) ** 3  # cm^3
q_zr = div_zr.sum() * V_CELL
q_ti = div_ti.sum() * V_CELL

print(f'\n=== Carga de polarización Q = Σ(∇·P)·V_celda, V_celda={V_CELL_A3:g} Å³ — x={COMPOSICION}, {SISTEMA}, {ETAPA} ===')
print(f'{"Q_Zr (uC)":>16}{"Q_Ti (uC)":>16}{"Q_total (uC)":>16}')
print(f'{q_zr:>16.4e}{q_ti:>16.4e}{q_zr + q_ti:>16.4e}')

print(f'\n=== Media de las componentes de ∇×P — x={COMPOSICION}, {SISTEMA}, {ETAPA} ===')
print(f'{"":>14}{"Zr":>14}{"Ti":>14}')
print(f'{"(∇×P)x":>14}{m_curlx_zr:>14.4e}{m_curlx_ti:>14.4e}')
print(f'{"(∇×P)y":>14}{m_curly_zr:>14.4e}{m_curly_ti:>14.4e}')
print(f'{"(∇×P)z":>14}{m_curlz_zr:>14.4e}{m_curlz_ti:>14.4e}')

# Rango de bins (divergencia y componentes del rotor pueden ser negativas)
vlim_div = np.percentile(np.abs(np.concatenate([div_zr, div_ti])), 99)
bins_div = np.linspace(-vlim_div, vlim_div, BINS + 1)

vlim_curlx = np.percentile(np.abs(np.concatenate([curlx_zr, curlx_ti])), 99)
bins_curlx = np.linspace(-vlim_curlx, vlim_curlx, BINS + 1)

vlim_curly = np.percentile(np.abs(np.concatenate([curly_zr, curly_ti])), 99)
bins_curly = np.linspace(-vlim_curly, vlim_curly, BINS + 1)

vlim_curlz = np.percentile(np.abs(np.concatenate([curlz_zr, curlz_ti])), 99)
bins_curlz = np.linspace(-vlim_curlz, vlim_curlz, BINS + 1)

# ── Figura ─────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=Fig_size)
fig.patch.set_facecolor('white')

# divergencia
ax = axes[0]
ax.hist(div_zr, bins=bins_div, color=COLOR_ZR, alpha=Hist_alpha,
        label=f'Zr (⟨⟩={m_div_zr:.2e})')
ax.hist(div_ti, bins=bins_div, color=COLOR_TI, alpha=Hist_alpha,
        label=f'Ti (⟨⟩={m_div_ti:.2e})')
ax.axvline(0, color='black', linewidth=Zero_line_width, linestyle=Zero_line_style)
ax.axvline(m_div_zr, color=COLOR_ZR, linewidth=Mean_line_width, linestyle='-')
ax.axvline(m_div_ti, color=COLOR_TI, linewidth=Mean_line_width, linestyle='-')
_style(ax, Title_div)
ax.set_xlabel(Xlabel_div, color='black', fontsize=Xlabel_size, fontweight=Xlabel_type)
ax.set_ylabel(Ylabel_div, color='black', fontsize=Ylabel_size, fontweight=Ylabel_type)
leg = ax.legend(facecolor='white', edgecolor='black', fontsize=legend_size, loc=legend_loc)
for text in leg.get_texts():
    text.set_color('black')
    text.set_fontweight(legend_type)

# rotor: una componente por panel
curl_components = [
    (curlx_zr, curlx_ti, m_curlx_zr, m_curlx_ti, bins_curlx, Title_curl_x, Xlabel_curl_x),
    (curly_zr, curly_ti, m_curly_zr, m_curly_ti, bins_curly, Title_curl_y, Xlabel_curl_y),
    (curlz_zr, curlz_ti, m_curlz_zr, m_curlz_ti, bins_curlz, Title_curl_z, Xlabel_curl_z),
]

for i, (c_zr, c_ti, m_zr, m_ti, bins_c, title, xlabel) in enumerate(curl_components, start=1):
    ax = axes[i]
    ax.hist(c_zr, bins=bins_c, color=COLOR_ZR, alpha=Hist_alpha,
            label=f'Zr (⟨⟩={m_zr:.2e})')
    ax.hist(c_ti, bins=bins_c, color=COLOR_TI, alpha=Hist_alpha,
            label=f'Ti (⟨⟩={m_ti:.2e})')
    ax.axvline(0, color='black', linewidth=Zero_line_width, linestyle=Zero_line_style)
    ax.axvline(m_zr, color=COLOR_ZR, linewidth=Mean_line_width, linestyle='-')
    ax.axvline(m_ti, color=COLOR_TI, linewidth=Mean_line_width, linestyle='-')
    _style(ax, title)
    ax.set_xlabel(xlabel, color='black', fontsize=Xlabel_size, fontweight=Xlabel_type)
    ax.set_ylabel(Ylabel_curl, color='black', fontsize=Ylabel_size, fontweight=Ylabel_type)
    leg = ax.legend(facecolor='white', edgecolor='black', fontsize=legend_size, loc=legend_loc)
    for text in leg.get_texts():
        text.set_color('black')
        text.set_fontweight(legend_type)

fig.suptitle(Suptitle, color='black', fontsize=Suptitle_size, fontweight=Suptitle_type)
fig.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig(Outfile, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"Guardado: {Outfile}")
