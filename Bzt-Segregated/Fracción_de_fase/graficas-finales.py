import numpy as np
import matplotlib.pyplot as plt
import csv

# INPUT FILES
FILE1='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.20_S0.txt'
FILE2='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.20_S50.txt'
FILE3='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.20_S100.txt'
FILE4='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.40_S0.txt'
FILE5='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.40_S50.txt'
FILE6='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Conteo_de_Fases/Fases_locales_0.40_S100.txt'

FILE11='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.20-S0.dat'
FILE22='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.20-S50.dat'
FILE33='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.20-S100.dat'
FILE44='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.40-S0.dat'
FILE55='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.40-S50.dat'
FILE66='/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/Datos_vecinos_Ti/vecinos-Ti-0.40-S100.dat'

# READING FUNCTIONS
def fases_transformadas(nombre_archivo):
    with open(nombre_archivo, 'r') as archivo_csv:
        lector_csv = csv.reader(archivo_csv, delimiter=',')
        next(lector_csv)  # Saltar la primera fila (encabezados)
        datos = list(lector_csv)

    # Convertir los datos a un array de NumPy
    datos_array = np.array(datos, dtype=float)

    # Separar las columnas en variables individuales
    T = datos_array[:, 0]  
    Tetra = datos_array[:, 1]  
    Rombo = datos_array[:, 2]
    Orto = datos_array[:, 3]
    Cubica = datos_array[:, 4]
    Zr = datos_array[:, 5]  
    
    return T, Rombo, Orto, Tetra, Cubica, Zr

def numeros_vecinos_ti(nombre_archivo):
    rows = []
    with open(nombre_archivo, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # Asegurarse de que hay al menos 4 columnas numéricas
            if len(parts) < 4:
                continue
            try:
                vals = [int(float(x)) for x in parts[:4]]
            except ValueError:
                # si alguna columna no es numérica, saltarla
                continue
            rows.append(vals)

    if not rows:
        return np.array([]), np.array([]), np.array([]), np.array([])

    data = np.array(rows, dtype=int)
    species_index = data[:, 0]
    neigh_same = data[:, 1]
    config_index = data[:, 2]
    neigh_B = data[:, 3]

    return species_index, neigh_same, config_index, neigh_B

# READING PHASES FRACTION
T1, Rombo1, Orto1, Tetra1, Cubica1, _ = fases_transformadas(FILE1)
T2, Rombo2, Orto2, Tetra2, Cubica2, _ = fases_transformadas(FILE2)
T3, Rombo3, Orto3, Tetra3, Cubica3, _ = fases_transformadas(FILE3)
T4, Rombo4, Orto4, Tetra4, Cubica4, _ = fases_transformadas(FILE4)
T5, Rombo5, Orto5, Tetra5, Cubica5, _ = fases_transformadas(FILE5)
T6, Rombo6, Orto6, Tetra6, Cubica6, _ = fases_transformadas(FILE6)

# REGOUPING DATA
T = np.array([T1, T2, T3, T4, T5, T6])
Rombo = np.array([Rombo1, Rombo2, Rombo3, Rombo4, Rombo5, Rombo6])
Orto = np.array([Orto1, Orto2, Orto3, Orto4, Orto5, Orto6])
Tetra = np.array([Tetra1, Tetra2, Tetra3, Tetra4, Tetra5, Tetra6])
Cubica = np.array([Cubica1, Cubica2, Cubica3, Cubica4, Cubica5, Cubica6])

# PHASE FRACTION CALCULATION    
T = 1.48*T  
N_Ti = np.array([Rombo[i] + Orto[i] + Tetra[i] + Cubica[i] for i in range(len(T))])
f_R = np.array([Rombo[i] / N_Ti[i] for i in range(len(T))])
f_O = np.array([Orto[i] / N_Ti[i] for i in range(len(T))])
f_T = np.array([Tetra[i] / N_Ti[i] for i in range(len(T))])
f_C = np.array([Cubica[i] / N_Ti[i] for i in range(len(T))])

# READING NEIGHBORS DATA
_, neigh_same_1, _, _ = numeros_vecinos_ti(FILE11)
_, neigh_same_2, _, _ = numeros_vecinos_ti(FILE22)
_, neigh_same_3, _, _ = numeros_vecinos_ti(FILE33)
_, neigh_same_4, _, _ = numeros_vecinos_ti(FILE44)
_, neigh_same_5, _, _ = numeros_vecinos_ti(FILE55)
_, neigh_same_6, _, _ = numeros_vecinos_ti(FILE66)








# PLOTTING SETTINGS
fig_size=(16,12)

Titles=['x = 0.20', 'x = 0.40']
Titles_size=24
Titles_type='bold'

Subtitles=['', '', '', '', '', '', '']
Subtitles_size=20
Subtitles_type='bold'

Xlabel=[r'Temperature $[K]$','Number of Ti Neighbors']
Xlabel_size=24
Xlabel_type='normal'
Xlimit=(0,800)
Xtick_step=200

Ylabel=[r'Ti Phase Fraction', r'$\%$ of Ti Atoms']
Ylabel_size=24
Ylabel_type='normal'
Ylimit=[(-0.05,1.1),(-1,80)]
Ytick_step=0.5

Marker_size=9
line_width=2.5

legend_size=14
legend_type='normal'
legend_loc=' right'

Grid_style='--'
Grid_alpha=0.0
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_labelsize=16
Hspace=0.25

# PLOTTING PHASE FRACTIONS AND HISTOGRAM OF TI NEIGHBORS FOR x = 0.20
fig1, axes1 = plt.subplots(3, 2, figsize=fig_size)

# fig1.suptitle('x = 0.20', fontsize=Titles_size, fontweight=Titles_type)

# Single Y-axis label for column 0 (left) and column 1 (right)
fig1.text(0.02, 0.5, Ylabel[0], va='center', ha='center', rotation='vertical',
          fontsize=Ylabel_size, fontweight=Ylabel_type)
fig1.text(0.98, 0.5, Ylabel[1], va='center', ha='center', rotation='vertical',
          fontsize=Ylabel_size, fontweight=Ylabel_type)

# --- Histogram data ---
bins    = np.arange(-0.5, 7.5, 1)
centers = np.arange(0, 7, 1)

counts_1, _ = np.histogram(neigh_same_1, bins=bins)
pct_1 = counts_1 / len(neigh_same_1) * 100.0

counts_2, _ = np.histogram(neigh_same_2, bins=bins)
pct_2 = counts_2 / len(neigh_same_2) * 100.0

counts_3, _ = np.histogram(neigh_same_3, bins=bins)
pct_3 = counts_3 / len(neigh_same_3) * 100.0

# -----------------------------------------------------------------------
# Row 0, Col 0 — Phase fractions S0
ax00 = axes1[0, 0]
ax00.plot(T[0], f_R[0], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax00.plot(T[0], f_O[0], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax00.plot(T[0], f_T[0], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax00.plot(T[0], f_C[0], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax00.set_title(Subtitles[0], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax00.set_xlim(Xlimit)
ax00.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax00.set_ylim(Ylimit[0])
ax00.set_yticks(np.arange(0, 1.1, Ytick_step))
ax00.tick_params(labelbottom=False)
ax00.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax00.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax00.minorticks_on()
ax00.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax00.spines.values():
    spine.set_linewidth(Spine_width)
ax00.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 0, Col 1 — Histogram S0
ax01 = axes1[0, 1]
ax01.bar(centers, pct_1, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_1):
    ax01.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax01.set_title(Subtitles[3], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax01.set_xlim(-0.5, 6.5)
ax01.set_ylim(Ylimit[1])
ax01.set_xticks(np.arange(0, 7, 1))
ax01.yaxis.set_label_position('right')
ax01.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax01.tick_params(labelbottom=False)
ax01.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax01.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax01.minorticks_on()
ax01.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax01.spines.values():
    spine.set_linewidth(Spine_width)

# -----------------------------------------------------------------------
# Row 1, Col 0 — Phase fractions S50
ax10 = axes1[1, 0]
ax10.plot(T[1], f_R[1], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax10.plot(T[1], f_O[1], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax10.plot(T[1], f_T[1], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax10.plot(T[1], f_C[1], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax10.set_title(Subtitles[1], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax10.set_xlim(Xlimit)
ax10.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax10.set_ylim(Ylimit[0])
ax10.set_yticks(np.arange(0, 1.1, Ytick_step))
ax10.tick_params(labelbottom=False)
ax10.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax10.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax10.minorticks_on()
ax10.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax10.spines.values():
    spine.set_linewidth(Spine_width)
# ax10.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 1, Col 1 — Histogram S50
ax11 = axes1[1, 1]
ax11.bar(centers, pct_2, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_2):
    ax11.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax11.set_title(Subtitles[4], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax11.set_xlim(-0.5, 6.5)
ax11.set_ylim(Ylimit[1])
ax11.set_xticks(np.arange(0, 7, 1))
ax11.yaxis.set_label_position('right')
ax11.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax11.tick_params(labelbottom=False)
ax11.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax11.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax11.minorticks_on()
ax11.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax11.spines.values():
    spine.set_linewidth(Spine_width)

# -----------------------------------------------------------------------
# Row 2, Col 0 — Phase fractions S100 (bottom row — x labels shown)
ax20 = axes1[2, 0]
ax20.plot(T[2], f_R[2], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax20.plot(T[2], f_O[2], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax20.plot(T[2], f_T[2], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax20.plot(T[2], f_C[2], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax20.set_title(Subtitles[2], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax20.set_xlim(Xlimit)
ax20.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax20.set_ylim(Ylimit[0])
ax20.set_yticks(np.arange(0, 1.1, Ytick_step))
ax20.set_xlabel(Xlabel[0], fontsize=Xlabel_size)
ax20.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax20.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax20.minorticks_on()
ax20.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax20.spines.values():
    spine.set_linewidth(Spine_width)
# ax20.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 2, Col 1 — Histogram S100 (bottom row — x labels shown)
ax21 = axes1[2, 1]
ax21.bar(centers, pct_3, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_3):
    ax21.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax21.set_title(Subtitles[5], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax21.set_xlim(-0.5, 6.5)
ax21.set_ylim(Ylimit[1])
ax21.set_xticks(np.arange(0, 7, 1))
ax21.set_xlabel(Xlabel[1], fontsize=Xlabel_size)
ax21.yaxis.set_label_position('right')
ax21.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax21.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax21.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax21.minorticks_on()
ax21.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax21.spines.values():
    spine.set_linewidth(Spine_width)

plt.tight_layout(rect=[0.06, 0.02, 0.94, 0.97])
plt.subplots_adjust(hspace=Hspace)

# Row labels centered above each row
for ax_l, ax_r, label in zip([ax00, ax10, ax20], [ax01, ax11, ax21], ['a) S0', 'b) S50', 'c) S100']):
    pos_l = ax_l.get_position()
    pos_r = ax_r.get_position()
    x_center = (pos_l.x0 + pos_r.x1) / 2
    fig1.text(x_center, pos_l.y1+0.01, label, ha='center', va='bottom',
              fontsize=Subtitles_size, fontweight=Subtitles_type)

plt.savefig('Fases-Histogramas-0.20.png', dpi=300, bbox_inches='tight')
plt.show()












# PLOTTING PHASE FRACTIONS AND HISTOGRAM OF TI NEIGHBORS FOR x = 0.40
fig2, axes2 = plt.subplots(3, 2, figsize=fig_size)

# fig2.suptitle('x = 0.40', fontsize=Titles_size, fontweight=Titles_type)

# Single Y-axis label for column 0 (left) and column 1 (right)
fig2.text(0.02, 0.5, Ylabel[0], va='center', ha='center', rotation='vertical',
          fontsize=Ylabel_size, fontweight=Ylabel_type)
fig2.text(0.98, 0.5, Ylabel[1], va='center', ha='center', rotation='vertical',
          fontsize=Ylabel_size, fontweight=Ylabel_type)

# --- Histogram data ---
bins    = np.arange(-0.5, 7.5, 1)
centers = np.arange(0, 7, 1)

counts_4, _ = np.histogram(neigh_same_4, bins=bins)
pct_4 = counts_4 / len(neigh_same_4) * 100.0

counts_5, _ = np.histogram(neigh_same_5, bins=bins)
pct_5 = counts_5 / len(neigh_same_5) * 100.0

counts_6, _ = np.histogram(neigh_same_6, bins=bins)
pct_6 = counts_6 / len(neigh_same_6) * 100.0

# -----------------------------------------------------------------------
# Row 0, Col 0 — Phase fractions S0
ax00 = axes2[0, 0]
ax00.plot(T[3], f_R[3], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax00.plot(T[3], f_O[3], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax00.plot(T[3], f_T[3], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax00.plot(T[3], f_C[3], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax00.set_title(Subtitles[0], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax00.set_xlim(Xlimit)
ax00.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax00.set_ylim(Ylimit[0])
ax00.set_yticks(np.arange(0, 1.1, Ytick_step))
ax00.tick_params(labelbottom=False)
ax00.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax00.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax00.minorticks_on()
ax00.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax00.spines.values():
    spine.set_linewidth(Spine_width)
ax00.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 0, Col 1 — Histogram S0
ax01 = axes2[0, 1]
ax01.bar(centers, pct_4, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_4):
    ax01.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax01.set_title(Subtitles[3], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax01.set_xlim(-0.5, 6.5)
ax01.set_ylim(Ylimit[1])
ax01.set_xticks(np.arange(0, 7, 1))
ax01.yaxis.set_label_position('right')
ax01.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax01.tick_params(labelbottom=False)
ax01.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax01.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax01.minorticks_on()
ax01.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax01.spines.values():
    spine.set_linewidth(Spine_width)

# -----------------------------------------------------------------------
# Row 1, Col 0 — Phase fractions S50
ax10 = axes2[1, 0]
ax10.plot(T[4], f_R[4], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax10.plot(T[4], f_O[4], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax10.plot(T[4], f_T[4], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax10.plot(T[4], f_C[4], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax10.set_title(Subtitles[1], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax10.set_xlim(Xlimit)
ax10.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax10.set_ylim(Ylimit[0])
ax10.set_yticks(np.arange(0, 1.1, Ytick_step))
ax10.tick_params(labelbottom=False)
ax10.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax10.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax10.minorticks_on()
ax10.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax10.spines.values():
    spine.set_linewidth(Spine_width)
# ax10.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 1, Col 1 — Histogram S50
ax11 = axes2[1, 1]
ax11.bar(centers, pct_5, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_5):
    ax11.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax11.set_title(Subtitles[4], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax11.set_xlim(-0.5, 6.5)
ax11.set_ylim(Ylimit[1])
ax11.set_xticks(np.arange(0, 7, 1))
ax11.yaxis.set_label_position('right')
ax11.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax11.tick_params(labelbottom=False)
ax11.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax11.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax11.minorticks_on()
ax11.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax11.spines.values():
    spine.set_linewidth(Spine_width)

# -----------------------------------------------------------------------
# Row 2, Col 0 — Phase fractions S100 (bottom row — x labels shown)
ax20 = axes2[2, 0]
ax20.plot(T[5], f_R[5], 's-', color='red', markersize=Marker_size, linewidth=line_width, label='R')
ax20.plot(T[5], f_O[5], '^-', color='blue', markersize=Marker_size, linewidth=line_width, label='O')
ax20.plot(T[5], f_T[5], 'o-', color='darkgreen', markersize=Marker_size, linewidth=line_width, label='T')
ax20.plot(T[5], f_C[5], 'd-', color='orange', markersize=Marker_size, linewidth=line_width, label='C')
ax20.set_title(Subtitles[2], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax20.set_xlim(Xlimit)
ax20.set_xticks(np.arange(Xlimit[0], Xlimit[1] + 1, Xtick_step))
ax20.set_ylim(Ylimit[0])
ax20.set_yticks(np.arange(0, 1.1, Ytick_step))
ax20.set_xlabel(Xlabel[0], fontsize=Xlabel_size)
ax20.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax20.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax20.minorticks_on()
ax20.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax20.spines.values():
    spine.set_linewidth(Spine_width)
# ax20.legend(fontsize=legend_size, loc=legend_loc.strip())

# Row 2, Col 1 — Histogram S100 (bottom row — x labels shown)
ax21 = axes2[2, 1]
ax21.bar(centers, pct_6, width=1.0, alpha=0.7, color='#E66101', edgecolor='black', align='center')
for x, pct in zip(centers, pct_6):
    ax21.text(x, pct + 0.5, f'{pct:.1f}%', ha='center', va='bottom', fontsize=16)
ax21.set_title(Subtitles[5], fontsize=Subtitles_size, fontweight=Subtitles_type)
ax21.set_xlim(-0.5, 6.5)
ax21.set_ylim(Ylimit[1])
ax21.set_xticks(np.arange(0, 7, 1))
ax21.set_xlabel(Xlabel[1], fontsize=Xlabel_size)
ax21.yaxis.set_label_position('right')
ax21.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
ax21.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_labelsize)
ax21.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
# ax21.minorticks_on()
ax21.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
for spine in ax21.spines.values():
    spine.set_linewidth(Spine_width)

plt.tight_layout(rect=[0.06, 0.02, 0.94, 0.97])
plt.subplots_adjust(hspace=Hspace)

# Row labels centered above each row
for ax_l, ax_r, label in zip([ax00, ax10, ax20], [ax01, ax11, ax21], ['a) S0', 'b) S50', 'c) S100']):
    pos_l = ax_l.get_position()
    pos_r = ax_r.get_position()
    x_center = (pos_l.x0 + pos_r.x1) / 2
    fig2.text(x_center, pos_l.y1+0.01, label, ha='center', va='bottom',
              fontsize=Subtitles_size, fontweight=Subtitles_type)

plt.savefig('Fases-Histogramas-0.40.png', dpi=300, bbox_inches='tight')
plt.show()


