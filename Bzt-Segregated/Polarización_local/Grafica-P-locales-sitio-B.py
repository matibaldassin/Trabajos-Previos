import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.lines import Line2D

# INPUT FILES - Local polarization by site B
FILE21 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S0.dat"
FILE22 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S50.dat"
FILE23 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S100.dat"
FILE24 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S0.dat"
FILE25 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S50.dat"
FILE26 = "/home/usuario/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S100.dat"

# INPUT FILES - Lattice parameters (cell volume)
FILE31 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S0"
FILE32 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S50"
FILE33 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S100"
FILE34 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S0"
FILE35 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S50"
FILE36 = "/home/usuario/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S100"

# OUTPUT FILE
OUTPUT_IMAGE = "Polarizacion-local-sitio-B.png"

# READING FUNCTION
def read_P_local_by_site(file):
    T, P_local_Ti, sigma_Ti, P_local_Zr, sigma_Zr = [], [], [], [], []
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            T.append(float(row[0]))
            P_local_Ti.append(float(row[1]))
            P_local_Zr.append(float(row[2]))
            sigma_Ti.append(float(row[3]))
            sigma_Zr.append(float(row[4]))
    return np.array(T), np.array(P_local_Ti), np.array(P_local_Zr), np.array(sigma_Ti), np.array(sigma_Zr)

# Columns: T  a  b  c  alpha  beta  gamma  Volume (whitespace separated, no header)
def read_Volume(file):
    T, V = [], []
    with open(file, 'r') as f:
        for line in f:
            row = line.split()
            if not row:
                continue
            T.append(float(row[0]))
            V.append(float(row[7]))
    return np.array(T), np.array(V)

# PLOT SETTINGS
fig_size = (16, 10)

Titles = ['a) x = 0.2', 'b) x = 0.4']
Titles_size = 20
Titles_type = 'bold'

Xlabel = r'Temperature $[K]$'
Xlabel_size = 20
Xlabel_type = 'normal'
Xlimit = (0, 850)

Ylabel_site = r'$P_L$ $[\mu C/cm^2]$'
Ylabel_volume = r'$V$ $[\AA^3]$'
Ylabel_size = 20
Ylabel_type = 'normal'

Ylimit_site = [(0, 26), (0, 26)]
Ylimit_volume = [(65.5, 67.5), (67.5, 69.5)]

Marker_size = 8
line_width = 2.0

legend_size = 16
legend_loc = ['upper right', 'upper right']

Grid_style = '--'
Grid_alpha = 0.8
Grid_width = 1.5

legend_color_handles = [
    Line2D([0], [0], color="#180081", linewidth=line_width, label='Ti'),
    Line2D([0], [0], color="#ff0077f4", linewidth=line_width, label='Zr'),
]
legend_marker_handles = [
    Line2D([0], [0], marker='^', color='black', linestyle='None', markersize=Marker_size, label='S0'),
    Line2D([0], [0], marker='d', color='black', linestyle='None', markersize=Marker_size, label='S50'),
    Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=Marker_size, label='S100'),
]

Spine_width = 1.5
Tick_major_width = 2.0
Tick_major_length = 7.0
Tick_minor_width = 1.5
Tick_minor_length = 3.0
Tick_label_size = 20

# READING FILES - Local polarization by site B
T11, P_local_Ti11, P_local_Zr11, sigma_Ti11, sigma_Zr11 = read_P_local_by_site(FILE21)
T12, P_local_Ti12, P_local_Zr12, sigma_Ti12, sigma_Zr12 = read_P_local_by_site(FILE22)
T13, P_local_Ti13, P_local_Zr13, sigma_Ti13, sigma_Zr13 = read_P_local_by_site(FILE23)
T14, P_local_Ti14, P_local_Zr14, sigma_Ti14, sigma_Zr14 = read_P_local_by_site(FILE24)
T15, P_local_Ti15, P_local_Zr15, sigma_Ti15, sigma_Zr15 = read_P_local_by_site(FILE25)
T16, P_local_Ti16, P_local_Zr16, sigma_Ti16, sigma_Zr16 = read_P_local_by_site(FILE26)

# SCALING TEMPERATURES
T11 = T11 * 1.48; T12 = T12 * 1.48; T13 = T13 * 1.48
T14 = T14 * 1.48; T15 = T15 * 1.48; T16 = T16 * 1.48

# READING FILES - Cell volume
T31, V31 = read_Volume(FILE31)
T32, V32 = read_Volume(FILE32)
T33, V33 = read_Volume(FILE33)
T34, V34 = read_Volume(FILE34)
T35, V35 = read_Volume(FILE35)
T36, V36 = read_Volume(FILE36)

# SCALING TEMPERATURES
T31 = T31 * 1.48; T32 = T32 * 1.48; T33 = T33 * 1.48
T34 = T34 * 1.48; T35 = T35 * 1.48; T36 = T36 * 1.48

# ===== FIGURA: Local polarization by site B + Cell volume (x=0.20 and x=0.40) =====
fig, ax = plt.subplots(2, 2, sharex=True, figsize=fig_size)

# Subplot (0,0): P_local by site B, x = 0.20
ax[0,0].plot(T11, P_local_Ti11, label=r'$S0$', marker='^', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,0].plot(T11, P_local_Zr11, label=r'$S0$', marker='^', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,0].plot(T12, P_local_Ti12, label=r'$S50$', marker='d', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,0].plot(T12, P_local_Zr12, label=r'$S50$', marker='d', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,0].plot(T13, P_local_Ti13, label=r'$S100$', marker='o', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,0].plot(T13, P_local_Zr13, label=r'$S100$', marker='o', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,0].set_title(Titles[0], size=Titles_size, fontweight=Titles_type)
ax[0,0].set_ylabel(Ylabel_site, size=Ylabel_size, fontweight=Ylabel_type)
ax[0,0].set_ylim(Ylimit_site[0])
ax[0,0].set_xlim(Xlimit)
ax[0,0].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax[0,0].minorticks_on()
ax[0,0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax[0,0].spines.values():
    spine.set_linewidth(Spine_width)
leg0_color = ax[0,0].legend(handles=legend_color_handles, fontsize=legend_size, loc=(0.6, 0.76))
ax[0,0].add_artist(leg0_color)
ax[0,0].legend(handles=legend_marker_handles, fontsize=legend_size, loc=legend_loc[0])

# Subplot (0,1): P_local by site B, x = 0.40
ax[0,1].plot(T14, P_local_Ti14, label=r'$S0$', marker='^', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,1].plot(T14, P_local_Zr14, label=r'$S0$', marker='^', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,1].plot(T15, P_local_Ti15, label=r'$S50$', marker='d', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,1].plot(T15, P_local_Zr15, label=r'$S50$', marker='d', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,1].plot(T16, P_local_Ti16, label=r'$S100$', marker='o', color="#180081", markersize=Marker_size, linewidth=line_width)
ax[0,1].plot(T16, P_local_Zr16, label=r'$S100$', marker='o', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax[0,1].set_title(Titles[1], size=Titles_size, fontweight=Titles_type)
# ax[0,1].set_ylabel(Ylabel_site, size=Ylabel_size, fontweight=Ylabel_type)
ax[0,1].set_ylim(Ylimit_site[1])
ax[0,1].set_xlim(Xlimit)
ax[0,1].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax[0,1].minorticks_on()
ax[0,1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax[0,1].spines.values():
    spine.set_linewidth(Spine_width)
leg1_color = ax[0,1].legend(handles=legend_color_handles, fontsize=legend_size, loc=(0.6, 0.76))
ax[0,1].add_artist(leg1_color)
ax[0,1].legend(handles=legend_marker_handles, fontsize=legend_size, loc=legend_loc[1])

# Subplot (1,0): Cell volume, x = 0.20
ax[1,0].plot(T31, V31, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax[1,0].plot(T32, V32, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax[1,0].plot(T33, V33, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax[1,0].set_xlabel(Xlabel, size=Xlabel_size, fontweight=Xlabel_type)
ax[1,0].set_ylabel(Ylabel_volume, size=Ylabel_size, fontweight=Ylabel_type)
ax[1,0].set_ylim(Ylimit_volume[0])
ax[1,0].set_xlim(Xlimit)
ax[1,0].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax[1,0].minorticks_on()
ax[1,0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax[1,0].spines.values():
    spine.set_linewidth(Spine_width)
ax[1,0].legend(fontsize=legend_size, loc='upper left')

# Subplot (1,1): Cell volume, x = 0.40
ax[1,1].plot(T34, V34, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax[1,1].plot(T35, V35, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax[1,1].plot(T36, V36, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax[1,1].set_xlabel(Xlabel, size=Xlabel_size, fontweight=Xlabel_type)
# ax[1,1].set_ylabel(Ylabel_volume, size=Ylabel_size, fontweight=Ylabel_type)
ax[1,1].set_ylim(Ylimit_volume[1])
ax[1,1].set_xlim(Xlimit)
ax[1,1].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax[1,1].minorticks_on()
ax[1,1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax[1,1].spines.values():
    spine.set_linewidth(Spine_width)
ax[1,1].legend(fontsize=legend_size, loc='upper left')

fig.tight_layout()
fig.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()
