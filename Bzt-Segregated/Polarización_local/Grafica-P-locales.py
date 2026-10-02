import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.lines import Line2D

# INPUT FILES - Macro polarization
FILE1="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S0"
FILE2="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S50"
FILE3="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S100"
FILE4="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S0"
FILE5="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S50"
FILE6="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S100"

# INPUT FILES - Local polarization average
FILE11 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.20-S0.dat"
FILE12 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.20-S50.dat"
FILE13 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.20-S100.dat"
FILE14 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.40-S0.dat"
FILE15 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.40-S50.dat"
FILE16= "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.40-S100.dat"

# INPUT FILES - Local polarization by site B
FILE21 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S0.dat"
FILE22 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S50.dat"
FILE23 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.20-S100.dat"
FILE24 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S0.dat"
FILE25 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S50.dat"
FILE26 = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Por-tipo-de-sitio-B/P-local-por-sitio-B-0.40-S100.dat"

# OUTPUT FILES
OUTPUT_IMAGE1 = "Polarizaciones-0.20.png"
OUTPUT_IMAGE2 = "Polarizaciones-0.40.png"

# READING FUNCTIONS
def read_P_macro(file):
    T, P_macro = [], []
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            aux = [float(x) for x in row if x]
            T.append(aux[0])
            P_macro.append(aux[4])
    return np.array(T), np.array(P_macro)

def read_P_local(file):
    T, P_local, sigma = [], [], []
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            T.append(float(row[0]))
            P_local.append(float(row[1]))
            sigma.append(float(row[2]))
    return np.array(T), np.array(P_local), np.array(sigma)

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

# PLOT SETTINGS
fig_size=[(8,14), (8,14)]

Titles=['a) x = 0.2','b) x = 0.4', 'c)', 'd)', 'e)', 'f)']
Titles_size=20
Titles_type='bold'

Xlabel=r'Temperature $[K]$'
Xlabel_size=20
Xlabel_type='normal'
Xlimit=(0,850)

Ylabel_macro=r'$P$ $[\mu C/cm^2]$'
Ylabel_local=r'$P_L$ $[\mu C/cm^2]$'
Ylabel_site=r'$P_L$ $[\mu C/cm^2]$'
Ylabel_size=20
Ylabel_type='normal'

Ylimit_macro=[(-1,18), (-1,18)]
Ylimit1=[(5,22),(0,26)]
Ylimit2=[(5,17),(0,26)]

Marker_size=8
line_width=2.0

legend_size=16
legend_type='normal'
legend_loc1=['upper right','upper right']
legend_loc2=['upper right','upper right']

Grid_style='--'
Grid_alpha=0.8
Grid_width=1.5

legend_color_handles = [
    Line2D([0], [0], color="#180081", linewidth=line_width, label='Ti'),
    Line2D([0], [0], color="#ff0077f4", linewidth=line_width, label='Zr'),
]
legend_marker_handles = [
    Line2D([0], [0], marker='^', color='black', linestyle='None', markersize=Marker_size, label='S0'),
    Line2D([0], [0], marker='d', color='black', linestyle='None', markersize=Marker_size, label='S50'),
    Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=Marker_size, label='S100'),
]

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=20

# READING FILES - Macro polarization
T_m1, P_m1 = read_P_macro(FILE1)
T_m2, P_m2 = read_P_macro(FILE2)
T_m3, P_m3 = read_P_macro(FILE3)
T_m4, P_m4 = read_P_macro(FILE4)
T_m5, P_m5 = read_P_macro(FILE5)
T_m6, P_m6 = read_P_macro(FILE6)

# READING FILES - Local polarization average
T1, P_local1, sigma1 = read_P_local(FILE11)
T2, P_local2, sigma2 = read_P_local(FILE12)
T3, P_local3, sigma3 = read_P_local(FILE13)
T4, P_local4, sigma4 = read_P_local(FILE14)
T5, P_local5, sigma5 = read_P_local(FILE15)
T6, P_local6, sigma6 = read_P_local(FILE16)

# READING FILES - Local polarization by site B
T11, P_local_Ti11, P_local_Zr11, sigma_Ti11, sigma_Zr11 = read_P_local_by_site(FILE21)
T12, P_local_Ti12, P_local_Zr12, sigma_Ti12, sigma_Zr12 = read_P_local_by_site(FILE22)
T13, P_local_Ti13, P_local_Zr13, sigma_Ti13, sigma_Zr13 = read_P_local_by_site(FILE23)
T14, P_local_Ti14, P_local_Zr14, sigma_Ti14, sigma_Zr14 = read_P_local_by_site(FILE24)
T15, P_local_Ti15, P_local_Zr15, sigma_Ti15, sigma_Zr15 = read_P_local_by_site(FILE25)
T16, P_local_Ti16, P_local_Zr16, sigma_Ti16, sigma_Zr16 = read_P_local_by_site(FILE26)

# SCALING TEMPERATURES
T_m1=T_m1*1.48; T_m2=T_m2*1.48; T_m3=T_m3*1.48
T_m4=T_m4*1.48; T_m5=T_m5*1.48; T_m6=T_m6*1.48
T1=T1*1.48; T2=T2*1.48; T3=T3*1.48
T4=T4*1.48; T5=T5*1.48; T6=T6*1.48
T11=T11*1.48; T12=T12*1.48; T13=T13*1.48
T14=T14*1.48; T15=T15*1.48; T16=T16*1.48

# ===== FIGURA 1: x = 0.20 =====
fig1, ax1 = plt.subplots(3, 1, sharex=True, figsize=fig_size[0])

# Subplot 1: Macro polarization
ax1[0].plot(T_m1, P_m1, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax1[0].plot(T_m2, P_m2, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax1[0].plot(T_m3, P_m3, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax1[0].set_title(Titles[0], size=Titles_size, fontweight=Titles_type)
ax1[0].set_ylabel(Ylabel_macro, size=Ylabel_size, fontweight=Ylabel_type)
if Ylimit_macro[0] is not None:
    ax1[0].set_ylim(Ylimit_macro[0])
ax1[0].set_xlim(Xlimit)
ax1[0].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax1[0].minorticks_on()
ax1[0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax1[0].spines.values():
    spine.set_linewidth(Spine_width)
ax1[0].legend(fontsize=legend_size, loc=legend_loc1[0])

# Subplot 2: Local polarization average
ax1[1].plot(T1, P_local1, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax1[1].plot(T2, P_local2, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax1[1].plot(T3, P_local3, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax1[1].set_title(Titles[2], size=Titles_size, fontweight=Titles_type)
ax1[1].set_ylabel(Ylabel_local, size=Ylabel_size, fontweight=Ylabel_type)
ax1[1].set_ylim(Ylimit1[0])
ax1[1].set_xlim(Xlimit)
ax1[1].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax1[1].minorticks_on()
ax1[1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax1[1].spines.values():
    spine.set_linewidth(Spine_width)
# ax1[1].legend(fontsize=legend_size, loc=legend_loc1[0])

# Subplot 3: Local polarization by site B
ax1[2].plot(T11, P_local_Ti11, label=r'$S0$', marker='^', color="#180081", markersize=Marker_size, linewidth=line_width)
ax1[2].plot(T11, P_local_Zr11, label=r'$S50$', marker='^', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax1[2].plot(T12, P_local_Ti12, label=r'$S50$', marker='d', color="#180081", markersize=Marker_size, linewidth=line_width)
ax1[2].plot(T12, P_local_Zr12, label=r'$S50$', marker='d', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax1[2].plot(T13, P_local_Ti13, label=r'$S100$', marker='o', color="#180081", markersize=Marker_size, linewidth=line_width)
ax1[2].plot(T13, P_local_Zr13, label=r'$S100$', marker='o', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax1[2].set_title(Titles[4], size=Titles_size, fontweight=Titles_type)
ax1[2].set_xlabel(Xlabel, size=Xlabel_size, fontweight=Xlabel_type)
ax1[2].set_ylabel(Ylabel_site, size=Ylabel_size, fontweight=Ylabel_type)
ax1[2].set_ylim(Ylimit1[1])
ax1[2].set_xlim(Xlimit)
ax1[2].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax1[2].minorticks_on()
ax1[2].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax1[2].spines.values():
    spine.set_linewidth(Spine_width)
leg1_color = ax1[2].legend(handles=legend_color_handles, fontsize=legend_size, loc=(0.6, 0.79))
ax1[2].add_artist(leg1_color)
ax1[2].legend(handles=legend_marker_handles, fontsize=legend_size, loc=legend_loc1[1])

fig1.tight_layout()
fig1.savefig(OUTPUT_IMAGE1, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()

# ===== FIGURA 2: x = 0.40 =====
fig2, ax2 = plt.subplots(3, 1, sharex=True, figsize=fig_size[1])

# Subplot 1: Macro polarization
ax2[0].plot(T_m4, P_m4, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax2[0].plot(T_m5, P_m5, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax2[0].plot(T_m6, P_m6, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax2[0].set_title(Titles[1], size=Titles_size, fontweight=Titles_type)
ax2[0].set_ylabel(Ylabel_macro, size=Ylabel_size, fontweight=Ylabel_type)
if Ylimit_macro[1] is not None:
    ax2[0].set_ylim(Ylimit_macro[1])
ax2[0].set_xlim(Xlimit)
ax2[0].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax2[0].minorticks_on()
ax2[0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax2[0].spines.values():
    spine.set_linewidth(Spine_width)
ax2[0].legend(fontsize=legend_size, loc=legend_loc1[0])

# Subplot 2: Local polarization average
ax2[1].plot(T4, P_local4, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax2[1].plot(T5, P_local5, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax2[1].plot(T6, P_local6, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax2[1].set_title(Titles[3], size=Titles_size, fontweight=Titles_type)
ax2[1].set_ylabel(Ylabel_local, size=Ylabel_size, fontweight=Ylabel_type)
ax2[1].set_ylim(Ylimit2[0])
ax2[1].set_xlim(Xlimit)
ax2[1].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax2[1].minorticks_on()
ax2[1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax2[1].spines.values():
    spine.set_linewidth(Spine_width)
# ax2[1].legend(fontsize=legend_size, loc=legend_loc1[0])

# Subplot 3: Local polarization by site B
ax2[2].plot(T14, P_local_Ti14, label=r'$S0$', marker='^', color="#180081", markersize=Marker_size, linewidth=line_width)
ax2[2].plot(T14, P_local_Zr14, label=r'$S50$', marker='^', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax2[2].plot(T15, P_local_Ti15, label=r'$S50$', marker='d', color="#180081", markersize=Marker_size, linewidth=line_width)
ax2[2].plot(T15, P_local_Zr15, label=r'$S50$', marker='d', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax2[2].plot(T16, P_local_Ti16, label=r'$S100$', marker='o', color="#180081", markersize=Marker_size, linewidth=line_width)
ax2[2].plot(T16, P_local_Zr16, label=r'$S100$', marker='o', color="#ff0077f4", markersize=Marker_size, linewidth=line_width)
ax2[2].set_title(Titles[5], size=Titles_size, fontweight=Titles_type)
ax2[2].set_xlabel(Xlabel, size=Xlabel_size, fontweight=Xlabel_type)
ax2[2].set_ylabel(Ylabel_site, size=Ylabel_size, fontweight=Ylabel_type)
ax2[2].set_ylim(Ylimit1[1])
ax2[2].set_xlim(Xlimit)
ax2[2].tick_params(axis='both', which='major', labelsize=Tick_label_size, width=Tick_major_width, length=Tick_major_length)
ax2[2].minorticks_on()
ax2[2].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax2[2].spines.values():
    spine.set_linewidth(Spine_width)
leg2_color = ax2[2].legend(handles=legend_color_handles, fontsize=legend_size, loc=(0.6, 0.79))
ax2[2].add_artist(leg2_color)
ax2[2].legend(handles=legend_marker_handles, fontsize=legend_size, loc=legend_loc2[1])

fig2.tight_layout()
fig2.savefig(OUTPUT_IMAGE2, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()
