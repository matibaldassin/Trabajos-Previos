import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import csv

# INPUT POLARIZATION FILES
FILE1="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S0"
FILE2="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S50"
FILE3="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S100"
FILE4="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S0"
FILE5="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S50"
FILE6="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S100"

# INPUT VOLUME FILES
FILE11="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S0"
FILE12="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S50"
FILE13="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S100"
FILE14="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S0"
FILE15="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S50"
FILE16="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S100"

# OUTPUT FILES
OUT_FILE_1 = "Volumen-0.20.png"
OUT_FILE_2 = "Volumen-0.40.png"

# READING FUNCTIONS
def read_polarizations(f):
    lector = csv.reader(f,delimiter=' ')
    T,px,py,pz,p=[],[],[],[],[]
    for fila in lector:
        aux = [float(x) for x in fila if x]
        T.append(aux[0])
        px.append(aux[1])
        py.append(aux[2])
        pz.append(aux[3])
        p.append(aux[4])
    return T,px,py,pz,p

def read_volumes(f):
    lector = csv.reader(f,delimiter=' ')
    T,a,b,c,alpha,beta,gamma,vol=[],[],[],[],[],[],[],[]
    for fila in lector:
        aux = [float(x) for x in fila if x]
        T.append(aux[0])
        a.append(aux[1])
        b.append(aux[2])
        c.append(aux[3])
        alpha.append(aux[4])
        beta.append(aux[5])
        gamma.append(aux[6])
        vol.append(aux[7])
    return T,a,b,c,alpha,beta,gamma,vol

# PLOTTING SETTINGS
fig_size=[(8,6), (8,6)]

Titles=['a) x = 0.2','b) x = 0.4']
Titles_size=20
Titles_type='bold'

Xlabel=r'Temperature $[K]$'
Xlabel_size=20
Xlabel_type='normal'
Xlimit=(0,800)

Ylabel=[r'Unit Cell Volume $[\AA^3]$',r'Polarization $[\mu C/cm^2]$']
Ylabel_size=20
Ylabel_type='normal'
Ylimit1=[(65.7,67.0),(-1,18)]
Ylimit2=[(67.6,69.0),(-1,14)]

Marker_size=8
line_width=2.0

legend_size=16
legend_type='normal'
legend_loc1=['upper left','upper right']
legend_loc2=['upper left','upper right']

Grid_style='--'
Grid_alpha=0.8
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=20

# READING THE FILES
with open(FILE1) as f1, open(FILE2) as f2, open(FILE3) as f3, open(FILE4) as f4, open(FILE5) as f5, open(FILE6) as f6:
    T1,px1,py1,pz1,p1 = read_polarizations(f1)
    T2,px2,py2,pz2,p2 = read_polarizations(f2)
    T3,px3,py3,pz3,p3 = read_polarizations(f3)
    T4,px4,py4,pz4,p4 = read_polarizations(f4)
    T5,px5,py5,pz5,p5 = read_polarizations(f5)
    T6,px6,py6,pz6,p6 = read_polarizations(f6)
with open(FILE11) as f11, open(FILE12) as f12, open(FILE13) as f13, open(FILE14) as f14, open(FILE15) as f15, open(FILE16) as f16:
    T11,a11,b11,c11,alpha11,beta11,gamma11,vol11 = read_volumes(f11)
    T12,a12,b12,c12,alpha12,beta12,gamma12,vol12 = read_volumes(f12)
    T13,a13,b13,c13,alpha13,beta13,gamma13,vol13 = read_volumes(f13)
    T14,a14,b14,c14,alpha14,beta14,gamma14,vol14 = read_volumes(f14)
    T15,a15,b15,c15,alpha15,beta15,gamma15,vol15 = read_volumes(f15)
    T16,a16,b16,c16,alpha16,beta16,gamma16,vol16 = read_volumes(f16)

# SCALING THE TEMPERATURES
T1 = 1.48*np.array(T1)
T2 = 1.48*np.array(T2)
T3 = 1.48*np.array(T3)
T4 = 1.48*np.array(T4)
T5 = 1.48*np.array(T5)
T6 = 1.48*np.array(T6)
T11 = 1.48*np.array(T11)
T12 = 1.48*np.array(T12)
T13 = 1.48*np.array(T13)
T14 = 1.48*np.array(T14)
T15 = 1.48*np.array(T15)
T16 = 1.48*np.array(T16)

# ===== FIGURA 1: x = 0.20 - Volume =====
fig, ax = plt.subplots(figsize=fig_size[0])
ax.plot(T11, vol11, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax.plot(T12, vol12, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax.plot(T13, vol13, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax.set_title(Titles[0], fontsize=Titles_size, fontweight=Titles_type)
ax.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax.set_ylabel(Ylabel[0], fontsize=Ylabel_size)
ax.set_xlim(Xlimit)
ax.set_ylim(Ylimit1[0])
ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
ax.legend(fontsize=legend_size, loc=legend_loc1[0])
ax.minorticks_on()
for spine in ax.spines.values():
    spine.set_linewidth(Spine_width)
ax.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
plt.savefig(OUT_FILE_1, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()

# ===== FIGURA 2: x = 0.40 - Volume =====
fig, ax = plt.subplots(figsize=fig_size[1])
ax.plot(T14, vol14, label=r'$S0$', marker='o', color="#1f77b4", markersize=Marker_size, linewidth=line_width)
ax.plot(T15, vol15, label=r'$S50$', marker='o', color="#ff7f0e", markersize=Marker_size, linewidth=line_width)
ax.plot(T16, vol16, label=r'$S100$', marker='o', color="#2ca02c", markersize=Marker_size, linewidth=line_width)
ax.set_title(Titles[1], fontsize=Titles_size, fontweight=Titles_type)
ax.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax.set_ylabel(Ylabel[0], fontsize=Ylabel_size)
ax.set_xlim(Xlimit)
ax.set_ylim(Ylimit2[0])
ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
ax.legend(fontsize=legend_size, loc=legend_loc1[0])
ax.minorticks_on()
for spine in ax.spines.values():
    spine.set_linewidth(Spine_width)
ax.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
plt.savefig(OUT_FILE_2, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()

