import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# INPUT FILES
FILE1="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.20-S0"
FILE2="/home/matias/Escritorio/CONICET/bzt/Pol_Macro/pol-macro/heating/polmacro-0.40-S0"
FILE3="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.20-S0"
FILE4="/home/matias/Escritorio/CONICET/bzt/Parám_de_Red/Parámetros-de-red-DL-POLY/parametros-0.40-S0"
FILE5="/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.20-S0.dat" 
FILE6="/home/matias/Escritorio/CONICET/bzt/Polarización_local/Polarización_local_datos/heating/Polarización-local-promedio-0.40-S0.dat" 

# OUTPUT FILE
OUT_FILE1="Pol-y-param-0.20-S0.png"
OUT_FILE2="Pol-y-param-0.40-S0.png"

# PLOT SETTINGS
fig_size=[(8,10), (8,10)]

Titles=['x = 0.2','x = 0.4']
Titles_size=20
Titles_type='bold'

Xlabel=r'Temperature $[K]$'
Xlabel_size=20
Xlabel_type='normal'
Xlimit=(0,600)

Ylabel=[r'Polarization $[\mu C/cm^2]$',r'Lattice Parameters $[\AA]$']
Ylabel_size=20
Ylabel_type='normal'
Ylimit1=[(-1,9),(4.035,4.050)]
Ylimit2=[(-1,9),(4.072,4.090)]

Marker_size=7
line_width=1.5

legend_size=16
legend_type='normal'
legend_loc1=['lower left','upper left']
legend_loc2=['upper left','upper left']

Grid_style='--'
Grid_alpha=0.8
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0

# INSET SETTINGS
Inset_bounds1=[0.5, 0.35, 0.44, 0.6]  # [x0, y0, ancho, alto] en fracción de ejes
Inset_bounds2=[0.4, 0.3, 0.44, 0.6]  # [x0, y0, ancho, alto] en fracción de ejes
Inset_Xlimit=(0,600)
Inset_Ylimit1=(4,18.5)
Inset_Ylimit2=(5,12.5)
Inset_Xlabel=r'$T$ $[K]$'
Inset_Ylabel=r'$P_{L}$ $[\mu C/cm^2]$'
Inset_label_size=14
Inset_tick_size=11
Inset_tick_major_width=2.0
Inset_tick_major_length=7.0
Inset_tick_minor_width=1.5
Inset_tick_minor_length=3.0
Inset_spine_width=1.5
Inset_legend_size=11
Inset_legend_loc='upper right'
Inset_marker='o'
Inset_color='#333333'
Inset_grid_style='--'
Inset_grid_alpha=0.8
Inset_grid_width=1.5

def variables(f):
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

def read_dat(filepath):
    T_dat, P_dat = [], []
    with open(filepath,'r') as f:
        for line in f:
            if line.strip().startswith('#'):
                continue
            cols = line.split()
            T_dat.append(float(cols[0]))
            P_dat.append(float(cols[1]))
    return 1.48*np.array(T_dat), np.array(P_dat)

with open(FILE1,"r") as f:
    T1,px1,py1,pz1,p1=variables(f)
with open(FILE2,"r") as g:
    T2,px2,py2,pz2,p2=variables(g)

with open(FILE3,"r") as h:
    _,a_F,b_F,c_F,_=variables(h)
with open(FILE4,"r") as i:
    _,a_R,b_R,c_R,_=variables(i)

T1=1.48*np.array(T1)
T2=1.48*np.array(T2)

T5, P_loc1 = read_dat(FILE5)
T6, P_loc2 = read_dat(FILE6)

# ===== FIGURA 1: x = 0.20 =====
fig1, ax1 = plt.subplots(2, 1, sharex=True, figsize=fig_size[0])

# Subplot 1: polarizaciones
ax1[0].plot(T1, px1, label=r'$p_x$', marker='d',color="#5E3C99", markersize=Marker_size, linewidth=line_width)
ax1[0].plot(T1, py1, label=r'$p_y$', marker='^',color="#E66101", markersize=Marker_size, linewidth=line_width)
ax1[0].plot(T1, pz1, label=r'$p_z$', marker='s',color="#1B9E77", markersize=Marker_size, linewidth=line_width)
ax1[0].set_title(Titles[0],size=Titles_size, fontweight=Titles_type)
ax1[0].set_ylabel(Ylabel[0],size=Ylabel_size, fontweight=Ylabel_type)
ax1[0].set_ylim(Ylimit1[0])
ax1[0].set_xlim(Xlimit)
ax1[0].tick_params(axis='both', which='major', labelsize=20, width=Tick_major_width, length=Tick_major_length)
ax1[0].minorticks_on()
ax1[0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax1[0].spines.values():
    spine.set_linewidth(Spine_width)
ax1[0].legend(fontsize=legend_size, loc=legend_loc1[0])
# ax1[0].grid(True,linestyle=Grid_style,alpha=Grid_alpha,linewidth=Grid_width)

# Inset: P_loc vs T
inset1 = ax1[0].inset_axes(Inset_bounds1)
inset1.plot(T5, P_loc1, marker=Inset_marker, color=Inset_color, label=r'$P_{loc}$')
inset1.set_xlabel(Inset_Xlabel, fontsize=Inset_label_size)
inset1.set_ylabel(Inset_Ylabel, fontsize=Inset_label_size)
inset1.set_xlim(Inset_Xlimit)
inset1.set_ylim(Inset_Ylimit1)
inset1.tick_params(axis='both', which='major', labelsize=Inset_tick_size, width=Inset_tick_major_width, length=Inset_tick_major_length)
inset1.minorticks_on()
inset1.tick_params(axis='both', which='minor', width=Inset_tick_minor_width, length=Inset_tick_minor_length)
for spine in inset1.spines.values():
    spine.set_linewidth(Inset_spine_width)
# inset1.grid(True, linestyle=Inset_grid_style, alpha=Inset_grid_alpha, linewidth=Inset_grid_width)

# Subplot 2: parámetros
ax1[1].plot(T1, a_F, label=r'$a$', marker='d',color="#5E3C99", markersize=Marker_size, linewidth=line_width)
ax1[1].plot(T1, b_F, label=r'$b$', marker='^',color="#E66101", markersize=Marker_size, linewidth=line_width)
ax1[1].plot(T1, c_F, label=r'$c$', marker='s',color="#1B9E77", markersize=Marker_size, linewidth=line_width)
ax1[1].set_xlabel(Xlabel,size=Xlabel_size, fontweight=Xlabel_type)
ax1[1].set_ylabel(Ylabel[1],size=Ylabel_size, fontweight=Ylabel_type)
ax1[1].set_ylim(Ylimit1[1])
ax1[1].set_xlim(Xlimit)
ax1[1].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
ax1[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax1[1].tick_params(axis='both', which='major', labelsize=20, width=Tick_major_width, length=Tick_major_length)
ax1[1].minorticks_on()
ax1[1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax1[1].spines.values():
    spine.set_linewidth(Spine_width)
ax1[1].legend(fontsize=legend_size, loc=legend_loc1[1])
# ax1[1].grid(True,linestyle=Grid_style,alpha=Grid_alpha,linewidth=Grid_width)

fig1.tight_layout()
fig1.savefig(OUT_FILE1, dpi=300,bbox_inches='tight')


# ===== FIGURA 2: x = 0.40 =====
fig2, ax2 = plt.subplots(2, 1, sharex=True, figsize=fig_size[1])

# Subplot 1: polarizaciones
ax2[0].plot(T2, px2, label=r'$p_x$', marker='d', color="#5E3C99", markersize=Marker_size, linewidth=line_width)
ax2[0].plot(T2, py2, label=r'$p_y$', marker='^', color="#E66101", markersize=Marker_size, linewidth=line_width)
ax2[0].plot(T2, pz2, label=r'$p_z$', marker='s', color="#1B9E77", markersize=Marker_size, linewidth=line_width)
ax2[0].set_title(Titles[1],size=Titles_size, fontweight=Titles_type)
ax2[0].set_ylabel(Ylabel[0],size=Ylabel_size, fontweight=Ylabel_type)
ax2[0].tick_params(axis='both', which='major', labelsize=20, width=Tick_major_width, length=Tick_major_length)
ax2[0].minorticks_on()
ax2[0].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax2[0].spines.values():
    spine.set_linewidth(Spine_width)
ax2[0].set_ylim(Ylimit2[0])
ax2[0].set_xlim(Xlimit)
ax2[0].legend(fontsize=legend_size, loc=legend_loc2[0])
# ax2[0].grid(True,linestyle=Grid_style,alpha=Grid_alpha,linewidth=Grid_width)

# Inset: P_loc vs T
inset2 = ax2[0].inset_axes(Inset_bounds2)
inset2.plot(T6, P_loc2, marker=Inset_marker, color=Inset_color, label=r'$P_{loc}$')
inset2.set_xlabel(Inset_Xlabel, fontsize=Inset_label_size)
inset2.set_ylabel(Inset_Ylabel, fontsize=Inset_label_size)
inset2.set_xlim(Inset_Xlimit)
inset2.set_ylim(Inset_Ylimit2)
inset2.tick_params(axis='both', which='major', labelsize=Inset_tick_size, width=Inset_tick_major_width, length=Inset_tick_major_length)
inset2.minorticks_on()
inset2.tick_params(axis='both', which='minor', width=Inset_tick_minor_width, length=Inset_tick_minor_length)
for spine in inset2.spines.values():
    spine.set_linewidth(Inset_spine_width)
# inset2.grid(True, linestyle=Inset_grid_style, alpha=Inset_grid_alpha, linewidth=Inset_grid_width)

# Subplot 2: parámetros
ax2[1].plot(T2, a_R, label=r'$a$', marker='d', color="#5E3C99", markersize=Marker_size, linewidth=line_width)
ax2[1].plot(T2, b_R, label=r'$b$', marker='^', color="#E66101", markersize=Marker_size, linewidth=line_width)
ax2[1].plot(T2, c_R, label=r'$c$', marker='s', color="#1B9E77", markersize=Marker_size, linewidth=line_width)
ax2[1].set_xlabel(Xlabel,size=Xlabel_size, fontweight=Xlabel_type)
ax2[1].set_ylabel(Ylabel[1],size=Ylabel_size, fontweight=Ylabel_type)
ax2[1].set_ylim(Ylimit2[1])
ax2[1].set_xlim(Xlimit)
ax2[1].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
ax2[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax2[1].tick_params(axis='both', which='major', labelsize=20, width=Tick_major_width, length=Tick_major_length)
ax2[1].minorticks_on()
ax2[1].tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
for spine in ax2[1].spines.values():
    spine.set_linewidth(Spine_width)
ax2[1].legend(fontsize=legend_size, loc=legend_loc2[1])
# ax2[1].grid(True,linestyle=Grid_style,alpha=Grid_alpha,linewidth=Grid_width)

fig2.tight_layout()
fig2.savefig(OUT_FILE2, dpi=300,bbox_inches='tight')

plt.show()

