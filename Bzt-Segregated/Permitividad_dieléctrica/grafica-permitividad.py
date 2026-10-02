import csv
import matplotlib.pyplot as plt
import numpy as np

# INPUT FILES
FILE1 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.20-S0.dat"
FILE2 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.20-S50.dat"
FILE3 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.20-S100.dat"
FILE4 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.40-S0.dat"
FILE5 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.40-S50.dat"
FILE6 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Permitividad_dieléctrica/Datos_permitividad/permittivity-e111-table-0.40-S100.dat"

def variables(f):
    lector = csv.reader(f,delimiter=' ')
    next(lector)  # Salta la primera línea
    T,e_rel=[],[]
    for fila in lector:
        aux = [float(x) for x in fila if x]
        T.append(aux[0])
        e_rel.append(aux[1])
    return np.array(T),np.array(e_rel)

with open(FILE1, 'r') as f1, open(FILE2, 'r') as f2, open(FILE3, 'r') as f3, open(FILE4, 'r') as f4, open(FILE5, 'r') as f5, open(FILE6, 'r') as f6:
    T1,eps_1 = variables(f1)
    T2,eps_2 = variables(f2)
    T3,eps_3 = variables(f3)
    T4,eps_4 = variables(f4)
    T5,eps_5 = variables(f5)
    T6,eps_6 = variables(f6)

T1 =1.48*T1
T2 =1.48*T2
T3 =1.48*T3
T4 =1.48*T4
T5 =1.48*T5
T6 =1.48*T6

# PLOT SETTINGS
fig_size=(16,5)

Titles=['a) x = 0.2','b) x = 0.4']
Titles_size=20
Titles_type='bold'

Xlabel=r'Temperature $[K]$'
Xlabel_size=18
Xlabel_type='normal'
Xlimit=(0,850)

Ylabel=[r'$\varepsilon_r$','']
Ylabel_size=24
Ylabel_type='normal'
Ylimit=[None,None]

Marker_size=9
line_width=2.5

legend_size=14
legend_type='normal'
legend_loc=['lower center',(0.2,0.05)]


Grid_style='--'
Grid_alpha=0.8
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0

# # Comparing both concentrations for S0 agregation
# plt.figure(figsize=(10,6))

# plt.xlabel(Xlabel, size=Xlabel_size)
# plt.ylabel(Ylabel, size=Ylabel_size)
# plt.xlim(Xlimit)
# plt.ylim(Ylimit[0])
# plt.minorticks_on()
# plt.tick_params(axis='both', which='major', labelsize=18)

# plt.plot(T1,e1,label=r"$Ferro$", marker='o', color="#CE1313")
# plt.plot(T4,e4,label=r"$Relaxor$", marker='o', color="#006824")

# plt.grid(True,linestyle='--', alpha=0.5)
# plt.legend(fontsize=14)
# plt.savefig('Permitividad-comparativa-S0.png', dpi=300, bbox_inches='tight')
# plt.show()

# Ploting permittivity vs temperature for both concentrations and all agregations
fig, axes = plt.subplots(1, 2, figsize=fig_size)

data = [
    [(T1, eps_1, r"$S0$"), (T2, eps_2, r"$S50$"), (T3[:30], eps_3[:30], r"$S100$")],
    [(T4, eps_4, r"$S0$"), (T5, eps_5, r"$S50$"), (T6[:30], eps_6[:30], r"$S100$")],
]

for i, ax in enumerate(axes):
    for T, eps, label in data[i]:
        ax.plot(T, eps, label=label, marker='o', markersize=Marker_size, linewidth=line_width)

    ax.set_title(Titles[i], fontsize=Titles_size, fontweight=Titles_type)
    ax.set_xlabel(Xlabel, fontsize=Xlabel_size, fontweight=Xlabel_type)
    ax.set_ylabel(Ylabel[i], fontsize=Ylabel_size, fontweight=Ylabel_type)
    ax.set_xticks(np.arange(0, 850, 200))
    ax.set_xlim(Xlimit)
    ax.set_ylim(Ylimit[i])
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=18,
                   width=Tick_major_width, length=Tick_major_length)
    ax.tick_params(axis='both', which='minor',
                   width=Tick_minor_width, length=Tick_minor_length)
    # ax.grid(True, linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
    ax.legend(loc=legend_loc[i],
              prop={'weight': legend_type, 'size': legend_size})
    for spine in ax.spines.values():
        spine.set_linewidth(Spine_width)

plt.tight_layout()
plt.savefig('Permitividad.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

