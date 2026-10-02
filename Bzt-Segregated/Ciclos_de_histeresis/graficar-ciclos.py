import matplotlib.pyplot as plt
import numpy as np

# Input file 
FILE1 = '/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-0.20-S0'
FILE2 = '/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-0.20-S50'
FILE3 = '/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-0.20-S100'
# FILE4 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-{CONC[1]}-S0'
# FILE5 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-{CONC[1]}-S50'
# FILE6 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_1/pol-histeresis-{CONC[1]}-S100'

FILE44 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_2/pol-histeresis-0.40-S0'
FILE55 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_2/pol-histeresis-0.40-S50'
FILE66 = f'/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_2/pol-histeresis-0.40-S100'

# Output file
Salida_1 = f'histeresis-{0.40}.png'
Salida_2 = f'histeresis-{0.20}.png'
# Salida_3 = f'histeresis-comparativa-0.40-S0.png'
# Salida_4 = f'histeresis-comparativa-0.40-S50.png'
# Salida_5 = f'histeresis-comparativa-0.40-S100.png'


# Function to read data from file
def read_data(file_path):
    E,Px,Py,Pz,P = [],[],[],[],[]
    with open(file_path, 'r') as f:
        for line in f:
            data = line.strip().split()
            if len(data) == 5:
                E.append(float(data[0]))
                Px.append(float(data[1]))
                Py.append(float(data[2]))
                Pz.append(float(data[3]))
                P.append(float(data[4]))
            else:
                print(f"Skipping line due to insufficient data: {line}")
    return E, Px, Py, Pz, P

# Read data from files
E1, Px1, Py1, Pz1, P1 = read_data(FILE1)
E2, Px2, Py2, Pz2, P2 = read_data(FILE2)
E3, Px3, Py3, Pz3, P3 = read_data(FILE3)
E4, Px4, Py4, Pz4, P4 = read_data(FILE44)
E5, Px5, Py5, Pz5, P5 = read_data(FILE55)
E6, Px6, Py6, Pz6, P6 = read_data(FILE66)

# PLOT SETTINGS
fig_size=[(8,14), (8,14)]

Titles=['a) x = 0.2','b) x = 0.4', 'c)', 'd)', 'e)', 'f)']
Titles_size=20
Titles_type='bold'

Xlabel=r'Electric Field $(E)$ [kV/mm]'
Xlabel_size=20
Xlabel_type='normal'
Xlimit=(-260, 260)

Ylabel_macro=r'Polarization $(P_z)$ [µC/cm²]'
Ylabel_size=20
Ylabel_type='normal'

Ylimit=(-20, 20)

Marker_size=8
line_width=2.0

legend_size=16
legend_type='normal'
legend_loc1=['upper left','upper right']


Grid_style='--'
Grid_alpha=0
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=20

# Plotting Histeresis curves - subplot 2x1
Salida = 'histeresis.png'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Top subplot: x = 0.20
ax1.plot(E1, Pz1, marker='o', markersize=Marker_size, linewidth=line_width, label='S0')
ax1.plot(E2, Pz2, marker='o', markersize=Marker_size, linewidth=line_width, label='S50')
ax1.plot(E3, Pz3, marker='o', markersize=Marker_size, linewidth=line_width, label='S100')
ax1.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.2)
ax1.axvline(0, color='black', linewidth=1, linestyle='--', alpha=0.2)
ax1.set_title(Titles[0], fontsize=Titles_size, fontweight=Titles_type)
ax1.set_xlim(Xlimit)
ax1.set_ylim(Ylimit)
ax1.set_yticks([-20, -10, 0, 10, 20])
ax1.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax1.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax1.minorticks_on()
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax1.grid(True, linestyle=Grid_style, linewidth=Grid_width, alpha=Grid_alpha)
ax1.legend(fontsize=legend_size, loc=legend_loc1[0])
for spine in ax1.spines.values():
    spine.set_linewidth(Spine_width)

# Bottom subplot: x = 0.40
ax2.plot(E4, Pz4, marker='o', markersize=Marker_size, linewidth=line_width, label='S0')
ax2.plot(E5, Pz5, marker='o', markersize=Marker_size, linewidth=line_width, label='S50')
ax2.plot(E6, Pz6, marker='o', markersize=Marker_size, linewidth=line_width, label='S100')
ax2.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.2)
ax2.axvline(0, color='black', linewidth=1, linestyle='--', alpha=0.2)
ax2.set_title(Titles[1], fontsize=Titles_size, fontweight=Titles_type)
ax2.set_xlim(Xlimit)
ax2.set_ylim(Ylimit)
ax2.set_yticks([-20, -10, 0, 10, 20])
ax2.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax2.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax2.minorticks_on()
ax2.grid(True, linestyle=Grid_style, linewidth=Grid_width, alpha=Grid_alpha)
# ax2.legend(fontsize=legend_size, loc=legend_loc1[1])
for spine in ax2.spines.values():
    spine.set_linewidth(Spine_width)

# Shared x- and y-axis labels
fig.supxlabel(Xlabel, fontsize=Xlabel_size, fontstyle=Xlabel_type)
fig.supylabel(Ylabel_macro, fontsize=Ylabel_size, fontstyle=Ylabel_type)

plt.tight_layout()
plt.subplots_adjust(hspace=0.15)
plt.savefig(Salida, dpi=300, bbox_inches='tight')
plt.show()
plt.close()


