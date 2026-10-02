import numpy as np
import matplotlib.pyplot as plt
import csv

# PLOT SETTINGS
fig_size=[(8,6), (8,6),(60,74)]

Titles=['c) Zr cells','d) Ti cells']
Titles_size=20
Titles_type='bold'

Xlabel=r'Unit cell volume $[\AA^3]$'
Xlabel_size=20
Xlabel_type='normal'
Xlimit=[(67,72.5),(63,70)]

Ylabel=[r'Polarization $[\mu C/cm^2]$',r'Lattice Parameters $[\AA]$']
Ylabel_size=20
Ylabel_type='normal'
Ylimit1=[(0,701), (0,801), (0,35)]

Marker_size=7
line_width=1.5

legend_size=16
legend_type='normal'
legend_loc=['upper left','upper right']


Grid_style='--'
Grid_alpha=0
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=16

# IMPORT FILES
FILE_1 = 'Resultados-0.40-S0/cell_lengths_pos-prom-20-0_kept.csv'
FILE_2 = 'Resultados-0.40-S50/cell_lengths_pos-prom-20-0_kept.csv'
FILE_3 = 'Resultados-0.40-S100/cell_lengths_pos-prom-20-0_kept.csv'

# OUTPUT FILE
OUTPUT_FILE = 'histograma_volumenes.png'

# READ CSV FILES
def read_csv(file):
    with open(file, 'r') as f:
        reader = csv.reader(f)
        next(reader) # Skip header
        data = [row for row in reader]
    return data

data_1 = read_csv(FILE_1)
data_2 = read_csv(FILE_2)
data_3 = read_csv(FILE_3)

# CALCULATE VOLUMES
a_x1 = np.array([float(data_1[i][3]) for i in range(len(data_1))])
a_y1 = np.array([float(data_1[i][4]) for i in range(len(data_1))])
a_z1 = np.array([float(data_1[i][5]) for i in range(len(data_1))])
volumes_1 = a_x1 * a_y1 * a_z1

a_x2 = np.array([float(data_2[i][3]) for i in range(len(data_2))])
a_y2 = np.array([float(data_2[i][4]) for i in range(len(data_2))])
a_z2 = np.array([float(data_2[i][5]) for i in range(len(data_2))])
volumes_2 = a_x2 * a_y2 * a_z2

a_x3 = np.array([float(data_3[i][3]) for i in range(len(data_3))])
a_y3 = np.array([float(data_3[i][4]) for i in range(len(data_3))])
a_z3 = np.array([float(data_3[i][5]) for i in range(len(data_3))])
volumes_3 = a_x3 * a_y3 * a_z3

# PLOT HISTOGRAM Zr Cells
fig, ax = plt.subplots(figsize=fig_size[0])
ax.hist(volumes_1[:8781], bins=70, alpha=0.5, label='S0')
ax.hist(volumes_2[:8781], bins=70, alpha=0.5, label='S50')
ax.hist(volumes_3[:8903], bins=70, alpha=0.5, label='S100')

ax.set_title(Titles[0], fontsize=Titles_size, fontweight=Titles_type)
ax.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax.set_ylabel('Counts', fontsize=Ylabel_size)
ax.set_xlim(Xlimit[0])
if Ylimit1[0] is not None:
    ax.set_ylim(Ylimit1[0])

ax.annotate('', xy=(69.68, 0.78), xytext=(69.68, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C0', lw=1.5), annotation_clip=False)
ax.annotate('', xy=(70.12, 0.78), xytext=(70.12, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C1', lw=1.5), annotation_clip=False)
ax.annotate('', xy=(70.69, 0.78), xytext=(70.69, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C2', lw=1.5), annotation_clip=False)

ax.legend(fontsize=legend_size, loc=legend_loc[0])
ax.grid(linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
ax.minorticks_on()

for spine in ax.spines.values():
    spine.set_linewidth(Spine_width)
ax.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)

plt.savefig('histograma_volumenes_Zr.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()



# PLOT HISTOGRAM Ti Cells
fig, ax = plt.subplots(figsize=fig_size[1])
ax.hist(volumes_1[8781:], bins=70, alpha=0.5, label='S0')
ax.hist(volumes_2[8781:], bins=70, alpha=0.5, label='S50')
ax.hist(volumes_3[8903:], bins=70, alpha=0.5, label='S100')

ax.set_title(Titles[1], fontsize=Titles_size, fontweight=Titles_type)
ax.set_xlabel(Xlabel, fontsize=Xlabel_size)
ax.set_ylabel('Counts', fontsize=Ylabel_size)
ax.set_xlim(Xlimit[1])
if Ylimit1[1] is not None:
    ax.set_ylim(Ylimit1[1])

ax.annotate('', xy=(66.39, 0.78), xytext=(66.39, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C0', lw=1.5), annotation_clip=False)
ax.annotate('', xy=(66.46, 0.78), xytext=(66.46, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C1', lw=1.5), annotation_clip=False)
ax.annotate('', xy=(66.57, 0.78), xytext=(66.57, 0.93),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='->', color='C2', lw=1.5), annotation_clip=False)

ax.legend(fontsize=legend_size, loc=legend_loc[1])
ax.grid(linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
ax.minorticks_on()

for spine in ax.spines.values():
    spine.set_linewidth(Spine_width)
ax.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)

plt.savefig('histograma_volumenes_Ti.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()



# VOLUMES VS POLARIZATION
FILE_POL_1 = 'Polarizaciones_locales/P-local-20-0.40-S0'
FILE_POL_2 = 'Polarizaciones_locales/P-local-20-0.40-S50'
FILE_POL_3 = 'Polarizaciones_locales/P-local-20-0.40-S100'

def read_polarization(file):
    return np.loadtxt(file)

pol_1 = read_polarization(FILE_POL_1)
pol_2 = read_polarization(FILE_POL_2)
pol_3 = read_polarization(FILE_POL_3)

polarization_1 = np.array([pol_1[i][4] for i in range(len(pol_1))])
polarization_2 = np.array([pol_2[i][4] for i in range(len(pol_2))])
polarization_3 = np.array([pol_3[i][4] for i in range(len(pol_3))])

# Ti cells only
polarization_1_Ti = polarization_1[8780:21953]
polarization_2_Ti = polarization_2[8780:21953]
polarization_3_Ti = polarization_3[8902:21953]
volumes_1_Ti = volumes_1[8780:21953]
volumes_2_Ti = volumes_2[8780:21953]
volumes_3_Ti = volumes_3[8902:21953]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 18), sharex=True, sharey=True)

# a) S0
ax1.scatter(volumes_1_Ti, polarization_1_Ti, alpha=0.5, s=Marker_size)
ax1.set_title('a) S0', fontsize=Titles_size, fontweight=Titles_type, loc='center')
ax1.minorticks_on()
for spine in ax1.spines.values():
    spine.set_linewidth(Spine_width)
ax1.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax1.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax1.tick_params(axis='y', which='both', labelleft=True)

# b) S50
ax2.scatter(volumes_2_Ti, polarization_2_Ti, alpha=0.5, s=Marker_size)
ax2.set_title('b) S50', fontsize=Titles_size, fontweight=Titles_type, loc='center')
ax2.minorticks_on()
for spine in ax2.spines.values():
    spine.set_linewidth(Spine_width)
ax2.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax2.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax2.tick_params(axis='y', which='both', labelleft=True)

# c) S100
ax3.scatter(volumes_3_Ti, polarization_3_Ti, alpha=0.5, s=Marker_size)
ax3.set_title('c) S100', fontsize=Titles_size, fontweight=Titles_type, loc='center')
ax3.minorticks_on()
for spine in ax3.spines.values():
    spine.set_linewidth(Spine_width)
ax3.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax3.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax3.set_xlabel(r'Unit cell volume $[\AA^3]$', fontsize=Xlabel_size)

fig.supylabel(r'Polarization magnitude $[\mu C/cm^2]$', fontsize=Ylabel_size)
plt.subplots_adjust(hspace=0.15)
plt.savefig('volumen_vs_polarizacion_Ti.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()
