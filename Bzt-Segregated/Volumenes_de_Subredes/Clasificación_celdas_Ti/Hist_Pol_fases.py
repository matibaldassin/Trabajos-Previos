import numpy as np
import matplotlib.pyplot as plt

# INPUT FILES
FILE1='Archivos_SOM_grouped/grouped-1-0.40-S0_extraido.txt'
FILE2='Archivos_SOM_grouped/grouped-1-0.40-S50_extraido.txt'
FILE3='Archivos_SOM_grouped/grouped-1-0.40-S100_extraido.txt'

# OUTPUT FILE
OUTPUT_FILE = 'Hist_Pol_fases.png'

# PLOT SETTINGS
fig_size=(18, 5)

Titles=['a) Tetragonal', 'b) Orthorhombic', 'c) Rhombohedral']
Titles_size=20
Titles_type='bold'

Xlabel=r'$P_L$ $[\mu C/cm^2]$'
Xlabel_size=20

Ylabel='Counts'
Ylabel_size=20

Xlimit=[(12,32), None, None]
Ylimit=[None, None, None]

legend_size=16
legend_loc=['upper left', 'upper right', 'upper right']

Grid_style='--'
Grid_alpha=0
Grid_width=1.5

Spine_width=1.5
Tick_major_width=2.0
Tick_major_length=7.0
Tick_minor_width=1.5
Tick_minor_length=3.0
Tick_label_size=16

# READING FILES
def read_file(filename):
    data = np.loadtxt(filename, skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2]  # DipoleMagnitude, Group, id

dipole1, group1, id1 = read_file(FILE1)
dipole2, group2, id2 = read_file(FILE2)
dipole3, group3, id3 = read_file(FILE3)

# SPLITTING DATA BY GROUP
def split_by_group(dipole, group):
    groups = np.unique(group)
    return {g: dipole[group == g] for g in groups}

data1 = split_by_group(dipole1, group1)
data2 = split_by_group(dipole2, group2)
data3 = split_by_group(dipole3, group3)

Tetra_1, Rombo_1, Orto_1 = data1[0], data1[1], data1[2]
Tetra_2, Rombo_2, Orto_2 = data2[0], data2[1], data2[2]
Tetra_3, Rombo_3, Orto_3 = data3[0], data3[1], data3[2]

# MEAN VALUES AND COUNTING
print("Mean Polarization Magnitude - Tetragonal Cells:")
print('S0: ' + str(np.mean(Tetra_1)) + r' $ ± $ ' + str(np.std(Tetra_1)) + str(' (N=' + str(len(Tetra_1)) + ')'))
print('S50: ' + str(np.mean(Tetra_2)) + r' $ ± $ ' + str(np.std(Tetra_2)) + str(' (N=' + str(len(Tetra_2)) + ')'))
print('S100: ' + str(np.mean(Tetra_3)) + r' $ ± $ ' + str(np.std(Tetra_3)) + str(' (N=' + str(len(Tetra_3)) + ')'))

print("Mean Polarization Magnitude - Rhombohedral Cells:")
print('S0: ' + str(np.mean(Rombo_1)) + r' $ ± $ ' + str(np.std(Rombo_1)) + str(' (N=' + str(len(Rombo_1)) + ')'))
print('S50: ' + str(np.mean(Rombo_2)) + r' $ ± $ ' + str(np.std(Rombo_2)) + str(' (N=' + str(len(Rombo_2)) + ')'))
print('S100: ' + str(np.mean(Rombo_3)) + r' $ ± $ ' + str(np.std(Rombo_3)) + str(' (N=' + str(len(Rombo_3)) + ')'))

print("Mean Polarization Magnitude - Orthorhombic Cells:")
print('S0: ' + str(np.mean(Orto_1)) + r' $ ± $ ' + str(np.std(Orto_1)) + str(' (N=' + str(len(Orto_1)) + ')'))
print('S50: ' + str(np.mean(Orto_2)) + r' $ ± $ ' + str(np.std(Orto_2)) + str(' (N=' + str(len(Orto_2)) + ')'))
print('S100: ' + str(np.mean(Orto_3)) + r' $ ± $ ' + str(np.std(Orto_3)) + str(' (N=' + str(len(Orto_3)) + ')'))

# PLOT HISTOGRAMS - 1 row x 3 columns
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=fig_size, sharey=True)

# a) Tetragonal
ax1.hist(Tetra_1, bins=60, alpha=0.5, label='S0')
ax1.hist(Tetra_2, bins=60, alpha=0.5, label='S50')
ax1.hist(Tetra_3, bins=60, alpha=0.5, label='S100')
ax1.set_title(Titles[0], fontsize=Titles_size, fontweight=Titles_type)
ax1.set_ylabel(Ylabel, fontsize=Ylabel_size)
if Xlimit[0] is not None:
    ax1.set_xlim(Xlimit[0])
if Ylimit[0] is not None:
    ax1.set_ylim(Ylimit[0])
ax1.legend(fontsize=legend_size, loc=legend_loc[0])
ax1.grid(linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
ax1.minorticks_on()
for spine in ax1.spines.values():
    spine.set_linewidth(Spine_width)
ax1.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax1.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)

# b) Orthorhombic
ax2.hist(Orto_1, bins=60, alpha=0.5, label='S0')
ax2.hist(Orto_2, bins=60, alpha=0.5, label='S50')
ax2.hist(Orto_3, bins=60, alpha=0.5, label='S100')
ax2.set_title(Titles[1], fontsize=Titles_size, fontweight=Titles_type)
ax2.set_xlabel(Xlabel, fontsize=Xlabel_size)
if Xlimit[1] is not None:
    ax2.set_xlim(Xlimit[1])
if Ylimit[1] is not None:
    ax2.set_ylim(Ylimit[1])
# ax2.legend(fontsize=legend_size, loc=legend_loc[1])
ax2.grid(linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
ax2.minorticks_on()
for spine in ax2.spines.values():
    spine.set_linewidth(Spine_width)
ax2.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax2.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)

# c) Rhombohedral
ax3.hist(Rombo_1, bins=60, alpha=0.5, label='S0')
ax3.hist(Rombo_2, bins=60, alpha=0.5, label='S50')
ax3.hist(Rombo_3, bins=60, alpha=0.5, label='S100')
ax3.set_title(Titles[2], fontsize=Titles_size, fontweight=Titles_type)
if Xlimit[2] is not None:
    ax3.set_xlim(Xlimit[2])
if Ylimit[2] is not None:
    ax3.set_ylim(Ylimit[2])
# ax3.legend(fontsize=legend_size, loc=legend_loc[2])
ax3.grid(linestyle=Grid_style, alpha=Grid_alpha, linewidth=Grid_width)
ax3.minorticks_on()
for spine in ax3.spines.values():
    spine.set_linewidth(Spine_width)
ax3.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax3.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)

plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.show()
plt.close()
