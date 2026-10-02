import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# INPUT FILES
FILE1 = "/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_150K/0.40-S50/P-total.txt"
FILE2 = "/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_150K/0.40-S100/P-total.txt"
FILE3 = "/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_200K/0.40-S100/P-total.txt"

# TIME STEP
dt = 0.0004  # Ps

# Read data from files
def read_data(file):
    data = np.loadtxt(file, skiprows=2)
    return np.array(data[:, 0]) * dt, np.array(data[:, 3])

time1, polarization1 = read_data(FILE1)
time2, polarization2 = read_data(FILE2)
time3, polarization3 = read_data(FILE3)

# Normalize polarization
polarization1 /= np.max(np.abs(polarization1))
polarization2 /= np.max(np.abs(polarization2))
polarization3 /= np.max(np.abs(polarization3))

# Exponential fitting function
# P0 is fixed at 1: the polarization is normalized, so the initial decay value is 1
# Asymptote is fixed at 0: full relaxation is assumed at long times
P0 = 1
Asymptote = 0
def exponential_fit(t, tau, Beta):
    return P0 * np.exp(-(t / tau) ** Beta) + Asymptote

# Fit exponential function to data
# tau and Beta are bounded to positive values to avoid negative-base powers (NaN)
# tau is also bounded to a reasonable order of the time range to avoid degenerate minima
t_max = max(time1[-1], time2[-1], time3[-1])
bounds = ([1e-6, 1e-6], [t_max * 1000, np.inf])
tau_1, beta_1 = curve_fit(exponential_fit, time1, polarization1, p0=[1, 1], bounds=bounds, maxfev=10000)[0]
tau_2, beta_2 = curve_fit(exponential_fit, time2, polarization2, p0=[1, 1], bounds=bounds, maxfev=10000)[0]
tau_3, beta_3 = curve_fit(exponential_fit, time3, polarization3, p0=[1, 1], bounds=bounds, maxfev=10000)[0]

# Saving fitting parameters to a text file
with open("fitting_parameters_fig_b.txt", "w") as f:
    f.write("Fitting parameters for 0.40-S50 (220K):\n")
    f.write(f"P0: {P0}, tau: {tau_1}, Beta: {beta_1}, Asymptote: {Asymptote}\n\n")
    f.write("Fitting parameters for 0.40-S100 (220K):\n")
    f.write(f"P0: {P0}, tau: {tau_2}, Beta: {beta_2}, Asymptote: {Asymptote}\n\n")
    f.write("Fitting parameters for 0.40-S100 (300K):\n")
    f.write(f"P0: {P0}, tau: {tau_3}, Beta: {beta_3}, Asymptote: {Asymptote}\n")

# PLOT SETTINGS (style from graficar-ciclos.py)
Title = 'b)'
Title_size = 20
Title_type = 'bold'

Xlabel = 'Time (ps)'
Xlabel_size = 20
Xlabel_type = 'normal'

Ylabel = r'$P/P_0$'
Ylabel_size = 20
Ylabel_type = 'normal'

Marker_size = 6
line_width = 2.0

legend_size = 16
legend_type = 'normal'
legend_loc = 'upper right'

Grid_style = '--'
Grid_alpha = 0
Grid_width = 1.5

Spine_width = 1.5
Tick_major_width = 2.0
Tick_major_length = 7.0
Tick_minor_width = 1.5
Tick_minor_length = 3.0
Tick_label_size = 20

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

line1, = ax.plot(time1, polarization1, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S50 - T = 220 K', color='#ff7f0e')
line2, = ax.plot(time2, polarization2, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S100 - T = 220 K', color='#2ca02c')
line3, = ax.plot(time3, polarization3, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S100 - T = 300 K', color='#005500')

time_fit1 = np.linspace(time1[0], time1[-1], 1000)
time_fit2 = np.linspace(time2[0], time2[-1], 1000)
time_fit3 = np.linspace(time3[0], time3[-1], 1000)
ax.plot(time_fit1, exponential_fit(time_fit1, tau_1, beta_1), linewidth=line_width, linestyle='--', color='#ff7f0e')
ax.plot(time_fit2, exponential_fit(time_fit2, tau_2, beta_2), linewidth=line_width, linestyle='--', color='#2ca02c')
ax.plot(time_fit3, exponential_fit(time_fit3, tau_3, beta_3), linewidth=line_width, linestyle='--', color='#005500')

ax.set_title(Title, fontsize=Title_size, fontweight=Title_type)
ax.set_xlabel(Xlabel, fontsize=Xlabel_size, fontstyle=Xlabel_type)
ax.set_ylabel(Ylabel, fontsize=Ylabel_size, fontstyle=Ylabel_type)
ax.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax.minorticks_on()
ax.grid(True, linestyle=Grid_style, linewidth=Grid_width, alpha=Grid_alpha)
ax.legend(fontsize=legend_size, loc=legend_loc)
for spine in ax.spines.values():
    spine.set_linewidth(Spine_width)

plt.tight_layout()
plt.savefig("polarization_relaxation_fig_b.png", dpi=300, bbox_inches='tight')
plt.show()
plt.close()