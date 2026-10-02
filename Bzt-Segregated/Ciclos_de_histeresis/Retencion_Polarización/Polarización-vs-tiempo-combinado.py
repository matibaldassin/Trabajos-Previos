import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

# INPUT FILES - Panel a) (100K)
FILE_A1 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_100K/0.40-S0/P-total.txt"
FILE_A2 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_100K/0.40-S50/P-total.txt"
FILE_A3 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_100K/0.40-S100/P-total.txt"

# INPUT FILES - Panel b) (220K / 300K)
FILE_B1 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_150K/0.40-S50/P-total.txt"
FILE_B2 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_150K/0.40-S100/P-total.txt"
FILE_B3 = "/home/matias/Escritorio/Trabajos_previos/Bzt-Segregated/Ciclos_de_histeresis/Retencion_Polarización/Datos_Retención_Polarización_200K/0.40-S100/P-total.txt"

# TIME STEP
dt = 0.0004  # Ps

# Read data from files
def read_data(file):
    data = np.loadtxt(file, skiprows=2)
    return np.array(data[:, 0]) * dt, np.array(data[:, 3])

time_a1, polarization_a1 = read_data(FILE_A1)
time_a2, polarization_a2 = read_data(FILE_A2)
time_a3, polarization_a3 = read_data(FILE_A3)

time_b1, polarization_b1 = read_data(FILE_B1)
time_b2, polarization_b2 = read_data(FILE_B2)
time_b3, polarization_b3 = read_data(FILE_B3)

# Normalize polarization
polarization_a1 /= np.max(np.abs(polarization_a1))
polarization_a2 /= np.max(np.abs(polarization_a2))
polarization_a3 /= np.max(np.abs(polarization_a3))

polarization_b1 /= np.max(np.abs(polarization_b1))
polarization_b2 /= np.max(np.abs(polarization_b2))
polarization_b3 /= np.max(np.abs(polarization_b3))

# Smoothing (moving average over 5 points, endpoints left untouched) — plotting only
def smooth5(y):
    y = np.asarray(y, dtype=float)
    smoothed = y.copy()
    smoothed[2:-2] = (y[:-4] + y[1:-3] + y[2:-2] + y[3:-1] + y[4:]) / 5
    return smoothed

polarization_a1_smooth = smooth5(polarization_a1)
polarization_a2_smooth = smooth5(polarization_a2)
polarization_a3_smooth = smooth5(polarization_a3)

polarization_b1_smooth = smooth5(polarization_b1)
polarization_b2_smooth = smooth5(polarization_b2)
polarization_b3_smooth = smooth5(polarization_b3)

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
t_max_a = max(time_a1[-1], time_a2[-1], time_a3[-1])
bounds_a = ([1e-6, 1e-6], [t_max_a * 10000, np.inf])
popt_a1, pcov_a1 = curve_fit(exponential_fit, time_a1, polarization_a1_smooth, p0=[1, 1], bounds=bounds_a, maxfev=10000)
popt_a2, pcov_a2 = curve_fit(exponential_fit, time_a2, polarization_a2_smooth, p0=[1, 1], bounds=bounds_a, maxfev=10000)
popt_a3, pcov_a3 = curve_fit(exponential_fit, time_a3, polarization_a3_smooth, p0=[1, 1], bounds=bounds_a, maxfev=10000)
tau_a1, beta_a1 = popt_a1
tau_a2, beta_a2 = popt_a2
tau_a3, beta_a3 = popt_a3
tau_a1_err, beta_a1_err = np.sqrt(np.diag(pcov_a1))
tau_a2_err, beta_a2_err = np.sqrt(np.diag(pcov_a2))
tau_a3_err, beta_a3_err = np.sqrt(np.diag(pcov_a3))

# Correlation coefficient between tau and Beta: rho = cov(tau,Beta) / (sigma_tau * sigma_Beta)
rho_a1 = pcov_a1[0, 1] / (tau_a1_err * beta_a1_err)
rho_a2 = pcov_a2[0, 1] / (tau_a2_err * beta_a2_err)
rho_a3 = pcov_a3[0, 1] / (tau_a3_err * beta_a3_err)

t_max_b = max(time_b1[-1], time_b2[-1], time_b3[-1])
bounds_b = ([1e-6, 1e-6], [t_max_b * 1000, np.inf])
popt_b1, pcov_b1 = curve_fit(exponential_fit, time_b1, polarization_b1_smooth, p0=[1, 1], bounds=bounds_b, maxfev=10000)
popt_b2, pcov_b2 = curve_fit(exponential_fit, time_b2, polarization_b2_smooth, p0=[1, 1], bounds=bounds_b, maxfev=10000)
popt_b3, pcov_b3 = curve_fit(exponential_fit, time_b3, polarization_b3_smooth, p0=[1, 1], bounds=bounds_b, maxfev=10000)
tau_b1, beta_b1 = popt_b1
tau_b2, beta_b2 = popt_b2
tau_b3, beta_b3 = popt_b3
tau_b1_err, beta_b1_err = np.sqrt(np.diag(pcov_b1))
tau_b2_err, beta_b2_err = np.sqrt(np.diag(pcov_b2))
tau_b3_err, beta_b3_err = np.sqrt(np.diag(pcov_b3))

rho_b1 = pcov_b1[0, 1] / (tau_b1_err * beta_b1_err)
rho_b2 = pcov_b2[0, 1] / (tau_b2_err * beta_b2_err)
rho_b3 = pcov_b3[0, 1] / (tau_b3_err * beta_b3_err)

# Saving fitting parameters (with 1-sigma uncertainties and tau-Beta correlation from the covariance matrix) to a text file
with open("fitting_parameters_fig_a.txt", "w") as f:
    f.write("Fitting parameters for 0.40-S0 (150K):\n")
    f.write(f"P0: {P0}, tau: {tau_a1} +/- {tau_a1_err}, Beta: {beta_a1} +/- {beta_a1_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_a1}\n\n")
    f.write("Fitting parameters for 0.40-S50 (150K):\n")
    f.write(f"P0: {P0}, tau: {tau_a2} +/- {tau_a2_err}, Beta: {beta_a2} +/- {beta_a2_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_a2}\n\n")
    f.write("Fitting parameters for 0.40-S100 (150K):\n")
    f.write(f"P0: {P0}, tau: {tau_a3} +/- {tau_a3_err}, Beta: {beta_a3} +/- {beta_a3_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_a3}\n")

with open("fitting_parameters_fig_b.txt", "w") as f:
    f.write("Fitting parameters for 0.40-S50 (220K):\n")
    f.write(f"P0: {P0}, tau: {tau_b1} +/- {tau_b1_err}, Beta: {beta_b1} +/- {beta_b1_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_b1}\n\n")
    f.write("Fitting parameters for 0.40-S100 (220K):\n")
    f.write(f"P0: {P0}, tau: {tau_b2} +/- {tau_b2_err}, Beta: {beta_b2} +/- {beta_b2_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_b2}\n\n")
    f.write("Fitting parameters for 0.40-S100 (300K):\n")
    f.write(f"P0: {P0}, tau: {tau_b3} +/- {tau_b3_err}, Beta: {beta_b3} +/- {beta_b3_err}, Asymptote: {Asymptote}, rho(tau,Beta): {rho_b3}\n")

# PLOT SETTINGS (style from graficar-ciclos.py)
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
fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Panel a)
ax_a.plot(time_a1, polarization_a1_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S0', color='#1f77b4')
ax_a.plot(time_a2, polarization_a2_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S50', color='#ff7f0e')
ax_a.plot(time_a3, polarization_a3_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S100', color='#2ca02c')

time_fit_a1 = np.linspace(time_a1[0], time_a1[-1], 1000)
time_fit_a2 = np.linspace(time_a2[0], time_a2[-1], 1000)
time_fit_a3 = np.linspace(time_a3[0], time_a3[-1], 1000)
ax_a.plot(time_fit_a1, exponential_fit(time_fit_a1, tau_a1, beta_a1), linewidth=line_width, linestyle='--', color='#1f77b4')
ax_a.plot(time_fit_a2, exponential_fit(time_fit_a2, tau_a2, beta_a2), linewidth=line_width, linestyle='--', color='#ff7f0e')
ax_a.plot(time_fit_a3, exponential_fit(time_fit_a3, tau_a3, beta_a3), linewidth=line_width, linestyle='--', color='#2ca02c')

ax_a.set_title('a)', fontsize=Title_size, fontweight=Title_type)
ax_a.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax_a.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax_a.minorticks_on()
ax_a.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_a.grid(True, linestyle=Grid_style, linewidth=Grid_width, alpha=Grid_alpha)
ax_a.legend(fontsize=legend_size, loc='center right')
for spine in ax_a.spines.values():
    spine.set_linewidth(Spine_width)

# Panel b)
ax_b.plot(time_b1, polarization_b1_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S50 - T = 220 K', color='#ff7f0e')
ax_b.plot(time_b2, polarization_b2_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S100 - T = 220 K', color='#2ca02c')
ax_b.plot(time_b3, polarization_b3_smooth, linestyle='-', marker='o', markersize=Marker_size, linewidth=line_width, label='S100 - T = 300 K', color='#005500')

time_fit_b1 = np.linspace(time_b1[0], time_b1[-1], 1000)
time_fit_b2 = np.linspace(time_b2[0], time_b2[-1], 1000)
time_fit_b3 = np.linspace(time_b3[0], time_b3[-1], 1000)
ax_b.plot(time_fit_b1, exponential_fit(time_fit_b1, tau_b1, beta_b1), linewidth=line_width, linestyle='--', color='#ff7f0e')
ax_b.plot(time_fit_b2, exponential_fit(time_fit_b2, tau_b2, beta_b2), linewidth=line_width, linestyle='--', color='#2ca02c')
ax_b.plot(time_fit_b3, exponential_fit(time_fit_b3, tau_b3, beta_b3), linewidth=line_width, linestyle='--', color='#005500')

ax_b.set_title('b)', fontsize=Title_size, fontweight=Title_type)
ax_b.yaxis.set_major_locator(MultipleLocator(0.2))
ax_b.tick_params(axis='both', which='major', width=Tick_major_width, length=Tick_major_length, labelsize=Tick_label_size)
ax_b.tick_params(axis='both', which='minor', width=Tick_minor_width, length=Tick_minor_length)
ax_b.minorticks_on()
ax_b.grid(True, linestyle=Grid_style, linewidth=Grid_width, alpha=Grid_alpha)
ax_b.legend(fontsize=legend_size, loc='upper right')
for spine in ax_b.spines.values():
    spine.set_linewidth(Spine_width)

fig.supxlabel(Xlabel, fontsize=Xlabel_size, fontstyle=Xlabel_type)
fig.supylabel(Ylabel, fontsize=Ylabel_size, fontstyle=Ylabel_type)

plt.tight_layout()
plt.subplots_adjust(hspace=0.09)
plt.savefig("polarization_relaxation_combinado.png", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

