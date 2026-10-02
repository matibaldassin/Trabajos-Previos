import numpy as np
import matplotlib.pyplot as plt

# INPUT FILES
INIT_POL1='/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_2/pol-histeresis-0.40-S0'
INIT_POL2='/home/matias/Escritorio/CONICET/bzt/Ciclos_de_histeresis/Datos_histeresis_ciclo_2/pol-histeresis-0.40-S50'

TEMPORAL_FILE1='0.40-S0/P-total.txt'
TEMPORAL_FILE2='0.40-S50/P-total.txt'

# TEMPORAL STEP
dt = 0.0004   # picosegundos

# READING INICIAL POLARIZATION STATE
initial_polarization1 = np.loadtxt(INIT_POL1, skiprows=12, max_rows=1)
initial_polarization2 = np.loadtxt(INIT_POL2, skiprows=12, max_rows=1)

# READING DATA
data1 = np.loadtxt(TEMPORAL_FILE1, skiprows=2)
data2 = np.loadtxt(TEMPORAL_FILE2, skiprows=2)

time1 = data1[:, 0]*dt
polarization1 = data1[:, 3]
time2 = data2[:, 0]*dt
polarization2 = data2[:, 3]

# PLOTTING
plt.figure(figsize=(10, 6))
plt.plot(0, initial_polarization1[3], 'bo')
plt.plot(0, initial_polarization2[3], 'go')
plt.plot(time1, polarization1, color='blue', marker='o', label='S0')
plt.plot(time2, polarization2, color='green', marker='s', label='S50')
plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización (C/m²)')
plt.title('Relajación de la Polarización en el tiempo')
plt.legend()
plt.grid(True)
plt.show()


# NORMALIZING POLARIZATION
normalized_polarization1 = polarization1 / initial_polarization1[3]
normalized_polarization2 = polarization2 / initial_polarization2[3]

# EXPONENTIAL FITTING
from scipy.optimize import curve_fit
def exponential_decay(t, A, tau, B, C):
    return A * np.exp(-(t / tau)**B) + C

curve_fit1, _ = curve_fit(exponential_decay, time1, normalized_polarization1, p0=(1, 100, 1, 0), maxfev=10000)
curve_fit2, _ = curve_fit(exponential_decay, time2, normalized_polarization2, p0=(1, 100, 1, 0), maxfev=10000)

A1, tau1, B1, C1 = curve_fit1
A2, tau2, B2, C2 = curve_fit2

print(f"Fitting parameters for S0: A={A1:.3f}, tau={tau1:.3f} ps, B={B1:.3f}, C={C1:.3f}")
print(f"Fitting parameters for S50: A={A2:.3f}, tau={tau2:.3f} ps, B={B2:.3f}, C={C2:.3f}")

# PLOTTING NORMALIZED POLARIZATION
plt.figure(figsize=(10, 6))
plt.plot(0, 1, 'o', color='black')
plt.plot(time1, normalized_polarization1, color='blue', marker='o', label='S0')
plt.plot(time2, normalized_polarization2, color='green', marker='s', label='S50')
plt.plot(time1, exponential_decay(time1, *curve_fit1), color='blue', linestyle='--')
plt.plot(time2, exponential_decay(time2, *curve_fit2), color='green', linestyle='--')
plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización Normalizada')
plt.title('Relajación de la Polarización Normalizada en el tiempo')
plt.legend()
plt.grid(True)
plt.show()
