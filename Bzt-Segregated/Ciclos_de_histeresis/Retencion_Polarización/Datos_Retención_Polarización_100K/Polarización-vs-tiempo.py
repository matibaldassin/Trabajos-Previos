import numpy as np
import matplotlib.pyplot as plt

# INPUT FILES

TEMPORAL_FILE1='0.40-S0/P-total.txt'
TEMPORAL_FILE2='0.40-S50/P-total.txt'
TEMPORAL_FILE3='0.40-S100/P-total.txt'

# TEMPORAL STEP
dt = 0.0004   # picosegundos

# READING DATA
data1 = np.loadtxt(TEMPORAL_FILE1, skiprows=2)
data2 = np.loadtxt(TEMPORAL_FILE2, skiprows=2)
data3 = np.loadtxt(TEMPORAL_FILE3, skiprows=2)

time1 = data1[:, 0]*dt
time2 = data2[:, 0]*dt
time3 = data3[:, 0]*dt

polarization1 = data1[:, 3]
polarization2 = data2[:, 3]
polarization3 = data3[:, 3]

# PLOTTING
plt.figure(figsize=(10, 6))
plt.plot(time1, polarization1, color='#1f77b4', marker='o', label='S0')
plt.plot(time2, polarization2, color='#ff7f0e', marker='o', label='S50')
plt.plot(time3, polarization3, color='#2ca02c', marker='o', label='S100')

plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización (C/m²)')
plt.title('Relajación de la Polarización en el tiempo')
plt.legend()
plt.grid(True)
plt.show()


# NORMALIZING POLARIZATION
normalized_polarization1 = polarization1 / polarization1[0]
normalized_polarization2 = polarization2 / polarization2[0]
normalized_polarization3 = polarization3 / polarization3[0]

# EXPONENTIAL FITTING
from scipy.optimize import curve_fit
def exponential_decay(t, A, tau, B, C):
    return A * np.exp(-(t / tau)**B) + C

curve_fit1, _ = curve_fit(exponential_decay, time1, normalized_polarization1, p0=(1, 100, 1, 0), maxfev=10000)
curve_fit2, _ = curve_fit(exponential_decay, time2, normalized_polarization2, p0=(1, 100, 1, 0), maxfev=10000)
curve_fit3, _ = curve_fit(exponential_decay, time3, normalized_polarization3, p0=(1, 1500, 1, 0.80), maxfev=1000000)

A1, tau1, Beta1, Asint1 = curve_fit1
A2, tau2, Beta2, Asint2 = curve_fit2
A3, tau3, Beta3, Asint3 = curve_fit3

print(f"Fitting parameters for S0: A={A1:.3f}, tau={tau1:.3f} ps, Beta={Beta1:.3f}, Asíntota={Asint1:.3f}")
print(f"Fitting parameters for S50: A={A2:.3f}, tau={tau2:.3f} ps, Beta={Beta2:.3f}, Asíntota={Asint2:.3f}")
print(f"Fitting parameters for S100: A={A2:.3f}, tau={tau3:.3f} ps, Beta={Beta3:.3f}, Asíntota={Asint3:.3f}")

# PLOTTING NORMALIZED POLARIZATION
plt.figure(figsize=(10, 6))
plt.plot(time1, normalized_polarization1, color='#1f77b4', marker='o', label='S0')
plt.plot(time2, normalized_polarization2, color='#ff7f0e', marker='o', label='S50')
plt.plot(time3, normalized_polarization3, color='#2ca02c', marker='o', label='S100')

plt.plot(time1, exponential_decay(time1, *curve_fit1), color='black', linestyle='--')
plt.plot(time2, exponential_decay(time2, *curve_fit2), color='black', linestyle='--')
plt.plot(time3, exponential_decay(time3, *curve_fit3), color='black', linestyle='--')

plt.hlines(Asint1, time1[1], time1[-1], linestyles='--', color='#1f77b4', alpha=0.5)
plt.hlines(Asint2, time2[1], time2[-1], linestyles='--', color='#ff7f0e', alpha=0.5)
plt.hlines(Asint3, time3[1], time3[-1], linestyles='--', color='#2ca02c', alpha=0.5)

plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización Normalizada')
plt.title('Relajación de la Polarización Normalizada en el tiempo')
plt.legend()
plt.savefig('Retención-pol-0.40-S0.png', dpi=300, bbox_inches='tight')
plt.show()
