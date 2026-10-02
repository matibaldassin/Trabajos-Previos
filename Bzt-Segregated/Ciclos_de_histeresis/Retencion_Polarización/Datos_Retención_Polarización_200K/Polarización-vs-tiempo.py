import numpy as np
import matplotlib.pyplot as plt

# INPUT FILES

TEMPORAL_FILE3='0.40-S100/P-total.txt'

# TEMPORAL STEP
dt = 0.0004   # picosegundos

# READING DATA
data3 = np.loadtxt(TEMPORAL_FILE3, skiprows=2)

time3 = data3[:, 0]*dt

polarization3 = data3[:, 3]

# PLOTTING
plt.figure(figsize=(10, 6))
plt.plot(time3, polarization3, color='#2ca02c', marker='o', label='S100')

plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización (C/m²)')
plt.title('Relajación de la Polarización en el tiempo')
plt.legend()
plt.grid(True)
plt.show()


# NORMALIZING POLARIZATION
normalized_polarization3 = polarization3 / polarization3[0]

# EXPONENTIAL FITTING
from scipy.optimize import curve_fit
def exponential_decay(t, A, tau, B, C):
    return A * np.exp(-(t / tau)**B) + C

curve_fit3, _ = curve_fit(exponential_decay, time3, normalized_polarization3, p0=(1, 1, 1, 0), maxfev=10000)

A3, tau3, Beta3, Asint3 = curve_fit3

print(f"Fitting parameters for S100: A={A3:.3f}, tau={tau3:.3f} ps, Beta={Beta3:.3f}, Asíntota={Asint3:.3f}")

# PLOTTING NORMALIZED POLARIZATION
plt.figure(figsize=(10, 6))
plt.plot(time3, normalized_polarization3, color='#2ca02c', marker='o', label='S100')

plt.plot(time3, exponential_decay(time3, *curve_fit3), color='black', linestyle='--')

plt.hlines(Asint3, time3[1], time3[-1], linestyles='--', color='#2ca02c', alpha=0.5)

plt.xlabel('Tiempo (ps)')
plt.ylabel('Polarización Normalizada')
plt.title('Relajación de la Polarización Normalizada en el tiempo')
plt.legend()
plt.savefig('Retención-pol-0.40-S0.png', dpi=300, bbox_inches='tight')
plt.show()
