import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.optimize as opt

CONC="0.20"

NOMBRE_ARCHIVO1 = f"{CONC}-S0/salida.csv"
NOMBRE_ARCHIVO2 = f"{CONC}-S50/salida.csv"
NOMBRE_ARCHIVO3 = f"{CONC}-S100/salida.csv"

grafica = f"Tangente-Hiperbolica-{CONC}.png"

def fases_transformadas(nombre_archivo):
    with open(nombre_archivo, 'r') as archivo_csv:
        lector_csv = csv.reader(archivo_csv, delimiter=',')
        next(lector_csv)  # Saltar la primera fila (encabezados)
        datos = list(lector_csv)

    # Convertir los datos a un array de NumPy
    datos_array = np.array(datos, dtype=float)

    # Separar las columnas en variables individuales
    T = datos_array[:, 0]  
    Tetra = datos_array[:, 1]  
    Rombo = datos_array[:, 2]
    Orto = datos_array[:, 3]
    Cubica = datos_array[:, 4]
    Zr = datos_array[:, 5]  
    
    return T, Rombo, Orto, Tetra, Cubica, Zr

T1,Rombo1,Orto1,Tetra1,Cubica1,Zr1 = fases_transformadas(NOMBRE_ARCHIVO1)
T2,Rombo2,Orto2,Tetra2,Cubica2,Zr2 = fases_transformadas(NOMBRE_ARCHIVO2)
T3,Rombo3,Orto3,Tetra3,Cubica3,Zr3 = fases_transformadas(NOMBRE_ARCHIVO3)

T1=1.45*T1
T2=1.45*T2
T3=1.45*T3

N_Ti1=Rombo1 + Orto1 + Tetra1 + Cubica1
N_Ti2=Rombo2 + Orto2 + Tetra2 + Cubica2
N_Ti3=Rombo3 + Orto3 + Tetra3 + Cubica3
f_C1=Cubica1/N_Ti1
f_C2=Cubica2/N_Ti2
f_C3=Cubica3/N_Ti3

# Ajuste de la tangente hiperbólica
def tangente_hiperbolica(T, T0, DT):
    return 0.5 * (1 + np.tanh((T - T0) / DT))

p_0 = [300, 10]  # Valores iniciales para T0 y DT
def ajustar_tangente_hiperbolica(T, f_C):
    popt, _ = opt.curve_fit(tangente_hiperbolica, T, f_C, p0=p_0)
    return popt

T0_1, DT_1 = ajustar_tangente_hiperbolica(T1, f_C1)
T0_2, DT_2 = ajustar_tangente_hiperbolica(T2, f_C2)
T0_3, DT_3 = ajustar_tangente_hiperbolica(T3, f_C3)

# Graficar los resultados
plt.figure(figsize=(8, 6))

x=np.linspace(min(T3), max(T3), 1000)

plt.xlabel(r'Temperatura $[K]$', fontsize=18)
plt.ylabel(r'Fracción de fase cúbica', fontsize=18)
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=18)
plt.grid(True, linestyle='--', alpha=0.5)

plt.plot(T1, f_C1, label='S0', marker='d',color="#1f77b4", markersize=8)
plt.plot(x, tangente_hiperbolica(x, T0_1, DT_1), color="#1f77b4", linewidth=2, linestyle='--')

plt.plot(T2, f_C2, label='S50', marker='d', color="#ff7f0e", markersize=8)
plt.plot(x, tangente_hiperbolica(x, T0_2, DT_2), color="#ff7f0e", linewidth=2, linestyle='--')

plt.plot(T3, f_C3, label='S100', marker='d', color="#2ca02c", markersize=8)
plt.plot(x, tangente_hiperbolica(x, T0_3, DT_3), color="#2ca02c", linewidth=2, linestyle='--')

plt.legend(fontsize=14)

plt.savefig(grafica, dpi=300, bbox_inches='tight')
plt.show()
plt.close()

with open(f"Tangente-Hiperbolica-{CONC}.txt", "w") as f:
    f.write(f"S0: T0 = {T0_1:.2f} K, DT = {DT_1:.2f} K\n")
    f.write(f"S50: T0 = {T0_2:.2f} K, DT = {DT_2:.2f} K\n")
    f.write(f"S100: T0 = {T0_3:.2f} K, DT = {DT_3:.2f} K\n")  
