import numpy as np
import matplotlib.pyplot as plt
import csv
import glob
import re

IN_FOLDER = "/home/matias/Escritorio/Tesina/graficos/Edward-Anderson/P-locales-prom-temp/heating/P-local-0.40-S0"

OUTPUT_FILE = "Edward-Anderson-Parameter-0.40-S0.dat"
OUTPUT_IMAGE = "Edward-Anderson-Parameter-0.40-S0.png"

def Edward_Anderson(f):
    lector = csv.reader(f, delimiter=' ')
    P_Ti = []
    P_Zr = []

    for fila in lector:
        aux = [float(x) for x in fila if x]

        if len(aux) <= 4:
            continue  # evita errores

        if aux[0] <= 8790:
            P_Zr.append(aux[4])
        else:
            P_Ti.append(aux[4])

    P_Ti = np.array(P_Ti)
    P_Zr = np.array(P_Zr)

    mean_Ti = np.mean(P_Ti**2)
    mean_Zr = np.mean(P_Zr**2)

    std_Ti = np.std(P_Ti**2)
    std_Zr = np.std(P_Zr**2)

    P_EA_Ti = np.sqrt(mean_Ti)
    P_EA_Zr = np.sqrt(mean_Zr)

    sigma_EA_Ti = std_Ti / (2 * P_EA_Ti)
    sigma_EA_Zr = std_Zr / (2 * P_EA_Zr)

    return P_EA_Ti, sigma_EA_Ti, P_EA_Zr, sigma_EA_Zr


def main():
    temperatures = list(range(10,401,10)) + list(range(410,791,20))
    T, EA_Ti, EA_Zr, sigma_EA_Ti, sigma_EA_Zr = [],[],[],[],[]

    for i in temperatures:
        with open(f"{IN_FOLDER}/P-local-{i}", "r", encoding="utf-8") as s:
            P_Ti, sigma_Ti, P_Zr, sigma_Zr = Edward_Anderson(s)

            T.append(1.48*i)
            EA_Ti.append(P_Ti)
            EA_Zr.append(P_Zr)
            sigma_EA_Ti.append(sigma_Ti)
            sigma_EA_Zr.append(sigma_Zr)

    with open(OUTPUT_FILE,"w") as r:
        r.write("# Temperature [K]\tEdward-Anderson Parameter (Ti)\tSigma (Ti)\tEdward-Anderson Parameter (Zr)\tSigma (Zr)\n")
        for i in range(len(T)):
            r.write(f"{T[i]} {EA_Ti[i]} {sigma_EA_Ti[i]} {EA_Zr[i]} {sigma_EA_Zr[i]}\n")

    plt.figure(figsize=(10,6))
    plt.title('Edward-Anderson Parameter vs Temperature')
    plt.xlabel(r"$Temperature \, [K]$")
    plt.ylabel(r"Edward-Anderson Parameter")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.minorticks_on()
    
    plt.errorbar(T, EA_Ti, yerr=sigma_EA_Ti, fmt='o', ecolor='red', capsize=5, label='Ti', color='red')
    plt.errorbar(T, EA_Zr, yerr=sigma_EA_Zr, fmt='s', ecolor='blue', capsize=5, label='Zr', color='blue')
    
    plt.legend()
    plt.savefig(OUTPUT_IMAGE, dpi=300)
    plt.show()
    plt.close()

if __name__ == "__main__":
    main()
