import numpy as np
import matplotlib.pyplot as plt
import csv

IN_FOLDER = "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.40-S100"

OUTPUT_FILE = "Polarización-local-promedio-0.40-S100.dat"
OUTPUT_IMAGE = "Polarización-local-promedio-0.40-S100.png"

def Edward_Anderson(f):
    lector = csv.reader(f, delimiter=' ')
    Pol = []

    for fila in lector:
        aux = [float(x) for x in fila if x]
        Pol.append(aux[4])

    Pol = np.array(Pol)

    Pol_mean_value = np.mean(Pol)
    Pol_std = np.std(Pol)

    return Pol_mean_value, Pol_std


def main():
    temperatures = list(range(10,401,10)) + list(range(410,591,20))
    T, P_loc, sigma_P_loc = [],[],[]

    for i in temperatures:
        with open(f"{IN_FOLDER}/P-local-{i}", "r", encoding="utf-8") as s:
            p_loc, s_loc = Edward_Anderson(s)

            T.append(i)
            P_loc.append(p_loc)
            sigma_P_loc.append(s_loc)

    with open(OUTPUT_FILE,"w") as r:
        r.write(f"{'# T [K]':>10}\t{'P_local':>15}\t{'Sigma':>15}\n")
        for i in range(len(T)):
            r.write(f"{T[i]:>10.2f}\t{P_loc[i]:>15.8f}\t{sigma_P_loc[i]:>15.8f}\n")

    T=1.48*np.array(T)
    plt.figure(figsize=(10,6))
    plt.title('Meanvalue of Local Polarization vs Temperature')
    plt.xlabel(r"$Temperature \, [K]$")
    plt.ylabel(r"$P_{local} \, [\mu C/cm^2]$")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.minorticks_on()

    plt.plot(T, P_loc, marker='o', label='P_local', color='red')
    
    plt.legend()
    plt.savefig(OUTPUT_IMAGE, dpi=300)
    plt.show()
    plt.close()

if __name__ == "__main__":
    main()
