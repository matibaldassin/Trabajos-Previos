import os
import numpy as np
import csv

# Numero de titanios y circonios
NumZr=[4390,4390,4364,8780,8780,8902]
NumTi=[17562,17562,17588,13172,13172,13050]

# INPUT FOLDERS
IN_FOLDERS=[
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.20-S0/",
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.20-S50/",
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.20-S100/",
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.40-S0/",
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.40-S50/",
    "/home/matias/Escritorio/CONICET/bzt/Polarización_local/P-locales-prom-temp/heating/P-local-0.40-S100/",
]

# OUTPUT FILES
OUT_FILES=[
    "P-local-por-sitio-B-0.20-S0.dat",
    "P-local-por-sitio-B-0.20-S50.dat",
    "P-local-por-sitio-B-0.20-S100.dat",
    "P-local-por-sitio-B-0.40-S0.dat",
    "P-local-por-sitio-B-0.40-S50.dat",
    "P-local-por-sitio-B-0.40-S100.dat",
]

# READING FUNCTIONS
def variables(f,i):
    lector = csv.reader(f,delimiter=' ')
    px_zr, py_zr, pz_zr, p_zr = [], [], [], []
    px_ti, py_ti, pz_ti, p_ti = [], [], [], []
    
    filas = []
    for fila in lector:
        aux = [float(x) for x in fila if x]
        if aux:
            filas.append(aux)
    
    if len(filas) < NumZr[i] + NumTi[i]:
        print("Erros: El archivo tiene menos filas de las esperadas.")
    
    for aux in filas[:NumZr[i]]:
        px_zr.append(aux[1])
        py_zr.append(aux[2])
        pz_zr.append(aux[3])
        p_zr.append(aux[4])
    
    for aux in filas[NumZr[i]:NumZr[i]+NumTi[i]]:
        px_ti.append(aux[1])
        py_ti.append(aux[2])
        pz_ti.append(aux[3])
        p_ti.append(aux[4])
    
    p_zr_prom = np.mean(p_zr)
    P_zr_std = np.std(p_zr)
    p_ti_prom = np.mean(p_ti)
    P_ti_std = np.std(p_ti)

    return p_zr_prom, P_zr_std, p_ti_prom, P_ti_std

Temp=list(range(10,401,10)) + list(range(410,591,20))

for i in range(len(IN_FOLDERS)):
    header = f"{'T [K]':>10}\t{'P_local_Ti':>15}\t{'P_local_Zr':>15}\t{'Sigma_Ti':>15}\t{'Sigma_Zr':>15}\n"
    with open(OUT_FILES[i], "w", encoding="utf-8") as r:
        r.write(header)
        for j in Temp:
            input_path = os.path.join(IN_FOLDERS[i], f"P-local-{j}")
            with open(input_path, "r", encoding="utf-8") as s:
                p_zr_prom, p_zr_std, p_ti_prom, p_ti_std = variables(s, i)
                r.write(f"{j:>10.2f}\t{p_ti_prom:>15.8f}\t{p_zr_prom:>15.8f}\t{p_ti_std:>15.8f}\t{p_zr_std:>15.8f}\n")



