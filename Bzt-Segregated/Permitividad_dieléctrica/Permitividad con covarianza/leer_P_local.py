"""
Lee los archivos P-local-*.txt (carpeta p-local/) y separa cada paso temporal
de las N_SITIOS polarizaciones locales.

Formato de cada archivo: bloques consecutivos, cada uno con
  - 1 linea de encabezado: numero de paso temporal (ej. 10000, 11000, ...)
  - N_SITIOS lineas de datos: indice_sitio Px Py Pz

Para cada archivo se obtiene:
  pasos : array (n_pasos,)               numero de paso de cada bloque
  P     : array (n_pasos, N_SITIOS, 3)    polarizaciones locales (Px, Py, Pz)
"""

import re
from pathlib import Path

import numpy as np

N_SITIOS = 21952
CARPETA = Path(__file__).parent / "p-local"


def leer_P_local(filepath, n_sitios=N_SITIOS):
    """Separa un archivo P-local-XX.txt en sus pasos temporales.

    Devuelve (pasos, P) con pasos.shape == (n_pasos,)
    y P.shape == (n_pasos, n_sitios, 3).
    """
    with open(filepath) as f:
        lineas = f.read().splitlines()

    stride = n_sitios + 1  # 1 encabezado (paso) + n_sitios lineas de datos
    if len(lineas) % stride != 0:
        raise ValueError(
            f"{filepath}: {len(lineas)} lineas no es multiplo de {stride} "
            f"(1 encabezado + {n_sitios} sitios)"
        )
    n_pasos = len(lineas) // stride

    idx_encabezados = set(range(0, len(lineas), stride))
    pasos = np.array([int(lineas[i]) for i in range(0, len(lineas), stride)], dtype=np.int64)

    lineas_datos = [l for i, l in enumerate(lineas) if i not in idx_encabezados]
    bloque = np.loadtxt(lineas_datos, dtype=np.float64)  # columnas: indice, Px, Py, Pz
    P = bloque[:, 1:4].reshape(n_pasos, n_sitios, 3)

    return pasos, P


# ──────────────────────────────────────────────
# Lectura de todos los archivos P-local-*.txt
# ──────────────────────────────────────────────
datos = {}
for filepath in sorted(CARPETA.glob("P-local-*.txt")):
    m = re.search(r"P-local-(\d+)\.txt", filepath.name)
    clave = int(m.group(1))
    pasos, P = leer_P_local(filepath)
    datos[clave] = {"pasos": pasos, "P": P}
    print(f"{filepath.name}: {len(pasos)} pasos temporales, P.shape={P.shape}")

# A partir de aca `datos` tiene, por cada corrida (ej. datos[70]):
#   datos[70]["pasos"] -> array (n_pasos,) con el paso temporal de cada bloque
#   datos[70]["P"]     -> array (n_pasos, 21952, 3) con las polarizaciones locales
