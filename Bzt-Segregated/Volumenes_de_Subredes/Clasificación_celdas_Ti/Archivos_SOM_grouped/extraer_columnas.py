import os

# Columnas a extraer (se pueden cambiar según el encabezado del archivo)
COLUMNAS = ["DipoleMagnitude", "Group", "id"]

# Directorio con los archivos .dump
DIRECTORIO = "."

def extraer_columnas(ruta_entrada, columnas):
    with open(ruta_entrada, "r") as f:
        lineas = f.readlines()

    # Encontrar la línea con los nombres de columna
    indices = {}
    datos = []
    en_atomos = False

    for linea in lineas:
        linea = linea.strip()
        if linea.startswith("ITEM: ATOMS"):
            encabezados = linea.split()[2:]  # quitar "ITEM:" y "ATOMS"
            for col in columnas:
                if col not in encabezados:
                    raise ValueError(f"Columna '{col}' no encontrada. Columnas disponibles: {encabezados}")
                indices[col] = encabezados.index(col)
            en_atomos = True
            continue
        if linea.startswith("ITEM:"):
            en_atomos = False
            continue
        if en_atomos and linea:
            valores = linea.split()
            fila = {col: valores[indices[col]] for col in columnas}
            datos.append(fila)

    return datos


def guardar_csv(datos, columnas, ruta_salida):
    with open(ruta_salida, "w") as f:
        f.write(" ".join(columnas) + "\n")
        for fila in datos:
            f.write(" ".join(fila[col] for col in columnas) + "\n")
    print(f"Guardado: {ruta_salida}  ({len(datos)} átomos)")


archivos = [f for f in os.listdir(DIRECTORIO) if f.endswith(".dump")]

for archivo in sorted(archivos):
    ruta_entrada = os.path.join(DIRECTORIO, archivo)
    nombre_salida = archivo.replace(".dump", "_extraido.txt")
    ruta_salida = os.path.join(DIRECTORIO, nombre_salida)

    datos = extraer_columnas(ruta_entrada, COLUMNAS)
    guardar_csv(datos, COLUMNAS, ruta_salida)
