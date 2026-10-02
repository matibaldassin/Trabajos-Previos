# ------ SCRIPT PARAMS ---------
CONC="0.40"
SEG="S0"

# Input files (histogram tables)
folder = f"/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/MultiSOM/Phases_Data/{CONC}-{SEG}"
files = [f"{folder}/Data-{i}.txt" for i in range(50)]

# Temperature values
T_vals = list(range(10,410,10))+list(range(410,600,20))

# Data output
data_file = f'Conteo_de_Fases/Fases_locales_{CONC}_{SEG}.txt'

# ------------------------------

# Script:
# Import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Group map
group_dict = {0: 'Tetra', 1: 'Rhombo',
              2: 'Ortho', 3: 'Cubic', 4: 'Zr'}

# Init data
df_dict = {'Temp':T_vals} 

for phase in group_dict.values():
    df_dict[phase] = []  # Crear una lista vacía para cada fase

# Read data
for f in files:
    data = pd.read_csv(f, comment='#', sep=' ',
                       names=['Group', 'Count'],
                       index_col=False)
    print(f"Data from file {f}:\n{data}")
   
 # Inicializar los conteos para este archivo
    counts = {phase: 0 for phase in group_dict.values()}
    
    # Para cada fase, obtener el conteo
    for i in range(4):
        phase = group_dict[i]
        count = data[data['Group'].sub(i).abs() < 0.5]['Count'].sum()
        counts[phase] = count
    
    # Agregar los valores al diccionario
    for phase, count in counts.items():
        df_dict[phase].append(count)
    print(f"Counts for file {f}: {counts}")

# Convertir a DataFrame
df = pd.DataFrame(df_dict)


# Write data
df.to_csv(data_file, index=False)
