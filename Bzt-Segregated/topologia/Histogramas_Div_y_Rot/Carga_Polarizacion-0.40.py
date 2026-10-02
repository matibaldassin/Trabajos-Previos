import numpy as np

# INPUT FILES
COMPOSICION = '0.40'
SISTEMAS = ['S0', 'S50', 'S100']
FILE_DIV = {s: f'../{COMPOSICION}-{s}/Cooling/div.npy' for s in SISTEMAS}
FILE_ID  = {s: f'../{COMPOSICION}-{s}/Cooling/id_grid.npy' for s in SISTEMAS}

# V_celda = 64 A^3 = 64*(1e-8 cm)^3, el volumen promedio de celda unidad usado para dx=4e-8.
V_CELL_A3 = 64.0
V_CELL = V_CELL_A3 * (1e-8) ** 3  # cm^3
N_CELLS = 21952  # 28x28x28, total de celdas del sistema

# ── Carga y separación por especie ────────────────
datos_div, medias_div = {}, {}
electrones_por_celda = {}
electrones_pos, electrones_neg = {}, {}
div_raw_todos, div_electrones_todos, id_todos = {}, {}, {}
for sistema in SISTEMAS:
    id_grid = np.load(FILE_ID[sistema])
    zr, ti = id_grid == 1, id_grid == 2

    div = np.load(FILE_DIV[sistema])
    datos_div[sistema] = (div[zr], div[ti])
    medias_div[sistema] = (div[zr].mean(), div[ti].mean())

    div_electrones = div * V_CELL / 1.6e-13  # uC/cm^3 -> electrones/celda
    electrones_por_celda[sistema] = (div_electrones).sum() / N_CELLS  # promedio de electrones de polarización por celda
    electrones_pos[sistema] = div_electrones[div_electrones > 0].sum()
    electrones_neg[sistema] = div_electrones[div_electrones < 0].sum()

    div_raw_todos[sistema] = div.ravel()
    div_electrones_todos[sistema] = div_electrones.ravel()
    id_todos[sistema] = id_grid.ravel()

# ── Volcado de todos los valores por celda (crudos y en electrones/celda) ──
OUT_VALORES = f'Div_valores_por_celda-{COMPOSICION}.txt'
with open(OUT_VALORES, 'w') as f:
    f.write(f'{"sistema":<8}{"id":>4}{"div_uC_cm3":>18}{"div_electrones":>18}\n')
    for sistema in SISTEMAS:
        for i, d, de in zip(id_todos[sistema], div_raw_todos[sistema], div_electrones_todos[sistema]):
            f.write(f'{sistema:<8}{int(i):>4d}{d:>18.6e}{de:>18.6e}\n')
print(f"Guardado: {OUT_VALORES}")

print(f'\n=== Media de la divergencia ∇·P — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Zr":>14}{"Ti":>14}')
for sistema in SISTEMAS:
    m_zr, m_ti = medias_div[sistema]
    print(f'{sistema:<8}{m_zr:>14.4e}{m_ti:>14.4e}')

# Carga de polarización: div.npy ya esta en uC/cm^3 (rho = div P), asi que
# integrar sobre el volumen es sumar cada celda * V_celda.
print(f'\n=== Carga de polarización por celda, V_celda={V_CELL_A3:g} Å³ — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Q_Zr":>16}{"Q_Ti":>16}{"Q_total (electrones por celda)":>16}')
for sistema in SISTEMAS:
    div_zr, div_ti = datos_div[sistema]
    q_zr = (div_zr / 1.6e-13).sum() * V_CELL
    q_ti = (div_ti / 1.6e-13).sum() * V_CELL
    print(f'{sistema:<8}{q_zr:>16.4e}{q_ti:>16.4e}{q_zr + q_ti:>24.17g}')

print(f'\n=== Electrones de polarización por celda (Σ(∇·P en e⁻/celda) / {N_CELLS}) — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"e- por celda":>16}')
for sistema in SISTEMAS:
    print(f'{sistema:<8}{electrones_por_celda[sistema]:>16.4e}')

print(f'\n=== Suma de electrones de polarización, positivos y negativos por separado — x={COMPOSICION} ===')
print(f'{"Sistema":<8}{"Σ positivos":>16}{"Σ negativos":>16}{"suma total":>16}')
for sistema in SISTEMAS:
    pos, neg = electrones_pos[sistema], electrones_neg[sistema]
    print(f'{sistema:<8}{pos:>16.4e}{neg:>16.4e}{pos + neg:>24.17g}')
