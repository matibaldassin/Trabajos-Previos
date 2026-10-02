"""
01_parse_and_compute.py
========================
Lee el archivo VTK de polarizaciones, extrae coordenadas y vectores,
y calcula rotor (curl) y divergencia sobre el campo de polarizacion
COMPLETO (Ti+Zr), sin enmascarar ni restringir a un sitio en particular.

Las derivadas espaciales (curl, div) se calculan con diferencias
centradas periódicas (np.roll), asumiendo condiciones periódicas de
contorno en las tres direcciones (celda simulada periódica).

Outputs (archivos .npy para los demás scripts):
  P_ti.npy       -- polarización enmascarada, Zr=0  (Nx, Ny, Nz, 3)
  ti_mask.npy    -- máscara booleana de sitios Ti   (Nx, Ny, Nz)
  id_grid.npy    -- campo identidad original         (Nx, Ny, Nz)
  curl_x/y/z.npy -- componentes del rotor
  curl_mag.npy   -- magnitud del rotor |∇×P|
  div.npy        -- divergencia ∇·P
  dims.npy       -- [Nx, Ny, Nz]
"""

import numpy as np

# ──────────────────────────────────────────────
# PARÁMETRO: ruta al archivo VTK
# ──────────────────────────────────────────────
filepath = "Pol-10.vtk"   # <-- cambiá esto al path de tu archivo

# ──────────────────────────────────────────────
# 1. PARSEO DEL ARCHIVO VTK
# ──────────────────────────────────────────────
with open(filepath, 'r') as f:
    lines = f.readlines()

dims = None
points_start = scalars_id_start = vectors_start = None

for i, line in enumerate(lines):
    if line.strip().startswith('DIMENSIONS'):
        dims = list(map(int, line.strip().split()[1:]))
    if line.strip().startswith('Points'):
        points_start = i + 1
    if 'SCALARS identidad' in line:
        scalars_id_start = i + 2   # saltar línea LOOKUP_TABLE
    if 'VECTORS polar_vector' in line:
        vectors_start = i + 1

Nx, Ny, Nz = dims
total = Nx * Ny * Nz
print(f"Dimensiones: {dims}  ({total} puntos)")

# Parsear vectores de polarización
vecs = []
i = vectors_start
while len(vecs) < total:
    vals = list(map(float, lines[i].split()))
    for j in range(0, len(vals), 3):
        if j + 2 < len(vals):
            vecs.append(vals[j:j+3])
    i += 1
vecs = np.array(vecs)

# Parsear campo identidad (1=Zr, 2=Ti)
identidad = []
i = scalars_id_start
while len(identidad) < total:
    vals = list(map(float, lines[i].split()))
    identidad.extend(vals)
    i += 1
identidad = np.array(identidad[:total])

# ──────────────────────────────────────────────
# 2. COORDENADAS → ÍNDICES ESPACIALES
# ──────────────────────────────────────────────
# Parsear las coordenadas (x, y, z) de cada punto para mapear
# correctamente cada dato a su posición, sin asumir orden de almacenamiento.
coords_list = []
i = points_start
while len(coords_list) < total:
    vals = list(map(float, lines[i].split()))
    for j in range(0, len(vals), 3):
        if j + 2 < len(vals):
            coords_list.append(vals[j:j+3])
    i += 1
coords = np.round(np.array(coords_list), 3)   # (total, 3): columnas = x, y, z

x_vals = np.sort(np.unique(coords[:, 0]))
y_vals = np.sort(np.unique(coords[:, 1]))
z_vals = np.sort(np.unique(coords[:, 2]))

assert len(x_vals) == Nx and len(y_vals) == Ny and len(z_vals) == Nz, \
    f"Dimensiones mismatch: {len(x_vals)}×{len(y_vals)}×{len(z_vals)} vs {Nx}×{Ny}×{Nz}"

print(f"x ∈ [{x_vals[0]:.1f}, {x_vals[-1]:.1f}]  "
      f"y ∈ [{y_vals[0]:.1f}, {y_vals[-1]:.1f}]  "
      f"z ∈ [{z_vals[0]:.1f}, {z_vals[-1]:.1f}]")

xi = np.searchsorted(x_vals, coords[:, 0])
yi = np.searchsorted(y_vals, coords[:, 1])
zi = np.searchsorted(z_vals, coords[:, 2])

P_full  = np.zeros((Nx, Ny, Nz, 3))
id_grid = np.zeros((Nx, Ny, Nz))
P_full[xi, yi, zi]  = vecs
id_grid[xi, yi, zi] = identidad

ti_mask = (id_grid == 2)   # True en sitios Ti

# Enmascarar: Zr → cero (no contribuye a gradientes)
P_ti = P_full * ti_mask[:, :, :, np.newaxis]

dx = 4e-8  # spacing de la grilla en cm (4.0 A = 4e-8 cm); P esta en uC/cm2,
           # asi que dx tiene que estar en cm para que rot/div P den uC/cm3
print(f"Sitios Ti: {ti_mask.sum()}  |  Sitios Zr: {(~ti_mask).sum()}")

# ──────────────────────────────────────────────
# 3. CAMPOS DERIVADOS (grilla completa, Ti+Zr)
# ──────────────────────────────────────────────
Px = P_full[:, :, :, 0]
Py = P_full[:, :, :, 1]
Pz = P_full[:, :, :, 2]


def _grad_periodic(f, dx, axis):
    """Derivada central con condiciones periódicas de contorno (np.roll)."""
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2 * dx)


# Rotor: ∇×P
curl_x = _grad_periodic(Pz, dx, axis=1) - _grad_periodic(Py, dx, axis=2)
curl_y = _grad_periodic(Px, dx, axis=2) - _grad_periodic(Pz, dx, axis=0)
curl_z = _grad_periodic(Py, dx, axis=0) - _grad_periodic(Px, dx, axis=1)
curl_mag = np.sqrt(curl_x**2 + curl_y**2 + curl_z**2)

# Divergencia: ∇·P
div = (_grad_periodic(Px, dx, axis=0) +
       _grad_periodic(Py, dx, axis=1) +
       _grad_periodic(Pz, dx, axis=2))

# Estadísticas sobre el campo completo (Ti+Zr)
pol_all  = np.linalg.norm(P_full, axis=-1).ravel()
curl_all = curl_mag.ravel()

print("\n=== ESTADÍSTICAS (campo completo, Ti+Zr) ===")
print(f"|P|   medio : {pol_all.mean():.4f}")
print(f"|∇×P| medio : {curl_all.mean():.4f}")
print(f"|∇×P|/|P|   : {curl_all.mean()/pol_all.mean()*100:.1f}%")
print(f"∇·P   medio : {div.mean():.4f}")


# ──────────────────────────────────────────────
# 4. GUARDAR PARA LOS DEMÁS SCRIPTS
# ──────────────────────────────────────────────
np.save('P_ti.npy',      P_ti)
np.save('P_full.npy',    P_full)
np.save('ti_mask.npy',   ti_mask)
np.save('id_grid.npy',   id_grid)
np.save('curl_x.npy',    curl_x)
np.save('curl_y.npy',    curl_y)
np.save('curl_z.npy',    curl_z)
np.save('curl_mag.npy',  curl_mag)
np.save('div.npy',       div)
np.save('dims.npy',      np.array(dims))

print("\nArchivos guardados correctamente.")
