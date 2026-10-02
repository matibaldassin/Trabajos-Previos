import numpy as np
import re

File1 = 'Resultados-0.40-S0/cell_lengths_summary_pos-prom-20-0.txt'
File2 = 'Resultados-0.40-S50/cell_lengths_summary_pos-prom-20-0.txt'
File3 = 'Resultados-0.40-S100/cell_lengths_summary_pos-prom-20-0.txt'

def parse_summary(filepath):
    """Extrae medias y stds de ax, ay, az para Ti y Zr del archivo summary."""
    with open(filepath) as f:
        text = f.read()

    def extract(prefix):
        vals = []
        for axis in ('a_x', 'a_y', 'a_z'):
            tag = prefix.replace(r'\[', '[').replace(r'\]', ']') + axis
            for line in text.splitlines():
                if tag in line and 'mean' in line and 'std' in line:
                    # extrae todos los números de la parte después de ':'
                    after_colon = line.split(':', 1)[1]
                    nums = re.findall(r'[\d]+\.[\d]+', after_colon)
                    if len(nums) < 2:
                        raise ValueError(f"No se pudieron leer valores en línea: {line!r}")
                    vals.append((float(nums[0]), float(nums[1])))
                    break
            else:
                raise ValueError(f"No se encontró '{tag}' en {filepath}")
        return vals  # [(ax_mean, ax_std), (ay_mean, ay_std), (az_mean, az_std)]

    ti = extract(r'\[Ti\] ')
    zr = extract(r'\[Zr\] ')
    return ti, zr

def volumen(params):
    """Calcula V = ax*ay*az y σ_V por propagación de errores."""
    (ax, sax), (ay, say), (az, saz) = params
    V = ax * ay * az
    sV = V * np.sqrt((sax/ax)**2 + (say/ay)**2 + (saz/az)**2)
    return V, sV

samples = [('S0', File1), ('S50', File2), ('S100', File3)]

header = f"{'Sample':<8}  {'V_Ti (Å³)':<18} {'σ_Ti (Å³)':<16}  {'V_Zr (Å³)':<18} {'σ_Zr (Å³)'}"
print(header)
print("-" * len(header))

for label, filepath in samples:
    ti_params, zr_params = parse_summary(filepath)
    V_Ti, sV_Ti = volumen(ti_params)
    V_Zr, sV_Zr = volumen(zr_params)
    print(f"{label:<8}  {V_Ti:<18.4f} {sV_Ti:<16.4f}  {V_Zr:<18.4f} {sV_Zr:.4f}")


