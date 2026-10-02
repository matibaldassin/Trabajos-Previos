#!/bin/bash

# Directorio origen
ORIGEN="p-locales/heating"

# Crear un arreglo con los archivos a copiar según conc y tipo
# Usar globbing correctamente: activar nullglob para que patrones sin coincidencias
# se expandan a lista vacía en lugar de quedarse sin expandir.
shopt -s nullglob
for conc in 0.20 0.40; do
    for tipo in sol esf esfsol; do
        dir="$ORIGEN/p-locales-${conc}-${tipo}"
        # iterar sobre todos los archivos que coincidan con el patrón
        for archivo in "$dir"/P-local-*; do
            # si no hay coincidencias, el bucle no se ejecuta (nullglob)
            if [ -e "$archivo" ]; then
                cp "$archivo" .
                python3 polarizaciones-rms.py
            else
                # esto normalmente no ocurrirá gracias a nullglob, pero queda como fallback
                echo "⚠ Archivo no encontrado: $archivo"
            fi
        done
    done
done
shopt -u nullglob

# Borrar los archivos copiados
rm -f P-local-*

