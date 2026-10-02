#!/bin/bash

python3 tamaño-celdas.py --config configs/CONFIG-0.40-S0 --tabpol tabpols/tabpol-0.40-S0.dat --lattice Averege_cell_vectors_20K/cell-vectors-20-0.40-S0 --pos Promedio_de_Posiciones_20K/pos-prom-20-0.40-S0 --outdir Resultados-0.40-S0

python3 tamaño-celdas.py --config configs/CONFIG-0.40-S50 --tabpol tabpols/tabpol-0.40-S50.dat --lattice Averege_cell_vectors_20K/cell-vectors-20-0.40-S50 --pos Promedio_de_Posiciones_20K/pos-prom-20-0.40-S50 --outdir Resultados-0.40-S50

python3 tamaño-celdas.py --config configs/CONFIG-0.40-S100 --tabpol tabpols/tabpol-0.40-S100.dat --lattice Averege_cell_vectors_20K/cell-vectors-20-0.40-S100 --pos Promedio_de_Posiciones_20K/pos-prom-20-0.40-S100 --outdir Resultados-0.40-S100


