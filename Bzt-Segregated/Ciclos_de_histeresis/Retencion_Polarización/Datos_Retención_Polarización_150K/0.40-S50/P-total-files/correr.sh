#!/bin/bash
cp ../../menos-0_kV_por_mm/P-total P-total-0

for i in $(seq 1 9); do
    cp ../block-000${i}/P-total P-total-$i
done

for i in $(seq 10 15); do  
    cp ../block-00${i}/P-total P-total-$i
done
