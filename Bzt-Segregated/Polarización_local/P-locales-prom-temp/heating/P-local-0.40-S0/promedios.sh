#!/bin/bash


for ii in `seq 590 20 790`


do


cp ../heating-2-$ii/P-local entrada.in

./loc

mv salida.out P-local-$ii 

rm entrada.in

done



