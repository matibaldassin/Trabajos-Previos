#/bin/bash


for ii in `seq 10 20 590`


do


cp ../cooling-2-$ii/P-local entrada.in

./loc

mv salida.out P-local-2-$ii 

rm entrada.in

done



