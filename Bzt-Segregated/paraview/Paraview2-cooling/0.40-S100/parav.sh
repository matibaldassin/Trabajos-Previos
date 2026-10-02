#!/bin/bash


for ii in `seq 10 20 400`


do


cp ../promedio-polarizaciones-locales/P-local-$ii entrada.dat

./parav

mv salida-no-angles-deltaP.vtk Pol-$ii.vtk 

done



