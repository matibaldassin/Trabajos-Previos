#!/bin/bash


for ii in `seq 10 20 590`


do


python3 reorder-vtk-points.py ../Pol-${ii}.vtk -o Pol-${ii}.vtk


done



