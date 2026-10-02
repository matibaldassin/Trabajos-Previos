#!/bin/bash
CONC="0.40"
SEG="S100"

cp ../CONFIG_files/CONFIG-$CONC-$SEG CONFIG

for i in $(seq 10 10 400)
do
  cp ../p-locales/p-locales-$CONC-$SEG/P-local-$i P-local
  gfortran ovito3.f90 -o ov3.x

  ./ov3.x
  mv ovito.xyz ovito-$i.xyz

done

for i in $(seq 410 20 590)
do
  cp ../p-locales/p-locales-$CONC-$SEG/P-local-$i P-local
  gfortran ovito3.f90 -o ov3.x

  ./ov3.x
  mv ovito.xyz ovito-$i.xyz

done

mv ovito-*.xyz $CONC-$SEG
rm -f ov3.x P-local CONFIG
