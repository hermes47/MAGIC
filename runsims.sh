#!/bin/bash

MOLECULES=`ls OutputFiles/ | sed s/'.*Files\/'/''/`

for mol in $MOLECULES;
do
cd ~/GitHub/MAGIC/OutputFiles/$mol
make_top @build $mol.mtb ../../GROMOSData/54A7.mtb @param ../../GROMOSData/54A7.ifp @seq $mol @solv H2O > $mol.top
min=$mol'_minimised.pdb'
pdb2g96 @topo $mol.top @pdb $min @lib 
done
