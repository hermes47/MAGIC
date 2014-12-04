#!/bin/bash

LIPIDS=`cat lipidstomake`

for lipid in $LIPIDS;
do
./magic.py -s "GLYM(PCHM$lipid)" >> mylogfile.txt 2>&1
./magic.py -s "GLYM(PET2$lipid)" >> mylogfile.txt 2>&1
./magic.py -s "GLYM(PSEM$lipid)" >> mylogfile.txt 2>&1
./magic.py -s "GLYM(PETM$lipid)" >> mylogfile.txt 2>&1
done

PRODUCEDFILES=`ls OutputFiles/ | sed s/'.*OutputFiles\/'/''/`

for file in $PRODUCEDFILES;
do
minimised=$file'_minimised.pdb'
if [ -e "OutputFiles/$file/$minimised" ]
then
echo $file
head -n 1 OutputFiles/$file/$file.mtb
fi
done
