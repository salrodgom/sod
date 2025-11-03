#!/bin/bash
export OMP_NUM_THREADS=8,1
export OMP_STACKSIZE=8G
ulimit -s unlimited

for file in *.vasp ; do
 rm -rf *.gen
 ase convert -i vasp -o gen $file geo.gen
 dftb+ | tee $file.log
 ase convert -i gen -o cif geo_end.gen $file.cif
done
exit 0
