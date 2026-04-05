#!/bin/bash
gfortran -O2 energy_stats.f90 -o energy_stats.exe
if [ -f energy_stats.txt ] ; then rm -rf energy_stats.txt ; touch energy_stats.txt ; fi
#loc=$(pwd)
#for i in ../hpm-18_mc_uniform_*/mc_samples_N* ; do
# cd $i
#  bash extract.sh ; mv ENERGIES ENERGIES.bak; awk '{if ($1<0) print $1}' ENERGIES.bak | sed '/#/d' > c; mv c ENERGIES.bak
# cd $loc
#done
for i in n0*/ENERGIES mc_samples_N*/ENERGIES mc_samples_N*/ENERGIES.bak ; do
# Archivo: ENERGIES
#  Entradas: 1000
#  <E> (media):     -3005.52538643
#  Desviacion estandar:         0.25276191
 ./energy_stats.exe -T 300.0 $i > tmp
# <E>_Boltzmann:     -3006.04870191
# Sigma_Boltzmann:         0.15280889
 name=$(grep 'Archivo' tmp | awk '{print $2}')
 E=$(grep '<E>_Boltzmann' tmp | awk '{print $2}')
 dE=$(grep 'Sigma_Boltzmann' tmp | awk '{print $2}')
 rm -rf tmp
 echo $E $dE '#' $name >> energy_stats.txt 
done
sed -i -e 's/\/ENERGIES//g' -e 's/# n//g' -e 's/# mc_samples_N00//g' -e 's/\.bak//g' energy_stats.txt 
sort -gk3 energy_stats.txt > c; mv c energy_stats.txt
cat energy_stats.txt
#echo "set autoscale; nGeO2 = -121.812200998519; nSiO2 = -128.799143746482; N=24.0;p '< sort -gk3 energy_stats.txt' u (\$3/N):(96.48533212331*(\$1-\$3*nGeO2-(N-\$3)*nSiO2)/N):(96.48533212331*sqrt(\$2)/N) w yerrorbars pt 7 ps 1 lc rgb 'red' notitle" > gp.tmp
