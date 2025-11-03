#!/bin/bash

function raspa_input_file {
 echo "#number of pseudo atoms
5
#type      print   as    chem  oxidation   mass        charge   polarization B-factor radii  connectivity anisotropic anisotropic-type   tinker-type
O1         yes     O1     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O2         yes     O2     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O3         yes     O3     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
Si         yes     Si    Si    0           28.0855     0.0      0.0          1.0      1.18   4            0           absolute           0
Ge         yes     Ge    Ge    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
" > pseudo_atoms.def 
echo "SimulationType                MC
NumberOfCycles                0
NumberOfInitializationCycles  0
PrintEvery                    10
     
Forcefield                    Local

RemoveAtomNumberCodeFromLabel yes

ModifyFrameworkAtomConnectedTo O O2 Si
ModifyFrameworkAtomConnectedTo O O1 Ge
ModifyFrameworkAtomConnectedTo O2 O1 Ge
ModifyFrameworkAtomConnectedTo O3 O1 Si
     
Framework 0
FrameworkName Local
UnitCells 1 1 1
ExternalTemperature 300.0
" > simulation.input
}

file=$1

# VASP 2 CIF
if [ -f Local.cif ] ; then rm -rf Local.cif ; fi
raspa_input_file
ase convert -i vasp -o cif $file Local.cif
sed -i -e 's/_space_group_IT_number/_symmetry_Int_Tables_number/g' -e '/space_group_name_H-M_alt/d' Local.cif 
simulate 
cp Movies/System_0/Framework_0_final_1_1_1_P1.cif ${file}.cif
rm -rf Local.cif VTK Restart Output Movies pseudo_atoms.def pseudo_atoms.def simulation.input

# CIF 2 GIN
echo "opti conp prop
name $file
cell" > $file.gin
cell1=$(grep '_cell_length_a' $file.cif | awk '{print $2}'); cell2=$(grep '_cell_length_b' $file.cif | awk '{print $2}'); cell3=$(grep '_cell_length_c' $file.cif | awk '{print $2}');
cell4=$(grep '_cell_angle_alpha' $file.cif | awk '{print $2}'); cell5=$(grep '_cell_angle_beta' $file.cif | awk '{print $2}'); cell6=$(grep '_cell_angle_gamma' $file.cif | awk '{print $2}')
echo " $cell1 $cell2 $cell3 $cell4 $cell5 $cell6
frac " >> $file.gin
sed -n '/^_atom_site_charge$/,$P' $file.cif | sed '/_atom_site_charge'/d | sed '/^[[:space:]]*$/d' | awk '{print $1 " core",$3,$4,$5}' >> $file.gin
echo "
space
P 1

species
Ge core  Ge core
Si core  Si core
O1 core  O1 core
O2 core  O2 core
O3 core  O3 core
O1 shel  O1 shel
O2 shel  O2 shel
O3 shel  O3 shel
end
library Germanate
switch_minimiser rfo gnorm 0.05
stepmx opt 0.1
dump every 1 $file.grs
output cif $file.cif" >> $file.gin

exit 0
