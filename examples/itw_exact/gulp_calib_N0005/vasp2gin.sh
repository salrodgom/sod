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
ModifyFrameworkAtomConnectedTo O1 O3 Si
     
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

function run_ase_convert {
	local src="$1"
	local dst="$2"
	local ase_bin

	if [ -x /home/salvador/miniforge3/bin/ase ] ; then
		/home/salvador/miniforge3/bin/ase convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
	fi

	if [ -x /home/salvador/.local/bin/ase ] ; then
		/home/salvador/.local/bin/ase convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
	fi

	if command -v ase >/dev/null 2>&1 ; then
		ase_bin=$(command -v ase)
		if [ -n "$ase_bin" ] ; then
			"$ase_bin" convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	fi

	for pycmd in python3 python ; do
		if command -v "$pycmd" >/dev/null 2>&1 ; then
			"$pycmd" -c "import ase" >/dev/null 2>&1 || continue
			"$pycmd" -m ase.cli.main convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	done

	echo "Error: no se pudo ejecutar ASE para convertir $src a CIF" >&2
	return 1
}

run_ase_convert "$file" Local.cif || exit 1
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
echo "# O1: Ge-O1-Ge
# O2: Si-O2-Si
# O3: Si-O3-Ge
species
Si    core  4.00000
O2    core  0.870733
O2    shel -2.870733
O3    core  1.330431
O3    shel -3.330431
O1    core  1.733957
O1    shel -3.733957
F     core  0.56
F     shel -1.56
Ge    core  4.0
buckingham
Si core    O1 shel     1315.2478 0.317759 10.141118 0.0 16.0
Si core    O2 shel     1315.2478 0.317759 10.141118 0.0 16.0
Si core    O3 shel     1315.2478 0.317759 10.141118 0.0 16.0
Ge core    O1 shel     1497.3996 0.325646 16.808599 0.0 16.0
Ge core    O2 shel     1497.3996 0.325646 16.808599 0.0 16.0
Ge core    O3 shel     1497.3996 0.325646 16.808599 0.0 16.0
O1 shel    O1 shel    22764.0000 0.149000 10.937044 0.0 16.0
O1 shel    O2 shel    22764.0000 0.149000 10.937044 0.0 16.0
O1 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0
O2 shel    O2 shel    22764.0000 0.149000 10.937044 0.0 16.0
O2 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0
O3 shel    O3 shel    22764.0000 0.149000 10.937044 0.0 16.0
Si core     F shel     976.82887 0.282000 0.0       0.0 16.0
Ge core     F shel     681.47288 0.320000 0.0       0.0 16.0
F  shel     F shel     540.39761 0.262490 0.0       0.0 16.0
O1 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0
O2 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0
O3 shel     F shel    1675.00000 0.268000 0.0       0.0 16.0
spring
O1  180.315770
O2   75.96980
O3  128.1427
F    33.452757
three
Si core O1 shel O1 shel 1.2614 109.47 1.9 1.9 3.5
Si core O1 shel O2 shel 1.2614 109.47 1.9 1.9 3.5
Si core O2 shel O2 shel 1.2614 109.47 1.9 1.9 3.5
" > Germanate.lib
exit 0
