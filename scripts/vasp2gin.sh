#!/bin/bash
#*******************************************************************************
#    Copyright (c) 2025, Salvador R.G. Balestra
#
#    This file is part of the SOD package.
#
#    SOD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#******************************************************************************

# Converts a POSCAR-like VASP file into the GULP input bundle required by SOD/GULP calibrations.
set -euo pipefail

# Validate inputs and capture the path to the structure that must be converted.
if [ $# -ne 1 ]; then
	echo "Uso: $0 fichero.vasp" >&2
	exit 1
fi
file="$1"
if [ ! -f "$file" ]; then
	echo "Error: no existe el fichero $file" >&2
	exit 1
fi

# Candidate ASE executables searched in order before falling back to PATH lookups.
ase_candidates=(
	"/home/salvador/miniforge3/bin/ase"
	"/home/salvador/.local/bin/ase"
)

# Write the minimal RASPA input files used to translate the CIF back into a relaxed structure.
raspa_input_file() {
	local target_dir="$1"
	cat >"${target_dir}/pseudo_atoms.def" <<'EOF'
#number of pseudo atoms
10
#type      print   as    chem  oxidation   mass        charge   polarization B-factor radii  connectivity anisotropic anisotropic-type   tinker-type
O1         yes     O1     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O2         yes     O2     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O3         yes     O3     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
Si         yes     Si    Si    0           28.0855     0.0      0.0          1.0      1.18   4            0           absolute           0
Ge         yes     Ge    Ge    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
F          yes     F     F     0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
C1         yes     C1    C1    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
C2         yes     C2    C2    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
N1         yes     N1    N1    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
H1         yes     H1    H1    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
EOF

	cat >"${target_dir}/simulation.input" <<'EOF'
SimulationType                MC
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
EOF
}

# Try a series of ASE entry points until one succeeds at producing a CIF snapshot of the VASP file.
run_ase_convert() {
	local src="$1"
	local dst="$2"
	local candidate

	for candidate in "${ase_candidates[@]}"; do
		if [ -x "$candidate" ]; then
			"$candidate" convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	done

	if command -v ase >/dev/null 2>&1; then
		candidate=$(command -v ase)
		"$candidate" convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
	fi

	for candidate in python3 python; do
		if command -v "$candidate" >/dev/null 2>&1; then
			"$candidate" -m ase.cli.main convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	done

	echo "Error: no se pudo ejecutar ASE para convertir $src a CIF" >&2
	return 1
}

# Run RASPA's "simulate" binary to relax the CIF and emit a final structure compatible with GULP.
run_raspa_convert() {
	local structure="$1"
	local workdir="$2"
	local output_cif
	local log_file

	if ! command -v simulate >/dev/null 2>&1; then
		echo "Error: no se encontró el binario simulate en PATH" >&2
		return 1
	fi

	output_cif="Movies/System_0/Framework_0_final_1_1_1_P1.cif"
	log_file="simulate.log"

	echo "Calculando la energía de ${structure} con simulate..."
	(
		cd "$workdir"
		if ! simulate >"$log_file" 2>&1; then
			echo "Error: simulate falló para ${structure}. Revisa ${workdir}/${log_file}" >&2
			exit 1
		fi
	) || return 1

	if [ ! -f "${workdir}/${output_cif}" ]; then
		echo "Error: simulate no generó $output_cif" >&2
		return 1
	fi

	cp "${workdir}/${output_cif}" "${structure}.cif"
	rm -f "${workdir}/${log_file}"
	rm -f "${workdir}/Local.cif" "${workdir}/pseudo_atoms.def" "${workdir}/simulation.input"
	rm -rf "${workdir}/Movies" "${workdir}/Restart" "${workdir}/Output" "${workdir}/VTK"
	return 0
}

# Main conversion pipeline: build auxiliary input, call ASE/RASPA, and emit the final *.gin package.
workspace=""
file_dir=$(dirname "$file")
workspace=$(mktemp -d "${file_dir}/vasp2gin_tmp.XXXXXX")
# Ensures the temporary workspace is removed on success and preserved on failure for debugging.
cleanup_workspace() {
	local status="$1"
	if [ -n "$workspace" ] && [ -d "$workspace" ]; then
		if [ "$status" -eq 0 ]; then
			rm -rf "$workspace"
		else
			echo "Aviso: se mantiene el directorio temporal para depuración: $workspace" >&2
		fi
	fi
}
trap 'cleanup_workspace $?' EXIT

raspa_input_file "$workspace"
rm -f "${workspace}/Local.cif"
run_ase_convert "$file" "${workspace}/Local.cif"
sed -i -e 's/_space_group_IT_number/_symmetry_Int_Tables_number/g' -e '/space_group_name_H-M_alt/d' "${workspace}/Local.cif"
run_raspa_convert "$file" "$workspace"

cif_path="${file}.cif"

if [ ! -f "$cif_path" ]; then
	echo "Error: no se encontró el CIF generado $cif_path" >&2
	exit 1
fi

read -r cell1 cell2 cell3 cell4 cell5 cell6 <<<"$(awk '
	/_cell_length_a/ {a=$2}
	/_cell_length_b/ {b=$2}
	/_cell_length_c/ {c=$2}
	/_cell_angle_alpha/ {al=$2}
	/_cell_angle_beta/ {be=$2}
	/_cell_angle_gamma/ {ga=$2}
	END {print a, b, c, al, be, ga}
' "$cif_path")"

if [ -z "$cell1" ]; then
	echo "Error: no se pudieron extraer los parámetros de celda de $cif_path" >&2
	exit 1
fi

gin_path="${file}.gin"

cat >"$gin_path" <<EOF
opti conp molmec
name $file
cell
 $cell1 $cell2 $cell3 $cell4 $cell5 $cell6
frac
EOF

sed -n '/^_atom_site_charge$/,$p' "$cif_path" | sed '/_atom_site_charge/d' | sed '/^[[:space:]]*$/d' | awk '{print $1 " core",$3,$4,$5}' >>"$gin_path"

cat >>"$gin_path" <<EOF
F  core 0.488968 0.476240 0.451806
F  core 0.996818 0.965041 0.422429
H1 core 0.410910 0.852615 0.163309
H1 core 0.573027 0.768482 0.426451
H1 core 0.543502 0.813574 0.599641
H1 core 0.404860 0.774465 0.449804
C1 core 0.510940 0.973266 0.474277
C1 core 0.437683 0.907460 0.242695
C2 core 0.503190 0.806491 0.474336
N1 core 0.486139 0.894237 0.401668
H1 core 0.913057 0.355881 0.186550
H1 core 0.907375 0.277853 0.475420
H1 core 0.926017 0.655554 0.485094
H1 core 0.954148 0.666324 0.297178
C1 core 0.934256 0.499336 0.241946
C1 core 0.939336 0.410718 0.266054
C2 core 0.989788 0.635305 0.410881
N1 core 0.988094 0.397409 0.424947
N1 core 0.980179 0.539264 0.386201
H1 core 0.902980 0.537821 0.136747
H1 core 0.631912 0.939609 0.700701
H1 core 0.976091 0.493288 0.720539
H1 core 1.073579 0.271293 0.446191
H1 core 1.049209 0.316198 0.622006
H1 core 1.094414 0.654046 0.461555
C1 core 1.012654 0.476428 0.497702
C2 core 1.005870 0.309461 0.497113
C2 core 1.060695 0.493015 0.667327
H1 core 1.110789 0.557246 0.688064
H1 core 1.132846 0.442661 0.724544
H1 core 0.429305 1.151864 0.467072
H1 core 0.447840 1.163373 0.275079
H1 core 0.594080 1.150953 0.432599
C1 core 0.433031 0.996066 0.218436
C2 core 0.488720 1.132103 0.387729
N1 core 0.478728 1.036054 0.362711
C2 core 0.559366 0.989840 0.643851
H1 core 0.401975 1.034503 0.113158
H1 core 0.609122 1.054186 0.664500
H1 core 0.475051 0.989761 0.697534

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
F  core  F  core
F  shel  F  shel
C2 core  C2 core
C1 core  C1 core
H1 core  H1 core
N1 core  N1 core
end

connect 75 80 single  0 0 0
connect 76 81 single  0 0 0
connect 77 81 single  0 0 0
connect 78 81 single  0 0 0
connect 79 82 resonant  0 0 0
connect 79 108 resonant  0 0 0
connect 79 109 single  0 0 0
connect 80 82 resonant  0 0 0
connect 80 106 resonant  0 0 0
connect 81 82 single  0 0 0
connect 83 88 single  0 0 0
connect 84 99 single  0 0 0
connect 85 89 single  0 0 0
connect 86 89 single  0 0 0
connect 87 88 resonant  0 0 0
connect 87 91 resonant  0 0 0
connect 87 92 single  0 0 0
connect 88 90 resonant  0 0 0
connect 89 91 single  0 0 0
connect 89 97 single  0 0 0
connect 90 98 resonant  0 0 0
connect 90 99 single  0 0 0
connect 91 98 resonant  0 0 0
connect 93 109 single  0 0 0
connect 94 100 single  0 0 0
connect 95 99 single  0 0 0
connect 96 99 single  0 0 0
connect 98 100 single  0 0 0
connect 100 101 single  0 0 0
connect 100 102 single  0 0 0
connect 103 107 single  0 0 0
connect 104 107 single  0 0 0
connect 105 107 single  0 0 0
connect 106 108 resonant  0 0 0
connect 106 110 single  0 0 0
connect 107 108 single  0 0 0
connect 109 111 single  0 0 0
connect 109 112 single  0 0 0

library Germanate_OSDA
switch_minimiser rfo gnorm 0.01
stepmx opt 0.05
cuts 1.0
dump every 1 $file.grs
output cif $file.cif
EOF

cat >Germanate_OSDA.lib <<'EOF'
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
C2    core  -0.3
C1    core  -0.1
H1    core   0.1
N1    core   0.55

epsilon/sigma
H1    core 0.659139E-03  3.19500
C2    core 0.412396E-02  3.89830
C1    core 0.412396E-02  3.89830
N1    core 0.335641E-02  3.66210
O1    shel 0.414998E-02  3.40460
O2    shel 0.414998E-02  3.40460
O3    shel 0.414998E-02  3.40460
F     shel 0.314392E-02  3.47200
Si    core 0.134430E-01  4.27000
Ge    core 0.173458E-01  4.27000
O1    core 0   3.40460
O2    core 0   3.40460
O3    core 0   3.40460
F     core 0   3.47200

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

harmonic intra bond
H1    core H1    core  30.3551     0.650000      0.00000
H1    core C2    core  30.3551      1.09000      0.00000
C2    core C2    core  30.3551      1.53000      0.00000
H1    core C1    core  30.3551      1.02000      0.00000
C1    core C2    core  30.3551      1.46000      0.00000
C1    core C1    core  45.5327      1.39000      0.00000
H1    core N1    core  30.3551     0.970000      0.00000
C2    core N1    core  30.3551      1.41000      0.00000
C1    core N1    core  45.5327      1.34000      0.00000
N1    core N1    core  45.5327      1.29000      0.00000

lenn epsilon geometric 12  6 x13
H1    core H1    core  0.000 12.500
H1    core C2    core  0.000 12.500
C2    core C2    core  0.000 12.500
H1    core C1    core  0.000 12.500
C1    core C2    core  0.000 12.500
C1    core C1    core  0.000 12.500
H1    core N1    core  0.000 12.500
C2    core N1    core  0.000 12.500
C1    core N1    core  0.000 12.500
N1    core N1    core  0.000 12.500
H1    core O1    shel  0.000 12.500
C2    core O1    shel  0.000 12.500
C1    core O1    shel  0.000 12.500
N1    core O1    shel  0.000 12.500
H1    core O2    shel  0.000 12.500
C2    core O2    shel  0.000 12.500
C1    core O2    shel  0.000 12.500
N1    core O2    shel  0.000 12.500
H1    core O3    shel  0.000 12.500
C2    core O3    shel  0.000 12.500
C1    core O3    shel  0.000 12.500
N1    core O3    shel  0.000 12.500
H1    core F     shel  0.000 12.500
C2    core F     shel  0.000 12.500
C1    core F     shel  0.000 12.500
N1    core F     shel  0.000 12.500
H1    core Si    core  0.000 12.500
C2    core Si    core  0.000 12.500
C1    core Si    core  0.000 12.500
N1    core Si    core  0.000 12.500
three bond intra regular
C2    core X     core X     core 4.3364     109.47
C1    core X     core X     core 4.3364     120.00
N1    core X     core X     core 4.3364     120.00
torsion bond intra single    dreiding
X     cor C2    cor C2    cor X     cor  0.43364E-01  -3 180.00
C2    cor C1    cor C2    cor X     cor  0.43364E-01  -3 180.00
C1    cor C1    cor C2    cor X     cor  0.21682E-01  -6 0.0000
N1    cor C1    cor C2    cor X     cor  0.21682E-01  -6 0.0000
H1    cor C1    cor C2    cor X     cor  0.43364E-01  -3 180.00
X     cor C1    cor C1    cor X     cor  0.54206      -2 0.0000
C2    cor N1    cor C2    cor X     cor  0.43364E-01  -3 180.00
C1    cor N1    cor C2    cor X     cor  0.21682E-01  -6 0.0000
N1    cor N1    cor C2    cor X     cor  0.21682E-01  -6 0.0000
H1    cor N1    cor C2    cor X     cor  0.43364E-01  -3 180.00
X     cor N1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor N1    cor X     cor  0.54206      -2 0.0000
torsion bond intra double    dreiding
X     cor C1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor N1    cor X     cor  0.54206      -2 0.0000
torsion bond intra resonant  dreiding
X     cor C1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor N1    cor X     cor  0.54206      -2 0.0000
torsion bond intra resonant  exocyclic dreiding
X     cor C1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor C1    cor X     cor  0.54206      -2 0.0000
X     cor N1    cor N1    cor X     cor  0.54206      -2 0.0000
torsion bond intra single    exocyclic dreiding
X     cor C1    cor C1    cor X     cor  0.21682      -2 0.0000
X     cor N1    cor C1    cor X     cor  0.21682      -2 0.0000
X     cor N1    cor N1    cor X     cor  0.21682      -2 0.0000
inversion bond intra only3
C1    cor X     cor X     cor X     cor 1.7346
N1    cor X     cor X     cor X     cor 1.7346

rtol 1.1
element
cova H     0.420
cova C     0.770
cova N     0.702
cova O     0.0
cova F     0.0
cova Ge    0.0
cova Si    0.0
end
cuts 1.0
EOF

exit 0
