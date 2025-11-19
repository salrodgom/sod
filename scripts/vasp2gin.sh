#!/bin/bash
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
5
#type      print   as    chem  oxidation   mass        charge   polarization B-factor radii  connectivity anisotropic anisotropic-type   tinker-type
O1         yes     O1     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O2         yes     O2     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
O3         yes     O3     O    0           15.9994     0.0      0.0          1.0      0.5    2            0           absolute           0
Si         yes     Si    Si    0           28.0855     0.0      0.0          1.0      1.18   4            0           absolute           0
Ge         yes     Ge    Ge    0           26.981539   0.0      0.0          1.0      1.18   4            0           absolute           0
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
opti conp prop
name $file
cell
 $cell1 $cell2 $cell3 $cell4 $cell5 $cell6
frac
EOF

sed -n '/^_atom_site_charge$/,$p' "$cif_path" | sed '/_atom_site_charge/d' | sed '/^[[:space:]]*$/d' | awk '{print $1 " core",$3,$4,$5}' >>"$gin_path"

cat >>"$gin_path" <<EOF

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
output cif $file.cif
EOF

cat >Germanate.lib <<'EOF'
# Force-field description tailored to the Si/Ge germanate framework used during calibration.
# O1: Ge-O1-Ge
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
EOF

exit 0
