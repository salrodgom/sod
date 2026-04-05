#!/bin/bash
#*******************************************************************************
#    Copyright (c) 2026, Salvador R.G. Balestra
#
#    This file is part of the SOD package.
#
#    SOD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#*******************************************************************************

# Unified VASP->GULP driver.
# Protocol 1.0 reproduces the classic vasp2gin conversion only.
# Protocol 2.0 runs the staged GULP workflow and promotes the final outputs.
set -euo pipefail
shopt -s nullglob

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

ase_candidates=(
	"/home/salvador/miniforge3/bin/ase"
	"/home/salvador/.local/bin/ase"
)

default_shell_cuts="${SOD_PROTOCOL_SHELL_CUTS:-1.4}"
default_recovery_cuts="${SOD_PROTOCOL_RECOVERY_CUTS:-1.5}"
default_shell_maxcyc="${SOD_PROTOCOL_SHELL_MAXCYC:-100}"
default_cell_maxcyc="${SOD_PROTOCOL_CELL_MAXCYC:-600}"
default_recovery_maxcyc="${SOD_PROTOCOL_RECOVERY_MAXCYC:-400}"

default_shell_stepmx_opt="${SOD_PROTOCOL_SHELL_STEPMX_OPT:-default}"
default_shell_stepmx_rfo="${SOD_PROTOCOL_SHELL_STEPMX_RFO:-default}"
default_shell_switch_stepmx="${SOD_PROTOCOL_SHELL_SWITCH_STEPMX:-default}"
default_shell_switch_gnorm="${SOD_PROTOCOL_SHELL_SWITCH_GNORM:-0.05}"

default_cell_stepmx_opt="${SOD_PROTOCOL_CELL_STEPMX_OPT:-${SOD_PROTOCOL_INTERNAL_STEPMX_OPT:-default}}"
default_cell_stepmx_rfo="${SOD_PROTOCOL_CELL_STEPMX_RFO:-${SOD_PROTOCOL_INTERNAL_STEPMX_RFO:-default}}"
default_cell_switch_stepmx="${SOD_PROTOCOL_CELL_SWITCH_STEPMX:-${SOD_PROTOCOL_INTERNAL_SWITCH_STEPMX:-default}}"
default_cell_switch_gnorm="${SOD_PROTOCOL_CELL_SWITCH_GNORM:-${SOD_PROTOCOL_INTERNAL_SWITCH_GNORM:-0.05}}"

default_recovery_stepmx_opt="${SOD_PROTOCOL_RECOVERY_STEPMX_OPT:-default}"
default_recovery_stepmx_rfo="${SOD_PROTOCOL_RECOVERY_STEPMX_RFO:-default}"
default_recovery_switch_stepmx="${SOD_PROTOCOL_RECOVERY_SWITCH_STEPMX:-default}"
default_recovery_switch_gnorm="${SOD_PROTOCOL_RECOVERY_SWITCH_GNORM:-0.02}"

print_usage() {
	cat <<'EOF'
Usage:
  vasp2gin.sh [--protocol 1.0|2.0] structure.vasp [template_payload_or_include]
  vasp2gin.sh [--protocol 2.0] [file1.vasp file2.vasp ...]

Protocols:
  1.0  Classic converter. Creates structure.vasp.cif, structure.vasp.gin, and Germanate_OSDA.lib.
  2.0  Default. Runs the staged GULP protocol and promotes the final *.gin/*.gout/*.grs/*.cif.

Options:
  --protocol <ver>    Select protocol 1.0 or 2.0 [default: 2.0]
  --protocole <ver>   Alias for --protocol
  -h, --help          Show this help

Environment used by protocol 2.0:
  SOD_PROTOCOL_TEMPLATE_GIN
  SOD_PROTOCOL_LIBRARY_FILE
  SOD_PROTOCOL_TEMPLATE_PAYLOAD
  SOD_GULP_SLOT_LABEL
  SOD_PROTOCOL_SHELL_CUTS
  SOD_PROTOCOL_RECOVERY_CUTS
  SOD_PROTOCOL_SHELL_MAXCYC
  SOD_PROTOCOL_CELL_MAXCYC
  SOD_PROTOCOL_RECOVERY_MAXCYC
  SOD_PROTOCOL_*_STEPMX_*
EOF
}

normalize_protocol_version() {
	local value="${1:-}"
	case "${value,,}" in
		1|1.0) printf '1.0\n' ;;
		2|2.0|"") printf '2.0\n' ;;
		*)
			echo "Error: unsupported protocol version '${value}'. Use 1.0 or 2.0." >&2
			return 1
			;;
	esac
}

resolve_existing_path() {
	local path="${1:-}"
	local dir_part=""
	local base_part=""

	[ -n "$path" ] || {
		printf '\n'
		return 0
	}

	if [ ! -e "$path" ]; then
		echo "Error: path ${path} does not exist." >&2
		return 1
	fi

	dir_part=$(dirname "$path")
	base_part=$(basename "$path")
	printf '%s/%s\n' "$(CDPATH= cd -- "$dir_part" && pwd)" "$base_part"
}

ase_candidate_works() {
	local candidate="$1"

	[ -n "$candidate" ] || return 1
	[ -x "$candidate" ] || return 1
	"$candidate" --version >/dev/null 2>&1
}

run_ase_convert() {
	local src="$1"
	local dst="$2"
	local candidate=""

	if command -v ase >/dev/null 2>&1; then
		candidate=$(command -v ase)
		if ase_candidate_works "$candidate"; then
			"$candidate" convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	fi

	for candidate in "${ase_candidates[@]}"; do
		if ase_candidate_works "$candidate"; then
			"$candidate" convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	done

	for candidate in python3 python; do
		if command -v "$candidate" >/dev/null 2>&1; then
			"$candidate" -m ase.cli.main convert -i vasp -o cif "$src" "$dst" && [ -f "$dst" ] && return 0
		fi
	done

	echo "Error: ASE could not convert $src into a CIF file." >&2
	return 1
}

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

run_raspa_convert() {
	local structure="$1"
	local workdir="$2"
	local output_cif="Movies/System_0/Framework_0_final_1_1_1_P1.cif"
	local log_file="simulate.log"

	if ! command -v simulate >/dev/null 2>&1; then
		echo "Error: the simulate binary was not found in PATH." >&2
		return 1
	fi

	(
		cd "$workdir"
		if ! simulate >"$log_file" 2>&1; then
			echo "Error: simulate failed for ${structure}. Check ${workdir}/${log_file}" >&2
			exit 1
		fi
	) || return 1

	if [ ! -f "${workdir}/${output_cif}" ]; then
		echo "Error: simulate did not generate ${output_cif}" >&2
		return 1
	fi

	cp "${workdir}/${output_cif}" "${structure}.cif"
	rm -f "${workdir}/${log_file}"
	rm -f "${workdir}/Local.cif" "${workdir}/pseudo_atoms.def" "${workdir}/simulation.input"
	rm -rf "${workdir}/Movies" "${workdir}/Restart" "${workdir}/Output" "${workdir}/VTK"
}

write_legacy_library() {
	local library_path="$1"

	cat >"${library_path}" <<'EOF'
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
cuts 1.2
EOF
}

run_protocol_1() {
	local file="${1:-}"
	local template_include="${2:-}"
	local file_dir=""
	local workspace=""
	local cif_path=""
	local gin_path=""
	local library_path=""
	local cell1="" cell2="" cell3="" cell4="" cell5="" cell6=""

	[ -n "$file" ] || {
		echo "Error: protocol 1.0 expects structure.vasp [template.include]" >&2
		return 1
	}
	[ -f "$file" ] || {
		echo "Error: missing VASP file ${file}" >&2
		return 1
	}
	if [ -n "$template_include" ] && [ ! -f "$template_include" ]; then
		echo "Error: template file not found: $template_include" >&2
		return 1
	fi

	file_dir=$(CDPATH= cd -- "$(dirname -- "$file")" && pwd)
	workspace=$(mktemp -d "${file_dir}/vasp2gin_tmp.XXXXXX")
	raspa_input_file "$workspace"
	run_ase_convert "$file" "${workspace}/Local.cif"
	sed -i -e 's/_space_group_IT_number/_symmetry_Int_Tables_number/g' -e '/space_group_name_H-M_alt/d' "${workspace}/Local.cif"
	run_raspa_convert "$file" "$workspace"

	cif_path="${file}.cif"
	[ -f "$cif_path" ] || {
		echo "Error: no typed CIF was generated at ${cif_path}" >&2
		return 1
	}

	read -r cell1 cell2 cell3 cell4 cell5 cell6 <<<"$(awk '
		/_cell_length_a/ {a=$2}
		/_cell_length_b/ {b=$2}
		/_cell_length_c/ {c=$2}
		/_cell_angle_alpha/ {al=$2}
		/_cell_angle_beta/ {be=$2}
		/_cell_angle_gamma/ {ga=$2}
		END {print a, b, c, al, be, ga}
	' "$cif_path")"
	[ -n "$cell1" ] || {
		echo "Error: failed to extract cell parameters from ${cif_path}" >&2
		return 1
	}

	gin_path="${file}.gin"
	cat >"$gin_path" <<EOF
opti conp molmec
name $file
cell
 $cell1 $cell2 $cell3 $cell4 $cell5 $cell6
frac
EOF
	sed -n '/^_atom_site_charge$/,$p' "$cif_path" | sed '/_atom_site_charge/d' | sed '/^[[:space:]]*$/d' | \
		awk '{print $1 " core",$3,$4,$5}' >>"$gin_path"
	if [ -n "$template_include" ]; then
		cat "$template_include" >>"$gin_path"
	fi
	cat >>"$gin_path" <<EOF
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

library Germanate_OSDA nodump
switch_minimiser rfo gnorm 0.01
stepmx opt 0.05
stepmx rfo 0.02
switch_stepmx 0.01 gnorm 0.05
cuts 1.0
dump every 1 $file.grs
output cif $file.cif
EOF

	library_path="${file_dir}/Germanate_OSDA.lib"
	write_legacy_library "$library_path"
	rm -rf "$workspace"
	return 0
}

protocol_detect_template_gin() {
	local template_override="${SOD_PROTOCOL_TEMPLATE_GIN:-}"

	if [ -n "$template_override" ]; then
		printf '%s\n' "$template_override"
		return 0
	fi
	if [ -f "./protocol_template.gin" ]; then
		printf '%s\n' "./protocol_template.gin"
		return 0
	fi
	if [ -f "${script_dir}/protocol_template.gin" ]; then
		printf '%s\n' "${script_dir}/protocol_template.gin"
		return 0
	fi
	if [ -f "${script_dir}/reference/protocol_template.gin" ]; then
		printf '%s\n' "${script_dir}/reference/protocol_template.gin"
		return 0
	fi

	echo "Error: no protocol template .gin was found in the working directory or under scripts/." >&2
	return 1
}

protocol_detect_library_file() {
	local lib_override="${SOD_PROTOCOL_LIBRARY_FILE:-}"

	if [ -n "$lib_override" ]; then
		printf '%s\n' "$lib_override"
		return 0
	fi
	if [ -f "./framework_template.lib" ]; then
		printf '%s\n' "./framework_template.lib"
		return 0
	fi
	if [ -f "${script_dir}/framework_template.lib" ]; then
		printf '%s\n' "${script_dir}/framework_template.lib"
		return 0
	fi

	echo "Error: no library file was found in the working directory or under scripts/." >&2
	return 1
}

protocol_detect_template_payload() {
	local payload_override="${1:-}"
	local template_override="${SOD_PROTOCOL_TEMPLATE_PAYLOAD:-}"

	if [ -n "$payload_override" ]; then
		printf '%s\n' "$payload_override"
		return 0
	fi
	if [ -n "$template_override" ]; then
		printf '%s\n' "$template_override"
		return 0
	fi
	if [ -f "./template_payload.gin" ]; then
		printf '%s\n' "./template_payload.gin"
		return 0
	fi
	if [ -f "./default_template.gin" ]; then
		printf '%s\n' "./default_template.gin"
		return 0
	fi
	if [ -f "./default_template.include" ]; then
		printf '%s\n' "./default_template.include"
		return 0
	fi
	if [ -f "${script_dir}/default_template.include" ]; then
		printf '%s\n' "${script_dir}/default_template.include"
		return 0
	fi

	printf '\n'
	return 0
}

extract_cell_from_cif() {
	local cif_path="$1"
	awk '
		/_cell_length_a/ {a=$2}
		/_cell_length_b/ {b=$2}
		/_cell_length_c/ {c=$2}
		/_cell_angle_alpha/ {al=$2}
		/_cell_angle_beta/ {be=$2}
		/_cell_angle_gamma/ {ga=$2}
		END {print a, b, c, al, be, ga}
	' "$cif_path"
}

append_cif_coords() {
	local cif_path="$1"
	local target_gin="$2"
	sed -n '/^_atom_site_charge$/,$p' "$cif_path" | sed '/_atom_site_charge/d' | sed '/^[[:space:]]*$/d' | \
		awk '{print $1 " core", $3, $4, $5}' >>"$target_gin"
}

append_template_payload_block() {
	local payload_gin="$1"
	local target_gin="$2"

	[ -n "$payload_gin" ] || return 1
	[ -f "$payload_gin" ] || return 1
	cat "$payload_gin" >>"$target_gin"
	return 0
}

count_gin_coords() {
	local gin_file="$1"
	awk '
		tolower($1)=="frac" || tolower($1)=="fractional" { in_coords = 1; next }
		in_coords && tolower($1)=="species" { in_coords = 0 }
		in_coords && NF >= 5 && tolower($2)=="core" { count++ }
		END { print count + 0 }
	' "$gin_file"
}

append_template_middle_block() {
	local template_gin="$1"
	local target_gin="$2"
	local coord_count="$3"

	awk -v coord_count="$coord_count" '
		BEGIN {
			after_coords = 0
			saw_explicit_coords = 0
		}
		tolower($1)=="frac" || tolower($1)=="fractional" {
			saw_explicit_coords = 1
			after_coords = 1
			next
		}
		!after_coords {
			if (!saw_explicit_coords) {
				if (NF >= 2 && (tolower($2)=="core" || tolower($2)=="shel")) {
					next
				}
				if (tolower($1)=="space" || tolower($1)=="connect") {
					after_coords = 1
				} else {
					next
				}
			} else {
				next
			}
		}
		tolower($1)=="species" { exit }
		NF >= 2 && (tolower($2)=="core" || tolower($2)=="shel") { next }
		tolower($1)=="connect" {
			if ($2 > coord_count || $3 > coord_count) {
				next
			}
		}
		{ print }
	' "$template_gin" >>"$target_gin"
}

select_middle_block_source() {
	local template_gin="$1"
	local payload_gin="$2"

	if [ -n "$payload_gin" ] && [ -f "$payload_gin" ]; then
		printf '%s\n' "$payload_gin"
	else
		printf '%s\n' "$template_gin"
	fi
}

append_template_species_block() {
	local template_gin="$1"
	local target_gin="$2"
	awk '
		/^[[:space:]]*species([[:space:]]|$)/ { in_species = 1 }
		in_species { print }
		in_species && /^[[:space:]]*end([[:space:]]|$)/ { exit }
	' "$template_gin" >>"$target_gin"
}

is_default_control_value() {
	local value="${1:-}"
	case "${value,,}" in
		""|default|defaults|auto|omit|none) return 0 ;;
	esac
	return 1
}

append_stage_controls() {
	local target_gin="$1"
	local stage_base="$2"
	local library_name="$3"
	local maxcyc_value="$4"
	local cuts_value="$5"
	local stepmx_opt="$6"
	local stepmx_rfo="$7"
	local switch_stepmx_value="$8"
	local switch_stepmx_gnorm="$9"

	{
		printf '\n'
		printf 'library %s nodump\n' "$library_name"
		printf 'maxcyc opt %s\n' "$maxcyc_value"
		printf 'switch_minimiser rfo gnorm 0.01\n'
		if ! is_default_control_value "$stepmx_opt"; then
			printf 'stepmx opt %s\n' "$stepmx_opt"
		fi
		if ! is_default_control_value "$stepmx_rfo"; then
			printf 'stepmx rfo %s\n' "$stepmx_rfo"
		fi
		if ! is_default_control_value "$switch_stepmx_value"; then
			printf 'switch_stepmx %s gnorm %s\n' "$switch_stepmx_value" "$switch_stepmx_gnorm"
		fi
		printf 'cuts %s\n' "$cuts_value"
		printf 'dump every 1 %s.grs\n' "$stage_base"
		printf 'output cif %s.cif\n' "$stage_base"
	} >>"$target_gin"
}

create_stage_from_cif() {
	local input_file="$1"
	local cif_path="$2"
	local template_gin="$3"
	local payload_gin="$4"
	local stage_base="$5"
	local keyword_line="$6"
	local library_name="$7"
	local maxcyc_value="$8"
	local cuts_value="$9"
	local stepmx_opt="${10}"
	local stepmx_rfo="${11}"
	local switch_stepmx_value="${12}"
	local switch_stepmx_gnorm="${13}"
	local target_gin="${stage_base}.gin"
	local cell1="" cell2="" cell3="" cell4="" cell5="" cell6=""
	local coord_count=0
	local middle_source=""

	read -r cell1 cell2 cell3 cell4 cell5 cell6 <<<"$(extract_cell_from_cif "$cif_path")"
	[ -n "${cell1:-}" ] || {
		echo "Error: failed to extract cell parameters from ${cif_path}" >&2
		return 1
	}

	cat >"$target_gin" <<EOF
${keyword_line}
name ${stage_base}
cell
 ${cell1} ${cell2} ${cell3} ${cell4} ${cell5} ${cell6}
frac
EOF
	append_cif_coords "$cif_path" "$target_gin"
	if ! append_template_payload_block "$payload_gin" "$target_gin"; then
		coord_count=$(count_gin_coords "$target_gin")
		middle_source=$(select_middle_block_source "$template_gin" "$payload_gin")
		printf '\n' >>"$target_gin"
		append_template_middle_block "$middle_source" "$target_gin" "$coord_count"
		printf '\n' >>"$target_gin"
	fi
	append_template_species_block "$template_gin" "$target_gin"
	append_stage_controls "$target_gin" "$stage_base" "$library_name" "$maxcyc_value" "$cuts_value" \
		"$stepmx_opt" "$stepmx_rfo" "$switch_stepmx_value" "$switch_stepmx_gnorm"
}

create_stage_from_grs() {
	local source_grs="$1"
	local template_gin="$2"
	local payload_gin="$3"
	local stage_base="$4"
	local keyword_line="$5"
	local coord_mode="$6"
	local library_name="$7"
	local maxcyc_value="$8"
	local cuts_value="$9"
	local stepmx_opt="${10}"
	local stepmx_rfo="${11}"
	local switch_stepmx_value="${12}"
	local switch_stepmx_gnorm="${13}"
	local target_gin="${stage_base}.gin"
	local cell_line=""
	local fractional_line=""
	local coord_count=0
	local middle_source=""

	[ -f "$source_grs" ] || {
		echo "Error: missing source restart file ${source_grs}" >&2
		return 1
	}

	cell_line=$(awk 'tolower($1)=="cell" { getline; if (NF) { print; exit } }' "$source_grs")
	fractional_line=$(awk 'tolower($1)=="fractional" { print; exit }' "$source_grs")
	[ -n "$cell_line" ] && [ -n "$fractional_line" ] || {
		echo "Error: failed to reconstruct header from ${source_grs}" >&2
		return 1
	}

	cat >"$target_gin" <<EOF
${keyword_line}
name ${stage_base}
cell
${cell_line}
${fractional_line}
EOF
	awk -v coord_mode="$coord_mode" '
		tolower($1)=="fractional" { in_coords = 1; next }
		in_coords {
			if (tolower($1)=="totalenergy" || tolower($1)=="species") {
				exit
			}
			if (NF >= 5 && coord_mode == "all" && (tolower($2)=="core" || tolower($2)=="shel")) {
				print $1, $2, $3, $4, $5
			}
			if (NF >= 5 && coord_mode == "core_only" && tolower($2)=="core") {
				print $1, $2, $3, $4, $5
			}
		}
	' "$source_grs" >>"$target_gin"
	coord_count=$(count_gin_coords "$target_gin")
	middle_source=$(select_middle_block_source "$template_gin" "$payload_gin")
	printf '\n' >>"$target_gin"
	append_template_middle_block "$middle_source" "$target_gin" "$coord_count"
	printf '\n' >>"$target_gin"
	append_template_species_block "$template_gin" "$target_gin"
	append_stage_controls "$target_gin" "$stage_base" "$library_name" "$maxcyc_value" "$cuts_value" \
		"$stepmx_opt" "$stepmx_rfo" "$switch_stepmx_value" "$switch_stepmx_gnorm"
}

extract_energy_value() {
	local gout_file="$1"
	local energy=""

	energy=$(awk '/l e/ { value = $4 } END { if (value != "") print value }' "$gout_file")
	if [ -z "$energy" ]; then
		energy=$(awk '/Total lattice energy[[:space:]]*=/ && $6 == "eV" { value = $5 } END { if (value != "") print value }' "$gout_file")
	fi
	if [ -z "$energy" ]; then
		energy=$(awk '/Final energy[[:space:]]*=/ { value = $4 } END { if (value != "") print value }' "$gout_file")
	fi
	if [[ "$energy" =~ ^[-+]?[0-9]+([.][0-9]+)?([Ee][-+]?[0-9]+)?$ ]]; then
		printf '%s\n' "$energy"
		return 0
	fi
	return 1
}

is_failed_output() {
	local gout_file="$1"
	grep -qF 'Conditions for a minimum have not been satisfied.' "$gout_file" || \
	grep -qF 'Too many failed attempts to optimise' "$gout_file" || \
	grep -qF 'Largest core-shell distance exceeds cutoff of cuts' "$gout_file"
}

extract_last_gnorm() {
	local gout_file="$1"
	local gnorm=""

	gnorm=$(awk '
		/Gnorm:/ {
			for (i = 1; i < NF; i++) {
				if ($i == "Gnorm:") {
					value = $(i + 1)
				}
			}
		}
		END { if (value != "") print value }
	' "$gout_file")
	if [[ "$gnorm" =~ ^[-+]?[0-9]+([.][0-9]+)?([Ee][-+]?[0-9]+)?$ ]]; then
		printf '%s\n' "$gnorm"
		return 0
	fi
	return 1
}

is_complete_output() {
	local gout_file="$1"
	[ -f "$gout_file" ] || return 1
	is_failed_output "$gout_file" && return 1
	extract_energy_value "$gout_file" >/dev/null
}

select_better_incomplete_output() {
	local current_gout="${1:-}"
	local candidate_gout="${2:-}"
	local current_energy=""
	local candidate_energy=""
	local current_gnorm=""
	local candidate_gnorm=""
	local decision=""

	if [ -z "$candidate_gout" ] || [ ! -f "$candidate_gout" ]; then
		printf '%s\n' "$current_gout"
		return 0
	fi
	if [ -z "$current_gout" ] || [ ! -f "$current_gout" ]; then
		printf '%s\n' "$candidate_gout"
		return 0
	fi
	if ! current_energy=$(extract_energy_value "$current_gout"); then
		printf '%s\n' "$candidate_gout"
		return 0
	fi
	if ! candidate_energy=$(extract_energy_value "$candidate_gout"); then
		printf '%s\n' "$current_gout"
		return 0
	fi

	current_gnorm=$(extract_last_gnorm "$current_gout" 2>/dev/null || true)
	candidate_gnorm=$(extract_last_gnorm "$candidate_gout" 2>/dev/null || true)

	if [ -n "$current_gnorm" ] && [ -n "$candidate_gnorm" ]; then
		decision=$(awk -v ga="$current_gnorm" -v gb="$candidate_gnorm" -v ea="$current_energy" -v eb="$candidate_energy" '
			BEGIN {
				gtol = 1.0e-8
				etol = 1.0e-10
				if (ga < gb - gtol) print "current"
				else if (gb < ga - gtol) print "candidate"
				else if (ea < eb - etol) print "current"
				else if (eb < ea - etol) print "candidate"
				else print "current"
			}
		')
	elif [ -n "$current_gnorm" ]; then
		decision="current"
	elif [ -n "$candidate_gnorm" ]; then
		decision="candidate"
	else
		decision=$(awk -v ea="$current_energy" -v eb="$candidate_energy" '
			BEGIN {
				tol = 1.0e-10
				if (ea <= eb + tol) print "current"
				else print "candidate"
			}
		')
	fi

	if [ "$decision" = "candidate" ]; then
		printf '%s\n' "$candidate_gout"
	else
		printf '%s\n' "$current_gout"
	fi
}

promote_protocol_result() {
	local canonical_base="$1"
	local selected_stage="$2"

	[ -n "$selected_stage" ] || return 1
	[ -f "${selected_stage}.gin" ] && cp "${selected_stage}.gin" "${canonical_base}.gin"
	[ -f "${selected_stage}.gout" ] && cp "${selected_stage}.gout" "${canonical_base}.gout"
	[ -f "${selected_stage}.grs" ] && cp "${selected_stage}.grs" "${canonical_base}.grs"
	[ -f "${selected_stage}.cif" ] && cp "${selected_stage}.cif" "${canonical_base}.cif"
}

archive_existing_workspace() {
	local workspace_path="$1"
	local archived_path=""
	local suffix=0

	[ -e "$workspace_path" ] || return 0
	archived_path="${workspace_path}.previous_$(date +%Y%m%d_%H%M%S)"
	while [ -e "$archived_path" ]; do
		suffix=$((suffix + 1))
		archived_path="${workspace_path}.previous_$(date +%Y%m%d_%H%M%S)_${suffix}"
	done
	mv "$workspace_path" "$archived_path"
	echo "[protocol] Preserved a previous workspace at ${archived_path}" >&2
}

select_workspace_path() {
	local file_dir="$1"
	local file_name="$2"
	local slot_label="${SOD_GULP_SLOT_LABEL:-}"

	if [ -n "$slot_label" ]; then
		printf '%s/%s\n' "$file_dir" "$slot_label"
	else
		printf '%s/%s\n' "$file_dir" "protocol_${file_name%.vasp}"
	fi
}

report_stage_status() {
	local stage_base="$1"
	local gout_file="${stage_base}.gout"
	local energy=""

	if [ ! -f "$gout_file" ]; then
		echo "[protocol] ${stage_base}: missing output file." >&2
		return
	fi
	if energy=$(extract_energy_value "$gout_file"); then
		echo "[protocol] ${stage_base}: energy = ${energy} eV"
	else
		echo "[protocol] ${stage_base}: no numeric final energy was found." >&2
	fi
	if is_failed_output "$gout_file"; then
		echo "[protocol] ${stage_base}: failure markers detected." >&2
	fi
}

run_stage() {
	local stage_base="$1"
	local gin_file="${stage_base}.gin"
	local gout_file="${stage_base}.gout"

	[ -f "$gin_file" ] || {
		echo "Error: missing input file ${gin_file}" >&2
		return 1
	}
	echo "[protocol] Running ${gin_file}"
	: >"$gout_file"
	{
		printf '# GULP stage: %s\n' "$stage_base"
		printf '# Input file: %s\n\n' "$gin_file"
	} >>"$gout_file"
	if ! gulp <"$gin_file" >>"$gout_file"; then
		echo "[protocol] GULP returned a non-zero status for ${gin_file}" >&2
	fi
	report_stage_status "$stage_base"
}

process_protocol_vasp() {
	local file="$1"
	local template_gin="$2"
	local library_file="$3"
	local payload_gin="$4"
	local file_dir=""
	local file_name=""
	local workspace=""
	local library_name=""
	local cif_path=""
	local stage1_base=""
	local stage2_base=""
	local stage3_base=""
	local final_stage=""
	local best_incomplete_stage=""
	local canonical_base=""
	local final_complete=0

	[ -f "$file" ] || {
		echo "Error: missing VASP file ${file}" >&2
		return 1
	}
	[ -f "$library_file" ] || {
		echo "Error: missing library file ${library_file}" >&2
		return 1
	}

	file_dir=$(CDPATH= cd -- "$(dirname -- "$file")" && pwd)
	file_name=$(basename "$file")
	library_name=$(basename "$library_file")
	canonical_base="${file_dir}/${file_name}"

	(
		workspace=$(select_workspace_path "$file_dir" "$file_name")
		archive_existing_workspace "$workspace"
		mkdir -p "$workspace"
		cp "$file" "${workspace}/${file_name}"
		cp "$library_file" "${workspace}/${library_name}"

		cd "$workspace"
		raspa_input_file "."
		run_ase_convert "$file_name" "Local.cif"
		sed -i -e 's/_space_group_IT_number/_symmetry_Int_Tables_number/g' -e '/space_group_name_H-M_alt/d' "Local.cif"
		run_raspa_convert "$file_name" "."

		cif_path="${file_name}.cif"
		[ -f "$cif_path" ] || {
			echo "Error: typed CIF ${cif_path} was not created." >&2
			exit 1
		}

		stage1_base="${file_name}.protocol01_shell"
		stage2_base="${file_name}.protocol02_cell"
		stage3_base="${file_name}.protocol03_shellfree"

		create_stage_from_cif "$file_name" "$cif_path" "$template_gin" "$payload_gin" "$stage1_base" "opti conv shell molmec" "$library_name" \
			"$default_shell_maxcyc" "$default_shell_cuts" "$default_shell_stepmx_opt" "$default_shell_stepmx_rfo" \
			"$default_shell_switch_stepmx" "$default_shell_switch_gnorm"
		run_stage "$stage1_base"
		[ -f "${stage1_base}.grs" ] || {
			echo "Error: stage 1 did not produce ${stage1_base}.grs" >&2
			exit 1
		}

		create_stage_from_grs "${stage1_base}.grs" "$template_gin" "$payload_gin" "$stage2_base" "opti conp molmec" "all" "$library_name" \
			"$default_cell_maxcyc" "$default_shell_cuts" "$default_cell_stepmx_opt" "$default_cell_stepmx_rfo" \
			"$default_cell_switch_stepmx" "$default_cell_switch_gnorm"
		run_stage "$stage2_base"
		[ -f "${stage2_base}.grs" ] || {
			echo "Error: stage 2 did not produce ${stage2_base}.grs" >&2
			exit 1
		}

		final_stage="$stage2_base"
		if [ -f "${stage2_base}.gout" ] && extract_energy_value "${stage2_base}.gout" >/dev/null 2>&1; then
			best_incomplete_stage="${stage2_base}.gout"
		fi

		if [ -f "${stage2_base}.gout" ] && (is_failed_output "${stage2_base}.gout" || ! extract_energy_value "${stage2_base}.gout" >/dev/null 2>&1); then
			if [ -f "${stage2_base}.grs" ]; then
				echo "[protocol] Stage 2 requested a shell-free recovery stage."
				create_stage_from_grs "${stage2_base}.grs" "$template_gin" "$payload_gin" "$stage3_base" "opti conp molmec" "core_only" "$library_name" \
					"$default_recovery_maxcyc" "$default_recovery_cuts" "$default_recovery_stepmx_opt" "$default_recovery_stepmx_rfo" \
					"$default_recovery_switch_stepmx" "$default_recovery_switch_gnorm"
				run_stage "$stage3_base"
				if is_complete_output "${stage3_base}.gout"; then
					final_stage="$stage3_base"
				elif is_complete_output "${stage2_base}.gout"; then
					final_stage="$stage2_base"
				elif [ -f "${stage3_base}.gout" ] && extract_energy_value "${stage3_base}.gout" >/dev/null 2>&1; then
					best_incomplete_stage=$(select_better_incomplete_output "$best_incomplete_stage" "${stage3_base}.gout")
					[ -n "$best_incomplete_stage" ] && final_stage="${best_incomplete_stage%.gout}"
				elif [ -n "$best_incomplete_stage" ]; then
					final_stage="${best_incomplete_stage%.gout}"
				fi
			else
				echo "[protocol] Stage 2 failed and no ${stage2_base}.grs was available for recovery." >&2
			fi
		fi

		promote_protocol_result "$canonical_base" "$final_stage"
		if is_complete_output "${final_stage}.gout"; then
			final_complete=1
		fi
		echo "[protocol] Final stage for ${file_name}: ${final_stage}"
		cd "$file_dir"
		if [ "$final_complete" -eq 1 ]; then
			rm -rf "$workspace"
		else
			echo "[protocol] Keeping workspace ${workspace} for inspection because the final result was incomplete." >&2
		fi
	) || return 1

	return 0
}

run_protocol_2() {
	local payload_override=""
	local template_gin=""
	local library_file=""
	local payload_gin=""
	local status=0
	local file=""
	local vasp_files=()

	if ! command -v gulp >/dev/null 2>&1; then
		echo "Error: gulp was not found in PATH." >&2
		return 1
	fi

	if [ $# -eq 2 ] && [[ "$1" == *.vasp ]] && [[ "$2" != *.vasp ]]; then
		payload_override="$2"
		set -- "$1"
	fi

	template_gin=$(protocol_detect_template_gin)
	library_file=$(protocol_detect_library_file)
	payload_gin=$(protocol_detect_template_payload "$payload_override")
	template_gin=$(resolve_existing_path "$template_gin")
	library_file=$(resolve_existing_path "$library_file")
	if [ -n "$payload_gin" ]; then
		payload_gin=$(resolve_existing_path "$payload_gin")
	fi

	if [ $# -gt 0 ]; then
		vasp_files=("$@")
	else
		vasp_files=( *.vasp )
	fi
	[ "${#vasp_files[@]}" -gt 0 ] || {
		echo "Error: no VASP files were provided and none were found in the current directory." >&2
		return 1
	}

	echo "[protocol] Using template GIN: ${template_gin}"
	echo "[protocol] Using library file: ${library_file}"
	if [ -n "$payload_gin" ]; then
		echo "[protocol] Using template payload: ${payload_gin}"
	else
		echo "[protocol] No template payload was found; framework-only coordinates will be used."
	fi

	for file in "${vasp_files[@]}"; do
		echo "[protocol] ============================================================"
		echo "[protocol] Processing ${file}"
		if ! process_protocol_vasp "$file" "$template_gin" "$library_file" "$payload_gin"; then
			echo "[protocol] Failed while processing ${file}" >&2
			status=1
		fi
	done

	return "$status"
}

dispatch_main() {
	local protocol_version="2.0"
	local raw_protocol=""

	while [ $# -gt 0 ]; do
		case "$1" in
			-h|--help)
				print_usage
				exit 0
				;;
			--protocol|--protocole)
				[ $# -ge 2 ] || {
					echo "Error: missing value after $1" >&2
					exit 1
				}
				raw_protocol="$2"
				shift 2
				;;
			--protocol=*|--protocole=*)
				raw_protocol="${1#*=}"
				shift
				;;
			--)
				shift
				break
				;;
			-*)
				echo "Error: unrecognized option $1" >&2
				exit 1
				;;
			*)
				break
				;;
		esac
	done

	protocol_version=$(normalize_protocol_version "${raw_protocol:-2.0}")
	if [ "$protocol_version" = "1.0" ]; then
		if [ $# -lt 1 ] || [ $# -gt 2 ]; then
			echo "Error: protocol 1.0 expects structure.vasp [template.include]" >&2
			exit 1
		fi
		run_protocol_1 "$@"
	else
		run_protocol_2 "$@"
	fi
}

dispatch_main "$@"
