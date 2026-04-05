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

set -euo pipefail
shopt -s nullglob

read_filer_value() {
	if [ -f "filer" ]; then
		awk 'NF { print $1; exit }' filer
		return 0
	fi
	if compgen -G '*.in.lammps' >/dev/null; then
		printf '2\n'
		return 0
	fi
	if [ -f "template_ase.py" ] && compgen -G '*.vasp' >/dev/null; then
		printf '14\n'
		return 0
	fi
	printf '1\n'
}

extract_energy_value() {
	local gout_file="$1"
	local energy=""

	energy=$(awk '
		/Final energy[[:space:]]*=/ && $5 == "eV" { value = $4 }
		END {
			if (value != "") {
				print value
			}
		}
	' "$gout_file")

	if [ -z "$energy" ]; then
		energy=$(awk '
			/Total lattice energy[[:space:]]*=/ && $6 == "eV" { value = $5 }
			END {
				if (value != "") {
					print value
				}
			}
		' "$gout_file")
	fi

	if [[ "$energy" =~ ^[-+]?[0-9]+([.][0-9]+)?([Ee][-+]?[0-9]+)?$ ]]; then
		printf '%s\n' "$energy"
		return 0
	fi

	return 1
}

extract_lammps_energy_value() {
	local output_file="$1"
	local energy=""

	energy=$(awk '
		/SOD_FINAL_ENERGY/ {
			for (i = 1; i <= NF; i++) {
				if ($i == "SOD_FINAL_ENERGY" && (i + 1) <= NF) {
					value = $(i + 1)
				}
			}
		}
		END {
			if (value != "") {
				print value
			}
		}
	' "$output_file")

	if [[ "$energy" =~ ^[-+]?[0-9]+([.][0-9]+)?([Ee][-+]?[0-9]+)?$ ]]; then
		printf '%s\n' "$energy"
		return 0
	fi

	return 1
}

extract_ase_energy_value() {
	local output_file="$1"
	local energy=""

	energy=$(awk '
		/SOD_FINAL_ENERGY/ {
			for (i = 1; i <= NF; i++) {
				if ($i == "SOD_FINAL_ENERGY" && (i + 1) <= NF) {
					value = $(i + 1)
				}
			}
		}
		END {
			if (value != "") {
				print value
			}
		}
	' "$output_file")

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

failure_reason() {
	local gout_file="$1"
	local reasons=()
	local joined_reason=""
	local idx

	if grep -qF 'Too many failed attempts to optimise' "$gout_file"; then
		reasons+=("Too many failed attempts to optimise")
	fi
	if grep -qF 'Conditions for a minimum have not been satisfied.' "$gout_file"; then
		reasons+=("Conditions for a minimum have not been satisfied")
	fi
	if grep -qF 'Largest core-shell distance exceeds cutoff of cuts' "$gout_file"; then
		reasons+=("Largest core-shell distance exceeds cutoff of cuts")
	fi

	if [ "${#reasons[@]}" -eq 0 ]; then
		return 1
	fi

	joined_reason="${reasons[0]}"
	for ((idx = 1; idx < ${#reasons[@]}; idx++)); do
		joined_reason="${joined_reason} | ${reasons[idx]}"
	done
	printf '%s\n' "$joined_reason"
}

last_gnorm_value() {
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
		END {
			if (value != "") {
				print value
			}
		}
	' "$gout_file")

	if [[ "$gnorm" =~ ^[-+]?[0-9]+([.][0-9]+)?([Ee][-+]?[0-9]+)?$ ]]; then
		printf '%s\n' "$gnorm"
		return 0
	fi

	return 1
}

temp_energies=$(mktemp ENERGIES.tmp.XXXXXX)
temp_relaxed=$(mktemp RELAXED_STRUCTURES.tmp.XXXXXX)
trap 'rm -f "$temp_energies" "$temp_relaxed"' EXIT

FILER_VALUE="$(read_filer_value)"

if [ "$FILER_VALUE" = "2" ]; then
	for input_file in *.in.lammps; do
		[ -e "$input_file" ] || continue
		base_name="${input_file%.in.lammps}"
		output_file="${base_name}.out.lammps"
		relaxed_file="${base_name}.min.data"

		if [ ! -f "$output_file" ]; then
			echo "[extract] Missing output file ${output_file}; writing N/A." >&2
			echo "N/A # $output_file # MISSING_OUTPUT" >> "$temp_energies"
			echo "N/A # $output_file # MISSING_OUTPUT" >> "$temp_relaxed"
			continue
		fi

		if extract_lammps_energy_value "$output_file" >/dev/null; then
			energy_value="$(extract_lammps_energy_value "$output_file")"
			echo "$energy_value # $output_file" >> "$temp_energies"
		else
			echo "[extract] No numeric LAMMPS energy could be extracted from ${output_file}; writing N/A." >&2
			echo "N/A # $output_file # NO_NUMERIC_ENERGY" >> "$temp_energies"
		fi

		if [ -f "$relaxed_file" ]; then
			echo "$relaxed_file # $output_file" >> "$temp_relaxed"
		else
			echo "N/A # $output_file # MISSING_RELAXED_STRUCTURE" >> "$temp_relaxed"
		fi
	done

	mv "$temp_energies" ENERGIES
	mv "$temp_relaxed" RELAXED_STRUCTURES
	trap - EXIT
	exit 0
fi

if [ "$FILER_VALUE" = "14" ]; then
	for vasp_file in c[0-9][0-9][0-9][0-9][0-9].vasp; do
		[ -e "$vasp_file" ] || continue
		base_name="${vasp_file%.vasp}"
		output_file="${base_name}.out.ase"
		relaxed_cif="${base_name}.relaxed.cif"
		relaxed_vasp="${base_name}.relaxed.vasp"

		if [ ! -f "$output_file" ]; then
			echo "[extract] Missing output file ${output_file}; writing N/A." >&2
			echo "N/A # $output_file # MISSING_OUTPUT" >> "$temp_energies"
			echo "N/A # $output_file # MISSING_OUTPUT" >> "$temp_relaxed"
			continue
		fi

		if extract_ase_energy_value "$output_file" >/dev/null; then
			energy_value="$(extract_ase_energy_value "$output_file")"
			echo "$energy_value # $output_file" >> "$temp_energies"
		else
			echo "[extract] No numeric ASE energy could be extracted from ${output_file}; writing N/A." >&2
			echo "N/A # $output_file # NO_NUMERIC_ENERGY" >> "$temp_energies"
		fi

		if [ -f "$relaxed_cif" ]; then
			echo "$relaxed_cif # $output_file" >> "$temp_relaxed"
		elif [ -f "$relaxed_vasp" ]; then
			echo "$relaxed_vasp # $output_file" >> "$temp_relaxed"
		else
			echo "N/A # $output_file # MISSING_RELAXED_STRUCTURE" >> "$temp_relaxed"
		fi
	done

	mv "$temp_energies" ENERGIES
	mv "$temp_relaxed" RELAXED_STRUCTURES
	trap - EXIT
	exit 0
fi

for vasp_file in *.vasp; do
	[ -e "$vasp_file" ] || continue
	gout_file="${vasp_file}.gout"
	gnorm_note=""

	if [ ! -f "$gout_file" ]; then
		echo "[extract] Missing output file ${gout_file}; writing N/A." >&2
		echo "N/A # $gout_file # MISSING_OUTPUT" >> "$temp_energies"
		continue
	fi

	if gnorm_value=$(last_gnorm_value "$gout_file"); then
		gnorm_note=" # LAST_GNORM $gnorm_value"
	fi

	if is_failed_output "$gout_file"; then
		failure_note=$(failure_reason "$gout_file" || printf '%s' 'Unknown incomplete optimisation')
		if ! energy_value=$(extract_energy_value "$gout_file"); then
			echo "[extract] Failed optimisation marker found in ${gout_file}, and no numeric energy could be recovered; writing N/A." >&2
			echo "N/A # $gout_file${gnorm_note} # INCOMPLETE_OPTIMISATION # $failure_note" >> "$temp_energies"
			continue
		fi
		echo "[extract] Incomplete optimisation detected in ${gout_file}; keeping the last numeric energy." >&2
		echo "$energy_value # $gout_file${gnorm_note} # INCOMPLETE_OPTIMISATION # $failure_note" >> "$temp_energies"
		continue
	fi

	if ! energy_value=$(extract_energy_value "$gout_file"); then
		echo "[extract] No numeric energy could be extracted from ${gout_file}; writing N/A." >&2
		echo "N/A # $gout_file${gnorm_note} # NO_NUMERIC_ENERGY" >> "$temp_energies"
		continue
	fi

	echo "$energy_value # $gout_file${gnorm_note}" >> "$temp_energies"

	if [ -f "${vasp_file}.cif" ]; then
		echo "${vasp_file}.cif # $gout_file" >> "$temp_relaxed"
	else
		echo "N/A # $gout_file # MISSING_RELAXED_STRUCTURE" >> "$temp_relaxed"
	fi
done

mv "$temp_energies" ENERGIES
mv "$temp_relaxed" RELAXED_STRUCTURES
trap - EXIT
