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
trap 'rm -f "$temp_energies"' EXIT

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
done

mv "$temp_energies" ENERGIES
trap - EXIT
