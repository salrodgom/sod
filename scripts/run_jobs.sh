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

# Orchestrates the end-to-end pipeline that converts *.vasp files to GULP inputs
# and executes the calculations while respecting per-core and global limits.

set -euo pipefail

# Returns the number of GULP processes currently running on the machine.
active_gulp_count() {
	if command -v pgrep >/dev/null 2>&1; then
		pgrep -x gulp 2>/dev/null | wc -l
	else
		ps -eo comm= 2>/dev/null | awk '$1=="gulp" {count++} END {print count+0}'
	fi
}

# Waits until the global concurrency limit leaves room for one more job.
wait_for_global_capacity() {
	if [ "$GLOBAL_GULP_LIMIT" -le 0 ]; then
		return
	fi
	while [ "$(active_gulp_count)" -ge "$GLOBAL_GULP_LIMIT" ]; do
		sleep 0.2
		# Loop until at least one global slot is free
	done
}

# Selects the next available per-core slot, cleaning up finished jobs on the way.
wait_for_slot() {
	while true; do
		for slot in "${!CPU_ARRAY[@]}"; do
			pid="${SLOT_PIDS[$slot]-}"
			if [ -z "${pid:-}" ]; then
				wait_for_global_capacity
				echo "$slot"
				return
			fi
			if ! kill -0 "$pid" 2>/dev/null; then
				wait "$pid" 2>/dev/null || true
				SLOT_PIDS[$slot]=''
				wait_for_global_capacity
				echo "$slot"
				return
			fi
		done
		sleep 0.2
	done
}

# main:
# Discover the list of CPU cores that will be used to pin each GULP process.
if [ -z "${SOD_GULP_CPUS:-}" ]; then
	total_cores=$(nproc)
	if [ "$total_cores" -lt 1 ]; then
		echo "[run_jobs] No CPU cores detected." >&2
		exit 1
	fi
	cpu_list=$(seq -s, 0 $((total_cores - 1)))
else
	cpu_list=$SOD_GULP_CPUS
fi
IFS=',' read -r -a CPU_ARRAY <<< "$cpu_list"
cpu_count=${#CPU_ARRAY[@]}
if [ "$cpu_count" -eq 0 ]; then
	echo "[run_jobs] Lista de CPUs vacía." >&2
	exit 1
fi

# Determine the global cap for concurrent GULP runs (defaults to local slots).
if [ -n "${SOD_GULP_GLOBAL_LIMIT:-}" ]; then
	GLOBAL_GULP_LIMIT=$SOD_GULP_GLOBAL_LIMIT
else
	GLOBAL_GULP_LIMIT=$cpu_count
fi
if ! [[ "$GLOBAL_GULP_LIMIT" =~ ^-?[0-9]+$ ]]; then
	echo "[run_jobs] Valor no numérico para SOD_GULP_GLOBAL_LIMIT: $GLOBAL_GULP_LIMIT" >&2
	GLOBAL_GULP_LIMIT=$cpu_count
fi
if [ "$GLOBAL_GULP_LIMIT" -eq 0 ]; then
	GLOBAL_GULP_LIMIT=-1
fi

declare -a SLOT_PIDS
for idx in "${!CPU_ARRAY[@]}"; do
	SLOT_PIDS[$idx]=''
done
# Convert every VASP file to a GULP input and launch the calculation pinned to one core.
for vasp in *.vasp; do
	[ -e "$vasp" ] || continue
	bash vasp2gin.sh "$vasp"
	slot=$(wait_for_slot)
	cpu=${CPU_ARRAY[$slot]}
	echo "[run_jobs] Ejecutando $vasp.gin en CPU $cpu"
	taskset -c "$cpu" gulp < "$vasp.gin" > "$vasp.gout" &
	SLOT_PIDS[$slot]=$!
done
for pid in "${SLOT_PIDS[@]}"; do
	if [ -n "${pid:-}" ]; then
		wait "$pid" || true
	fi
done