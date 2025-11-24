#!/bin/bash
# Orchestrates GULP single-point jobs for the VASP structures in the current directory.
set -euo pipefail

active_gulp_count() {
	if command -v pgrep >/dev/null 2>&1; then
		pgrep -x gulp 2>/dev/null | wc -l
	else
		ps -eo comm= 2>/dev/null | awk '$1=="gulp" {count++} END {print count+0}'
	fi
}

wait_for_global_capacity() {
	if [ "$GLOBAL_GULP_LIMIT" -le 0 ]; then
		return
	fi
	while [ "$(active_gulp_count)" -ge "$GLOBAL_GULP_LIMIT" ]; do
		sleep 0.2
		# Loop until at least one global slot is free
	done
}

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