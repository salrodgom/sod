#!/bin/bash
# Orchestrates GULP single-point jobs for the VASP structures in the current directory.
set -euo pipefail

wait_for_slot() {
	while true; do
		for slot in "${!CPU_ARRAY[@]}"; do
			pid="${SLOT_PIDS[$slot]-}"
			if [ -z "${pid:-}" ]; then
				echo "$slot"
				return
			fi
			if ! kill -0 "$pid" 2>/dev/null; then
				wait "$pid" 2>/dev/null || true
				SLOT_PIDS[$slot]=''
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