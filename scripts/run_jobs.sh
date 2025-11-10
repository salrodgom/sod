#!/bin/bash
# Orchestrates GULP single-point jobs for the VASP structures in the current directory.
set -euo pipefail
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
idx=0
# Convert every VASP file to a GULP input and launch the calculation pinned to one core.
for vasp in *.vasp; do
	[ -e "$vasp" ] || continue
	bash vasp2gin.sh "$vasp"
	cpu=${CPU_ARRAY[idx]}
	echo "[run_jobs] Ejecutando $vasp.gin en CPU $cpu"
	taskset -c "$cpu" gulp < "$vasp.gin" > "$vasp.gout" &
	idx=$(((idx + 1) % cpu_count))
done
wait