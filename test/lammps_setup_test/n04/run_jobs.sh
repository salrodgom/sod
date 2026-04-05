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
#ESP Orquesta el flujo completo que convierte ficheros *.vasp en entradas de GULP
#ESP y ejecuta los cálculos respetando los límites por núcleo y globales.

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
	printf '1\n'
}

# Returns the number of GULP processes currently running on the machine.
#ESP Devuelve el número de procesos GULP que se están ejecutando en la máquina.
active_gulp_count() {
	if command -v pgrep >/dev/null 2>&1; then
		pgrep -x gulp 2>/dev/null | wc -l
	else
		ps -eo comm= 2>/dev/null | awk '$1=="gulp" {count++} END {print count+0}'
	fi
}

resolve_lammps_executable() {
	if [ -n "${SOD_LAMMPS_EXECUTABLE:-}" ] && [ -x "${SOD_LAMMPS_EXECUTABLE}" ]; then
		printf '%s\n' "${SOD_LAMMPS_EXECUTABLE}"
		return 0
	fi
	if [ -x "/home/salvador/bin/lmp_fftw" ]; then
		printf '%s\n' "/home/salvador/bin/lmp_fftw"
		return 0
	fi
	if command -v lmp_fftw >/dev/null 2>&1; then
		command -v lmp_fftw
		return 0
	fi
	if command -v lmp >/dev/null 2>&1; then
		command -v lmp
		return 0
	fi
	return 1
}

# Extracts a numeric energy from a GULP output file.
#ESP Extrae una energía numérica de un fichero de salida de GULP.
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

# Checks whether a GULP output explicitly reports a failed optimisation.
#ESP Comprueba si una salida de GULP informa explícitamente de una optimización fallida.
is_failed_output() {
	local gout_file="$1"
	grep -qF 'Conditions for a minimum have not been satisfied.' "$gout_file" || \
	grep -qF 'Too many failed attempts to optimise' "$gout_file" || \
	grep -qF 'Largest core-shell distance exceeds cutoff of cuts' "$gout_file"
}

# Detects the specific core-shell cutoff failure that requires dropping shell coordinates.
#ESP Detecta el fallo específico del corte core-shell que requiere eliminar las coordenadas shell.
has_core_shell_cutoff_error() {
	local gout_file="$1"
	grep -qF 'Largest core-shell distance exceeds cutoff of cuts' "$gout_file"
}

# Decides whether a calculation needs a retry because it failed or produced no usable energy.
#ESP Decide si un cálculo necesita reintento porque falló o no produjo una energía utilizable.
output_requires_retry() {
	local gout_file="$1"

	if [ ! -f "$gout_file" ]; then
		return 0
	fi

	if is_failed_output "$gout_file"; then
		return 0
	fi

	if ! extract_energy_value "$gout_file" >/dev/null; then
		return 0
	fi

	return 1
}

# Extracts the last Gnorm reported by GULP, if present.
#ESP Extrae el último Gnorm informado por GULP, si existe.
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

# Returns success when an output is both numeric and free of explicit failure markers.
#ESP Devuelve éxito cuando una salida es numérica y no contiene marcadores explícitos de fallo.
is_complete_output() {
	local gout_file="$1"
	[ -f "$gout_file" ] || return 1
	is_failed_output "$gout_file" && return 1
	extract_energy_value "$gout_file" >/dev/null
}

# Copies a selected output into the canonical .gout name while preserving the first attempt.
#ESP Copia una salida seleccionada al nombre canónico .gout conservando el primer intento.
promote_output_as_final() {
	local source_gout="$1"
	local canonical_gout="$2"
	local initial_gout="$3"

	if [ "$source_gout" != "$canonical_gout" ]; then
		if [ -f "$canonical_gout" ]; then
			mv -f "$canonical_gout" "$initial_gout"
		fi
		cp "$source_gout" "$canonical_gout"
	fi
}

# Chooses the better of two incomplete outputs using Gnorm first and energy as tiebreaker.
#ESP Elige la mejor de dos salidas incompletas usando primero Gnorm y la energía como desempate.
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
				if (ga < gb - gtol) {
					print "current"
				} else if (gb < ga - gtol) {
					print "candidate"
				} else if (ea < eb - etol) {
					print "current"
				} else if (eb < ea - etol) {
					print "candidate"
				} else {
					print "current"
				}
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
				if (ea <= eb + tol) {
					print "current"
				} else {
					print "candidate"
				}
			}
		')
	fi

	if [ "$decision" = "candidate" ]; then
		printf '%s\n' "$candidate_gout"
	else
		printf '%s\n' "$current_gout"
	fi
}

# Rebuilds a restart input from the previous .grs using a staged retry strategy.
#ESP Reconstruye una entrada de reinicio a partir del .grs usando una estrategia de reintento por etapas.
build_retry_input() {
	local vasp_file="$1"
	local restart_base="$2"
	local coord_mode="$3"
	local retry_style="$4"
	local source_grs_file="$5"
	local gin_file="${vasp_file}.gin"
	local grs_file="$source_grs_file"
	local restart_gin="${restart_base}.gin"
	local keywords_line=""
	local cell_line=""
	local fractional_line=""

	if [ ! -f "$gin_file" ]; then
		echo "[run_jobs] Missing original input ${gin_file} for restart." >&2
		return 1
	fi

	if [ ! -f "$grs_file" ]; then
		echo "[run_jobs] Missing restart file ${grs_file} for ${vasp_file}." >&2
		return 1
	fi

	keywords_line=$(awk 'NF { print; exit }' "$gin_file")
	cell_line=$(awk 'tolower($1)=="cell" { getline; if (NF) { print; exit } }' "$grs_file")
	fractional_line=$(awk 'tolower($1)=="fractional" { print; exit }' "$grs_file")

	if [ -z "$keywords_line" ] || [ -z "$cell_line" ] || [ -z "$fractional_line" ]; then
		echo "[run_jobs] Could not reconstruct restart header for ${vasp_file}." >&2
		return 1
	fi

	if ! awk '/^[[:space:]]*species([[:space:]]|$)/ { found=1; exit } END { exit !found }' "$gin_file"; then
		echo "[run_jobs] Missing species block in ${gin_file}; restart cannot be generated." >&2
		return 1
	fi

	{
		printf '%s\n' "$keywords_line"
		printf 'name %s\n' "$restart_base"
		printf 'cell\n'
		printf '%s\n' "$cell_line"
		printf '%s\n' "$fractional_line"
		awk -v coord_mode="$coord_mode" '
			tolower($1)=="fractional" { in_coords=1; next }
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
		' "$grs_file"
		printf '\n'
		awk '
			tolower($1)=="frac" || tolower($1)=="fractional" {
				after_coords = 1
				next
			}
			!after_coords {
				next
			}
			tolower($1)=="species" {
				exit
			}
			NF >= 2 && (tolower($2)=="core" || tolower($2)=="shel") {
				next
			}
			{
				print
			}
		' "$gin_file"
		printf '\n'
		awk -v restart_base="$restart_base" -v retry_style="$retry_style" '
			BEGIN {
				in_tail = 0
			}
			/^[[:space:]]*species([[:space:]]|$)/ {
				in_tail = 1
			}
			!in_tail {
				next
			}
			/^[[:space:]]*library([[:space:]]|$)/ {
				line = $0
				if (line !~ /(^|[[:space:]])nodump([[:space:]]|$)/) {
					line = line " nodump"
				}
				print line
				next
			}
			/^[[:space:]]*stepmx([[:space:]]|$)/ {
				if (retry_style != "conservative") {
					next
				}
			}
			/^[[:space:]]*switch_stepmx([[:space:]]|$)/ {
				if (retry_style != "conservative") {
					next
				}
			}
			/^[[:space:]]*cuts([[:space:]]|$)/ {
				if (retry_style != "conservative") {
					next
				}
			}
			/^[[:space:]]*dump([[:space:]]|$)/ { next }
			/^[[:space:]]*output([[:space:]]|$)/ { next }
			{
				print
			}
			END {
				if (retry_style == "aggressive") {
					print "stepmx opt 0.01"
					print "stepmx rfo 0.01"
					print "switch_stepmx 0.005 gnorm 0.02"
					print "cuts 1.2"
				} else if (retry_style == "very_aggressive") {
					print "stepmx opt 0.005"
					print "stepmx rfo 0.005"
					print "switch_stepmx 0.002 gnorm 0.01"
					print "cuts 1.3"
				}
				print "dump every 1 " restart_base ".grs"
				print "output cif " restart_base ".cif"
			}
		' "$gin_file"
	} > "$restart_gin"
}

# Runs one restart stage and classifies its result for the caller.
#ESP Ejecuta una etapa de reinicio y clasifica su resultado para quien la llama.
run_restart_stage() {
	local vasp_file="$1"
	local cpu="$2"
	local restart_base="$3"
	local coord_mode="$4"
	local retry_style="$5"
	local source_grs_file="$6"
	local restart_gin="${restart_base}.gin"
	local restart_gout="${restart_base}.gout"
	local stage_label=""

	STAGE_STATUS="build_failed"
	STAGE_GOUT="$restart_gout"

	if [ "$retry_style" = "conservative" ]; then
		stage_label="conservative restart"
	elif [ "$retry_style" = "aggressive" ]; then
		stage_label="aggressive restart"
	else
		stage_label="third restart"
	fi

	if ! build_retry_input "$vasp_file" "$restart_base" "$coord_mode" "$retry_style" "$source_grs_file"; then
		echo "[run_jobs] Could not build ${stage_label} input ${restart_gin} for ${vasp_file}." >&2
		return 0
	fi

	echo "[run_jobs] Running ${stage_label} for ${vasp_file} on CPU ${cpu} with ${restart_gin}."
	if ! taskset -c "$cpu" gulp < "$restart_gin" > "$restart_gout"; then
		echo "[run_jobs] ${stage_label^} execution failed for ${vasp_file}." >&2
	fi

	if is_complete_output "$restart_gout"; then
		STAGE_STATUS="complete"
	elif extract_energy_value "$restart_gout" >/dev/null; then
		STAGE_STATUS="incomplete_numeric"
	else
		STAGE_STATUS="invalid"
	fi
}

# Runs the three-stage restart strategy and keeps the best available result.
#ESP Ejecuta la estrategia de reinicio en tres etapas y conserva el mejor resultado disponible.
retry_from_restart() {
	local vasp_file="$1"
	local cpu="$2"
	local gout_file="${vasp_file}.gout"
	local initial_gout="${vasp_file}.initial.gout"
	local best_incomplete_gout=""
	local fallback_gout=""
	local conservative_base="${vasp_file}.restart"
	local aggressive_base="${vasp_file}.restart2"
	local third_base="${vasp_file}.restart3"
	local original_grs="${vasp_file}.grs"
	local aggressive_source_grs=""
	local third_source_grs=""
	local selected_label=""
	local skip_conservative=0

	if [ -f "$gout_file" ] && extract_energy_value "$gout_file" >/dev/null; then
		best_incomplete_gout="$gout_file"
	fi

	if [ -f "$gout_file" ] && has_core_shell_cutoff_error "$gout_file"; then
		skip_conservative=1
		echo "[run_jobs] ${gout_file} exceeded the core-shell cutoff; skipping the conservative restart and switching directly to a shell-free restart." >&2
	fi

	if [ "$skip_conservative" -eq 0 ]; then
		run_restart_stage "$vasp_file" "$cpu" "$conservative_base" "all" "conservative" "$original_grs"
		case "$STAGE_STATUS" in
			complete)
				echo "[run_jobs] Conservative restart converged for ${vasp_file}; promoting ${STAGE_GOUT}."
				promote_output_as_final "$STAGE_GOUT" "$gout_file" "$initial_gout"
				return 0
				;;
			incomplete_numeric)
				echo "[run_jobs] Conservative restart remained incomplete for ${vasp_file}; storing it as a candidate result." >&2
				best_incomplete_gout=$(select_better_incomplete_output "$best_incomplete_gout" "$STAGE_GOUT")
				fallback_gout="$STAGE_GOUT"
				;;
			invalid)
				echo "[run_jobs] Conservative restart produced no numeric energy for ${vasp_file}." >&2
				fallback_gout="$STAGE_GOUT"
				;;
		esac

		aggressive_source_grs="${conservative_base}.grs"
		if [ ! -f "$aggressive_source_grs" ]; then
			aggressive_source_grs="$original_grs"
		fi
	else
		aggressive_source_grs="$original_grs"
	fi

	run_restart_stage "$vasp_file" "$cpu" "$aggressive_base" "core_only" "aggressive" "$aggressive_source_grs"
	case "$STAGE_STATUS" in
		complete)
			echo "[run_jobs] Aggressive restart converged for ${vasp_file}; promoting ${STAGE_GOUT}."
			promote_output_as_final "$STAGE_GOUT" "$gout_file" "$initial_gout"
			return 0
			;;
		incomplete_numeric)
			echo "[run_jobs] Aggressive restart remained incomplete for ${vasp_file}; comparing it with the other candidates." >&2
			best_incomplete_gout=$(select_better_incomplete_output "$best_incomplete_gout" "$STAGE_GOUT")
			fallback_gout="$STAGE_GOUT"
			;;
		invalid)
			echo "[run_jobs] Aggressive restart produced no numeric energy for ${vasp_file}." >&2
			fallback_gout="$STAGE_GOUT"
			;;
	esac

	third_source_grs="${aggressive_base}.grs"
	if [ ! -f "$third_source_grs" ]; then
		third_source_grs="$aggressive_source_grs"
	fi

	run_restart_stage "$vasp_file" "$cpu" "$third_base" "core_only" "very_aggressive" "$third_source_grs"
	case "$STAGE_STATUS" in
		complete)
			echo "[run_jobs] Third restart converged for ${vasp_file}; promoting ${STAGE_GOUT}."
			promote_output_as_final "$STAGE_GOUT" "$gout_file" "$initial_gout"
			return 0
			;;
		incomplete_numeric)
			echo "[run_jobs] Third restart remained incomplete for ${vasp_file}; comparing it with the other candidates." >&2
			best_incomplete_gout=$(select_better_incomplete_output "$best_incomplete_gout" "$STAGE_GOUT")
			fallback_gout="$STAGE_GOUT"
			;;
		invalid)
			echo "[run_jobs] Third restart produced no numeric energy for ${vasp_file}." >&2
			fallback_gout="$STAGE_GOUT"
			;;
	esac

	if [ -n "$best_incomplete_gout" ] && [ -f "$best_incomplete_gout" ]; then
		if [ "$best_incomplete_gout" = "$gout_file" ]; then
			echo "[run_jobs] Keeping the initial incomplete optimisation for ${vasp_file}; it remains the best candidate." >&2
		else
			if [ "$best_incomplete_gout" = "${conservative_base}.gout" ]; then
				selected_label="conservative restart"
			elif [ "$best_incomplete_gout" = "${aggressive_base}.gout" ]; then
				selected_label="aggressive restart"
			else
				selected_label="third restart"
			fi
			echo "[run_jobs] Keeping the ${selected_label} result for ${vasp_file}; it is the best incomplete candidate." >&2
			promote_output_as_final "$best_incomplete_gout" "$gout_file" "$initial_gout"
		fi
		return 0
	fi

	if [ -n "$fallback_gout" ] && [ -f "$fallback_gout" ] && [ "$fallback_gout" != "$gout_file" ]; then
		echo "[run_jobs] No numeric energy recovered for ${vasp_file}; promoting ${fallback_gout} so extraction can emit N/A." >&2
		promote_output_as_final "$fallback_gout" "$gout_file" "$initial_gout"
	fi
	return 0
}

# Runs one VASP->GULP job, performing staged restarts when the first pass fails.
#ESP Ejecuta un trabajo VASP->GULP y realiza reinicios por etapas cuando falla la primera pasada.
run_job_pipeline() {
	local vasp_file="$1"
	local cpu="$2"
	local gin_file="${vasp_file}.gin"
	local gout_file="${vasp_file}.gout"

	echo "[run_jobs] Running ${gin_file} on CPU ${cpu}"
	if taskset -c "$cpu" gulp < "$gin_file" > "$gout_file"; then
		if ! output_requires_retry "$gout_file"; then
			return 0
		fi
		echo "[run_jobs] ${gout_file} requested a retry." >&2
	else
		echo "[run_jobs] Initial execution failed for ${vasp_file}; attempting restart." >&2
	fi

	retry_from_restart "$vasp_file" "$cpu"
}

# Runs one VASP job through protocol 2.0 while keeping the slot-level scheduler.
#ESP Ejecuta un trabajo VASP mediante el protocolo 2.0 manteniendo el planificador por slots.
run_job_protocol_2() {
	local vasp_file="$1"
	local cpu="$2"
	local slot="$3"
	local slot_label="slot_$((slot + 1))"

	if [ ! -f "./vasp2gin.sh" ]; then
		echo "[run_jobs] Missing canonical converter ./vasp2gin.sh." >&2
		return 1
	fi

	echo "[run_jobs] Running protocol 2.0 for ${vasp_file} on CPU ${cpu} using workspace ${slot_label}"
	if [ -n "${TEMPLATE_GIN_FILE:-}" ]; then
		SOD_GULP_SLOT_INDEX="$((slot + 1))" SOD_GULP_SLOT_LABEL="$slot_label" \
			taskset -c "$cpu" bash ./vasp2gin.sh --protocol 2.0 "$vasp_file" "$TEMPLATE_GIN_FILE"
	else
		SOD_GULP_SLOT_INDEX="$((slot + 1))" SOD_GULP_SLOT_LABEL="$slot_label" \
			taskset -c "$cpu" bash ./vasp2gin.sh --protocol 2.0 "$vasp_file"
	fi
}

run_job_lammps() {
	local input_file="$1"
	local cpu="$2"
	local output_file="${input_file%.in.lammps}.out.lammps"
	local log_file="${input_file%.in.lammps}.log.lammps"

	echo "[run_jobs] Running ${input_file} on CPU ${cpu} with ${LAMMPS_EXECUTABLE}"
	taskset -c "$cpu" "$LAMMPS_EXECUTABLE" -log "$log_file" -in "$input_file" > "$output_file" 2>&1
}

# Waits until the global concurrency limit leaves room for one more job.
#ESP Espera hasta que el límite global de concurrencia deje hueco para un trabajo más.
wait_for_global_capacity() {
	if [ "$GLOBAL_GULP_LIMIT" -le 0 ]; then
		return
	fi
	while [ "$(active_gulp_count)" -ge "$GLOBAL_GULP_LIMIT" ]; do
		sleep 0.2
		# Loop until at least one global slot is free
		#ESP Bucle hasta que se libere al menos un hueco global.
	done
}

# Reaps a finished slot and records whether the background pipeline failed.
#ESP Recoge un slot terminado y registra si el flujo en segundo plano falló.
reap_slot() {
	local slot="$1"
	local pid="${SLOT_PIDS[$slot]-}"

	if [ -z "${pid:-}" ]; then
		return
	fi

	if ! wait "$pid"; then
		HAD_FAILURE=1
	fi
	SLOT_PIDS[$slot]=''
}

# Selects the next available per-core slot, cleaning up finished jobs on the way.
#ESP Selecciona el siguiente slot libre por núcleo, limpiando los trabajos terminados por el camino.
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
				reap_slot "$slot"
				wait_for_global_capacity
				echo "$slot"
				return
			fi
		done
		sleep 0.2
	done
}

# main:
FILER_VALUE="$(read_filer_value)"

# Discover the list of CPU cores that will be used to pin each GULP process.
#ESP Descubre la lista de núcleos de CPU que se usarán para fijar cada proceso GULP.
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
	echo "[run_jobs] Empty CPU list." >&2
	exit 1
fi

# Determine the global cap for concurrent GULP runs (defaults to local slots).
#ESP Determina el límite global de ejecuciones GULP concurrentes (por defecto, los slots locales).
if [ -n "${SOD_GULP_GLOBAL_LIMIT:-}" ]; then
	GLOBAL_GULP_LIMIT=$SOD_GULP_GLOBAL_LIMIT
else
	GLOBAL_GULP_LIMIT=$cpu_count
fi
if ! [[ "$GLOBAL_GULP_LIMIT" =~ ^-?[0-9]+$ ]]; then
	echo "[run_jobs] Non-numeric value for SOD_GULP_GLOBAL_LIMIT: $GLOBAL_GULP_LIMIT" >&2
	GLOBAL_GULP_LIMIT=$cpu_count
fi
if [ "$GLOBAL_GULP_LIMIT" -eq 0 ]; then
	GLOBAL_GULP_LIMIT=-1
fi

if [ "$FILER_VALUE" = "2" ]; then
	GLOBAL_GULP_LIMIT=-1
	if ! LAMMPS_EXECUTABLE="$(resolve_lammps_executable)"; then
		echo "[run_jobs] LAMMPS executable not found. Set SOD_LAMMPS_EXECUTABLE or install lmp_fftw/lmp." >&2
		exit 1
	fi
fi

declare -a SLOT_PIDS
for idx in "${!CPU_ARRAY[@]}"; do
	SLOT_PIDS[$idx]=''
done
HAD_FAILURE=0

TEMPLATE_GIN_FILE="${SOD_TEMPLATE_GIN_FILE:-}"
if [ -z "$TEMPLATE_GIN_FILE" ]; then
	if [ -f "template_payload.gin" ]; then
		TEMPLATE_GIN_FILE="./template_payload.gin"
	elif [ -f "default_template.gin" ]; then
		TEMPLATE_GIN_FILE="./default_template.gin"
	elif [ -f "default_template.include" ]; then
		TEMPLATE_GIN_FILE="./default_template.include"
	fi
fi
if [ "$FILER_VALUE" = "2" ]; then
	for input_file in *.in.lammps; do
		[ -e "$input_file" ] || continue
		slot=$(wait_for_slot)
		cpu=${CPU_ARRAY[$slot]}
		run_job_lammps "$input_file" "$cpu" &
		SLOT_PIDS[$slot]=$!
	done
else
	PROTOCOL_VERSION="${SOD_GULP_PROTOCOL_VERSION:-2.0}"
	case "${PROTOCOL_VERSION}" in
		1|1.0) PROTOCOL_VERSION="1.0" ;;
		2|2.0) PROTOCOL_VERSION="2.0" ;;
		*)
			echo "[run_jobs] Unknown protocol version '${PROTOCOL_VERSION}'. Supported values: 1.0 or 2.0." >&2
			exit 1
			;;
	esac
	# Convert every VASP file to a GULP input and launch the calculation pinned to one core.
	#ESP Convierte cada fichero VASP en una entrada de GULP y lanza el cálculo fijado a un núcleo.
	for vasp in *.vasp; do
		[ -e "$vasp" ] || continue
		if [ "$PROTOCOL_VERSION" = "2.0" ]; then
			slot=$(wait_for_slot)
			cpu=${CPU_ARRAY[$slot]}
			run_job_protocol_2 "$vasp" "$cpu" "$slot" &
		else
			if [ -n "$TEMPLATE_GIN_FILE" ]; then
				bash vasp2gin.sh --protocol 1.0 "$vasp" "$TEMPLATE_GIN_FILE"
			else
				bash vasp2gin.sh --protocol 1.0 "$vasp"
			fi
			slot=$(wait_for_slot)
			cpu=${CPU_ARRAY[$slot]}
			run_job_pipeline "$vasp" "$cpu" &
		fi
		SLOT_PIDS[$slot]=$!
	done
fi
for slot in "${!SLOT_PIDS[@]}"; do
	reap_slot "$slot"
done

exit "$HAD_FAILURE"
