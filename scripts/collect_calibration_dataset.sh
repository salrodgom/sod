#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_SOD_BIN="${REPO_ROOT}/bin/sod"
DEFAULT_SCRIPTS_DIR="${REPO_ROOT}/scripts"

SOD_BIN="${DEFAULT_SOD_BIN}"
SCRIPTS_DIR="${DEFAULT_SCRIPTS_DIR}"
LEVEL_SPEC=""
TOTAL_SITES=""
TEMPERATURE="1000"
UNIFORM_SAMPLES="300"
EXACT_MAX_COMBOS="150000"
SEED="-1"
PROTOCOL_VERSION="${SOD_GULP_PROTOCOL_VERSION:-2.0}"
TEMPLATE_GIN=""
FORCE=0
DRY_RUN=0
RESET_SUMMARY=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Collects calibration data over several N levels, choosing a cheap mode per level:
  - exact      for levels whose combinatorics are still moderate
  - mc uniform for larger levels, to cover more unique states per unit cost

Run this from the calculation directory that contains INSOD and SGO.

Options:
  --levels <spec>            Levels to process, e.g. 5:10 or 5,6,8,10
                             Default: auto -> 5:min(npos-5,10)
  --total-sites <npos>       Override inferred number of substitutable sites
  -T, --temperature <K>      Temperature for MC uniform runs [default: 1000]
  -C, --samples <N>          Uniform MC samples per sampled level [default: 300]
  --exact-max-combos <N>     Use exact when C(npos,N) <= this threshold [default: 150000]
  -s, --seed <value>         Seed for MC runs [default: -1]
  --sod-bin <path>           Path to sod [default: ${DEFAULT_SOD_BIN}]
  --scripts-dir <path>       Path to scripts/ for GULP helpers [default: ${DEFAULT_SCRIPTS_DIR}]
  --template-gin <file>      Pass through to exact/mc
  --no-template-gin          Skip template fragments explicitly (default behavior)
  --protocol <ver>           Pass through to exact/mc [default: 2.0]
  --force                    Recompute levels even if a per-level report already exists
  --reset-summary            Remove cumulative CSV summaries before starting
  --dry-run                  Print commands without executing them
  -h, --help                 Show this help

Outputs collected in the calculation root:
  CALIBRATION_COEFFICIENTS_SUMMARY.csv
  CALIBRATION_BLEND_OVERRIDES_SUMMARY.csv

Per-level detailed reports:
  xNN/CALIBRATION_COEFFICIENTS.txt
  mcNN/CALIBRATION_COEFFICIENTS.txt

Example:
  $(basename "$0") --levels 5:10 --temperature 1000 --no-template-gin
EOF
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

info() {
    printf '%s\n' "$*"
}

run_cmd() {
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        printf '[dry-run] '
        printf '%q ' "$@"
        printf '\n'
    else
        "$@"
    fi
}

looks_like_scripts_dir() {
    local dir="$1"
    [[ -d "${dir}" ]] || return 1
    [[ -f "${dir}/run_jobs.sh" ]] || return 1
    [[ -f "${dir}/extract.sh" ]] || return 1
    return 0
}

infer_total_sites() {
    local output
    local npos

    output="$("${SOD_BIN}" exact --just-outsod -N 0 --no-template-gin 2>&1 || true)"
    npos="$(printf '%s\n' "${output}" | awk '/Substitutable sites \(npos\):/ {print $NF; exit}')"
    [[ -n "${npos}" ]] || return 1
    printf '%s\n' "${npos}"
}

binomial_int() {
    local n="$1"
    local k="$2"
    local i
    local result=1
    local kk

    if (( k < 0 || k > n )); then
        printf '0\n'
        return
    fi
    if (( k == 0 || k == n )); then
        printf '1\n'
        return
    fi
    kk="$k"
    if (( kk > n - kk )); then
        kk=$((n - kk))
    fi
    for ((i = 1; i <= kk; i++)); do
        result=$(( result * (n - kk + i) / i ))
    done
    printf '%s\n' "${result}"
}

expand_level_spec() {
    local spec="$1"
    local token start end value
    local -n out_ref="$2"

    out_ref=()
    IFS=',' read -r -a tokens <<< "${spec}"
    for token in "${tokens[@]}"; do
        token="${token// /}"
        [[ -n "${token}" ]] || continue
        if [[ "${token}" == *:* ]]; then
            start="${token%%:*}"
            end="${token##*:}"
            [[ "${start}" =~ ^-?[0-9]+$ && "${end}" =~ ^-?[0-9]+$ ]] || die "invalid level range: ${token}"
            if (( start <= end )); then
                for ((value = start; value <= end; value++)); do
                    out_ref+=("${value}")
                done
            else
                for ((value = start; value >= end; value--)); do
                    out_ref+=("${value}")
                done
            fi
        else
            [[ "${token}" =~ ^-?[0-9]+$ ]] || die "invalid level value: ${token}"
            out_ref+=("${token}")
        fi
    done
}

format_level_dir() {
    local prefix="$1"
    local level="$2"
    printf '%s0%s\n' "${prefix}" "${level}"
}

choose_mode_for_level() {
    local npos="$1"
    local level="$2"
    local combos

    combos="$(binomial_int "${npos}" "${level}")"
    if (( combos <= EXACT_MAX_COMBOS )); then
        printf 'exact;%s\n' "${combos}"
    else
        printf 'uniform;%s\n' "${combos}"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --levels)
            [[ $# -ge 2 ]] || die "missing value after $1"
            LEVEL_SPEC="$2"
            shift 2
            ;;
        --total-sites)
            [[ $# -ge 2 ]] || die "missing value after $1"
            TOTAL_SITES="$2"
            shift 2
            ;;
        -T|--temperature)
            [[ $# -ge 2 ]] || die "missing value after $1"
            TEMPERATURE="$2"
            shift 2
            ;;
        -C|--samples)
            [[ $# -ge 2 ]] || die "missing value after $1"
            UNIFORM_SAMPLES="$2"
            shift 2
            ;;
        --exact-max-combos)
            [[ $# -ge 2 ]] || die "missing value after $1"
            EXACT_MAX_COMBOS="$2"
            shift 2
            ;;
        -s|--seed)
            [[ $# -ge 2 ]] || die "missing value after $1"
            SEED="$2"
            shift 2
            ;;
        --sod-bin)
            [[ $# -ge 2 ]] || die "missing value after $1"
            SOD_BIN="$2"
            shift 2
            ;;
        --scripts-dir)
            [[ $# -ge 2 ]] || die "missing value after $1"
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        --template-gin)
            [[ $# -ge 2 ]] || die "missing value after $1"
            TEMPLATE_GIN="$2"
            shift 2
            ;;
        --no-template-gin)
            TEMPLATE_GIN=""
            shift
            ;;
        --protocol|--protocole)
            [[ $# -ge 2 ]] || die "missing value after $1"
            PROTOCOL_VERSION="$2"
            shift 2
            ;;
        --protocol=*|--protocole=*)
            PROTOCOL_VERSION="${1#*=}"
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --reset-summary)
            RESET_SUMMARY=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unrecognized argument: $1"
            ;;
    esac
done

[[ -x "${SOD_BIN}" ]] || die "cannot execute sod at ${SOD_BIN}"
looks_like_scripts_dir "${SCRIPTS_DIR}" || die "scripts directory looks incomplete: ${SCRIPTS_DIR}"
[[ -f INSOD ]] || die "INSOD was not found in $(pwd)"
[[ -f SGO ]] || die "SGO was not found in $(pwd)"
[[ "${TEMPERATURE}" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "invalid temperature: ${TEMPERATURE}"
[[ "${UNIFORM_SAMPLES}" =~ ^[0-9]+$ ]] || die "invalid sample count: ${UNIFORM_SAMPLES}"
[[ "${EXACT_MAX_COMBOS}" =~ ^[0-9]+$ ]] || die "invalid exact threshold: ${EXACT_MAX_COMBOS}"
[[ "${SEED}" =~ ^-?[0-9]+$ ]] || die "invalid seed: ${SEED}"

export SOD_SCRIPTS="${SCRIPTS_DIR}"

if [[ -z "${TOTAL_SITES}" ]]; then
    info "Inferring npos from sod exact --just-outsod -N 0"
    TOTAL_SITES="$(infer_total_sites)" || die "failed to infer the number of substitutable sites"
fi
[[ "${TOTAL_SITES}" =~ ^[0-9]+$ ]] || die "invalid total_sites value: ${TOTAL_SITES}"

if [[ -z "${LEVEL_SPEC}" ]]; then
    if (( TOTAL_SITES <= 10 )); then
        die "auto level window requires npos > 10; pass --levels explicitly"
    fi
    auto_end=$(( TOTAL_SITES - 5 ))
    if (( auto_end > 10 )); then
        auto_end=10
    fi
    if (( auto_end < 5 )); then
        die "no default calibrable window found; pass --levels explicitly"
    fi
    LEVEL_SPEC="5:${auto_end}"
fi

declare -a LEVELS
expand_level_spec "${LEVEL_SPEC}" LEVELS
[[ "${#LEVELS[@]}" -gt 0 ]] || die "no levels to process"

if [[ "${RESET_SUMMARY}" -eq 1 ]]; then
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        info "[dry-run] would remove CALIBRATION_COEFFICIENTS_SUMMARY.csv and CALIBRATION_BLEND_OVERRIDES_SUMMARY.csv"
    else
        rm -f CALIBRATION_COEFFICIENTS_SUMMARY.csv CALIBRATION_BLEND_OVERRIDES_SUMMARY.csv
    fi
fi

info "Calculation directory: $(pwd)"
info "sod: ${SOD_BIN}"
info "scripts dir: ${SCRIPTS_DIR}"
info "npos: ${TOTAL_SITES}"
info "levels: ${LEVEL_SPEC}"
info "temperature for MC: ${TEMPERATURE} K"
info "uniform samples: ${UNIFORM_SAMPLES}"
info "exact max combos: ${EXACT_MAX_COMBOS}"
info "protocol: ${PROTOCOL_VERSION}"
if [[ -n "${TEMPLATE_GIN}" ]]; then
    info "Template payload: ${TEMPLATE_GIN}"
else
    info "Template payload: none"
fi
info ""

processed=0
skipped=0

for level in "${LEVELS[@]}"; do
    if (( level <= 4 || level >= TOTAL_SITES - 4 )); then
        info "Skipping N=${level}: not in the dual-calibration window (needs N>4 and holes>4)."
        skipped=$((skipped + 1))
        continue
    fi

    IFS=';' read -r mode combos <<< "$(choose_mode_for_level "${TOTAL_SITES}" "${level}")"
    if [[ "${mode}" == "exact" ]]; then
        level_dir="$(format_level_dir "x" "${level}")"
    else
        level_dir="$(format_level_dir "mc" "${level}")"
    fi
    report_file="${level_dir}/CALIBRATION_COEFFICIENTS.txt"

    if [[ -f "${report_file}" && "${FORCE}" -eq 0 ]]; then
        info "Skipping N=${level}: found ${report_file}"
        skipped=$((skipped + 1))
        continue
    fi

    info "N=${level}  combinations=${combos}  mode=${mode}"
    if [[ "${mode}" == "exact" ]]; then
        cmd=( "${SOD_BIN}" exact -N "${level}" --protocol "${PROTOCOL_VERSION}" )
        if [[ -n "${TEMPLATE_GIN}" ]]; then
            cmd+=( --template-gin "${TEMPLATE_GIN}" )
        else
            cmd+=( --no-template-gin )
        fi
    else
        cmd=( "${SOD_BIN}" mc -N "${level}" -a uniform -C "${UNIFORM_SAMPLES}" -T "${TEMPERATURE}" -s "${SEED}" --force-mc --protocol "${PROTOCOL_VERSION}" )
        if [[ -n "${TEMPLATE_GIN}" ]]; then
            cmd+=( --template-gin "${TEMPLATE_GIN}" )
        else
            cmd+=( --no-template-gin )
        fi
    fi

    run_cmd "${cmd[@]}"
    processed=$((processed + 1))
    info ""
done

info "Finished."
info "Processed levels: ${processed}"
info "Skipped levels: ${skipped}"
info "Coefficient summary: $(pwd)/CALIBRATION_COEFFICIENTS_SUMMARY.csv"
info "Blend summary: $(pwd)/CALIBRATION_BLEND_OVERRIDES_SUMMARY.csv"
