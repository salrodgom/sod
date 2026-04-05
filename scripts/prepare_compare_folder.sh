#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_SOD_BIN="${REPO_ROOT}/bin/sod"

LEVEL_SPEC="-1"
SOD_BIN="${DEFAULT_SOD_BIN}"
FORCE_OUTSOD=0
FORCE_ENTROPY=0
FORCE_MC=0
SKIP_MC=0
LABEL=""
MC_TEMPERATURE=""
MC_SAMPLER="metropolis"
MC_SEED="-1"
MC_TEMPLATE_GIN=""
PROTOCOL_VERSION="${SOD_GULP_PROTOCOL_VERSION:-2.0}"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Prepares the current folder so it can later be used by:
  sod compare --system <folder> --reference <folder> -T <K>

What it does:
  1. Generates exact OUTSOD files with: sod exact --just-outsod
  2. Generates sod_entropy_summary.csv with: sod entropy
  3. Generates sod_ensemble_summary.csv with: sod mc when needed
  4. Writes compare_folder_status.txt with the folder status

Options:
  -N, --levels <spec>      Level specification passed to exact/entropy [default: -1]
  -T, --temperature <K>    Temperature for MC if sod_ensemble_summary.csv must be generated
  --sod-bin <path>         Path to the sod executable
  --label <text>           Optional label stored in compare_folder_status.txt
  -a, --sampler <name>     MC sampler [default: metropolis]
  -s, --seed <value>       MC seed [default: -1]
  --template-gin <file>    Passed through to sod mc
  --protocol <ver>         Passed through to sod mc [default: 2.0]
  --force-outsod           Regenerate exact OUTSOD files even if some already exist
  --force-entropy          Regenerate sod_entropy_summary.csv
  --force-mc               Regenerate sod_ensemble_summary.csv even if it already exists
  --skip-mc                Do not run sod mc
  -h, --help               Show this help

Examples:
  $(basename "$0") -T 1500
  $(basename "$0") -T 1500 -N -1 --template-gin template_payload.gin --protocol 2.0
  $(basename "$0") --sod-bin /path/to/bin/sod -T 1500 --force-outsod --force-mc
EOF
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

info() {
    printf '%s\n' "$*"
}

status_bool() {
    if [[ "$1" -eq 0 ]]; then
        printf 'no'
    else
        printf 'yes'
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -N|--levels)
            [[ $# -ge 2 ]] || die "missing value after $1"
            LEVEL_SPEC="$2"
            shift 2
            ;;
        -T|--temperature)
            [[ $# -ge 2 ]] || die "missing value after $1"
            MC_TEMPERATURE="$2"
            shift 2
            ;;
        --sod-bin)
            [[ $# -ge 2 ]] || die "missing value after $1"
            SOD_BIN="$2"
            shift 2
            ;;
        --label)
            [[ $# -ge 2 ]] || die "missing value after $1"
            LABEL="$2"
            shift 2
            ;;
        -a|--sampler)
            [[ $# -ge 2 ]] || die "missing value after $1"
            MC_SAMPLER="$2"
            shift 2
            ;;
        -s|--seed)
            [[ $# -ge 2 ]] || die "missing value after $1"
            MC_SEED="$2"
            shift 2
            ;;
        --template-gin)
            [[ $# -ge 2 ]] || die "missing value after $1"
            MC_TEMPLATE_GIN="$2"
            shift 2
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
        --force-outsod)
            FORCE_OUTSOD=1
            shift
            ;;
        --force-entropy)
            FORCE_ENTROPY=1
            shift
            ;;
        --force-mc)
            FORCE_MC=1
            shift
            ;;
        --skip-mc)
            SKIP_MC=1
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

[[ -x "${SOD_BIN}" ]] || die "cannot execute sod at ${SOD_BIN}. Build it first or pass --sod-bin."
[[ -f INSOD ]] || die "INSOD was not found in $(pwd)"
[[ -f SGO ]] || die "SGO was not found in $(pwd)"

if [[ -z "${LABEL}" ]]; then
    LABEL="$(basename "$(pwd)")"
fi

OUTSOD_FOUND=0
shopt -s nullglob
outsod_files=(OUTSOD_N[0-9][0-9][0-9][0-9] x0*/OUTSOD)
if [[ ${#outsod_files[@]} -gt 0 ]]; then
    OUTSOD_FOUND=1
fi
shopt -u nullglob

if [[ "${FORCE_OUTSOD}" -eq 1 || "${OUTSOD_FOUND}" -eq 0 ]]; then
    info "Generating exact OUTSOD files with exact --just-outsod (-N ${LEVEL_SPEC})"
    "${SOD_BIN}" exact --just-outsod -N "${LEVEL_SPEC}"
    OUTSOD_GENERATED=1
else
    info "Reusing existing OUTSOD files. Use --force-outsod to regenerate them."
    OUTSOD_GENERATED=0
fi

if [[ "${FORCE_ENTROPY}" -eq 1 || "${OUTSOD_GENERATED}" -eq 1 || ! -f sod_entropy_summary.csv ]]; then
    info "Generating sod_entropy_summary.csv with entropy (-N ${LEVEL_SPEC})"
    "${SOD_BIN}" entropy -N "${LEVEL_SPEC}"
    ENTROPY_GENERATED=1
else
    info "Reusing existing sod_entropy_summary.csv. Use --force-entropy to regenerate it."
    ENTROPY_GENERATED=0
fi

MC_READY=0
MC_GENERATED=0
if [[ "${FORCE_MC}" -eq 1 || ! -f sod_ensemble_summary.csv ]]; then
    if [[ "${SKIP_MC}" -eq 1 ]]; then
        info "Skipping MC generation because --skip-mc was requested."
    else
        [[ -n "${MC_TEMPERATURE}" ]] || die "MC generation requires -T/--temperature <K>"
        info "Generating sod_ensemble_summary.csv with mc (-T ${MC_TEMPERATURE}, -N ${LEVEL_SPEC}, -a ${MC_SAMPLER}, -s ${MC_SEED})"
        mc_cmd=( "${SOD_BIN}" mc -T "${MC_TEMPERATURE}" -s "${MC_SEED}" -N "${LEVEL_SPEC}" -a "${MC_SAMPLER}" --protocol "${PROTOCOL_VERSION}" )
        if [[ -n "${MC_TEMPLATE_GIN}" ]]; then
            mc_cmd+=( --template-gin "${MC_TEMPLATE_GIN}" )
        fi
        "${mc_cmd[@]}"
        MC_GENERATED=1
    fi
else
    info "Reusing existing sod_ensemble_summary.csv. Use --force-mc to regenerate it."
fi

if [[ -f sod_ensemble_summary.csv ]]; then
    MC_READY=1
    info "Found sod_ensemble_summary.csv"
else
    info "Warning: sod_ensemble_summary.csv is missing."
    info "The folder is not fully ready for compare until the MC summary exists."
fi

if [[ -f sod_entropy_summary.csv && "${MC_READY}" -eq 1 ]]; then
    READY_FOR_COMPARE=1
else
    READY_FOR_COMPARE=0
fi

cat > compare_folder_status.txt <<EOF
folder=$(pwd)
label=${LABEL}
levels=${LEVEL_SPEC}
sod_bin=${SOD_BIN}
outsod_generated=$(status_bool "${OUTSOD_GENERATED}")
entropy_generated=$(status_bool "${ENTROPY_GENERATED}")
mc_generated=$(status_bool "${MC_GENERATED}")
has_mc_summary=$(status_bool "${MC_READY}")
ready_for_compare=$(status_bool "${READY_FOR_COMPARE}")
mc_summary_file=$(pwd)/sod_ensemble_summary.csv
entropy_summary_file=$(pwd)/sod_entropy_summary.csv
EOF

info ""
info "Folder status written to compare_folder_status.txt"
if [[ "${READY_FOR_COMPARE}" -eq 1 ]]; then
    info "This folder is ready for sod compare."
else
    info "This folder is only partially prepared. You still need sod_ensemble_summary.csv from MC."
fi
info ""
info "Suggested compare command later:"
info "  ${SOD_BIN} compare --system /path/to/system --reference /path/to/reference -T 1500 -o compare_1500K"
