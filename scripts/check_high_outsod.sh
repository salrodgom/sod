#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  check_high_outsod.sh [--npos N] [--max-hole-order M] [base_dir]

Checks the high-side OUTSOD files used by the complementary expansion.

Arguments:
  base_dir            Directory containing nXX/nXXX folders. Default: current directory.

Options:
  --npos N            Total number of substitution sites. If omitted, infer it from the
                      highest existing nXX/nXXX directory.
  --max-hole-order M  Highest hole order to inspect. Default: 4.
  -h, --help          Show this help message.

The script verifies, for each high-side OUTSOD:
  - expected number of numeric fields per configuration line
  - repeated Ge indices inside a configuration
  - hole count consistent with npos and the target hole order
EOF
}

find_cluster_dir() {
  local base_dir="$1"
  local level="$2"
  local candidate

  printf -v candidate "n%02d" "$level"
  if [[ -f "$base_dir/$candidate/OUTSOD" || -f "$base_dir/$candidate/ENERGIES" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  printf -v candidate "n%03d" "$level"
  if [[ -f "$base_dir/$candidate/OUTSOD" || -f "$base_dir/$candidate/ENERGIES" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  candidate="n${level}"
  if [[ -f "$base_dir/$candidate/OUTSOD" || -f "$base_dir/$candidate/ENERGIES" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  return 1
}

infer_npos() {
  local base_dir="$1"
  local path dir name level max_level

  max_level=-1
  shopt -s nullglob
  for path in "$base_dir"/n*; do
    [[ -d "$path" ]] || continue
    dir=$(basename "$path")
    name="${dir#n}"
    [[ "$name" =~ ^[0-9]+$ ]] || continue
    level=$((10#$name))
    if (( level > max_level )); then
      max_level=$level
    fi
  done
  shopt -u nullglob

  if (( max_level < 0 )); then
    return 1
  fi

  printf '%d\n' "$max_level"
}

check_one_outsod() {
  local file="$1"
  local npos="$2"
  local ge_count="$3"
  local hole_count="$4"

  awk -v npos="$npos" -v ge_count="$ge_count" -v hole_count="$hole_count" '
    BEGIN {
      expected_nf = ge_count + 2
      bad = 0
      checked = 0
    }
    NF > 0 && $1 !~ /#/ && $0 !~ /configuration|substitutions/ {
      checked++
      if (NF != expected_nf) {
        printf("  bad field count at line %d: got %d, expected %d\n", NR, NF, expected_nf)
        bad = 1
        next
      }

      delete seen
      dup = 0
      for (i = 3; i <= NF; i++) {
        if ($i < 1 || $i > npos) {
          printf("  out-of-range index at line %d: %s (expected 1..%d)\n", NR, $i, npos)
          bad = 1
        }
        if (seen[$i]++) dup = 1
      }
      if (dup) {
        printf("  duplicate indices at line %d: %s\n", NR, $0)
        bad = 1
      }

      holes = 0
      for (i = 1; i <= npos; i++) {
        if (!(i in seen)) holes++
      }
      if (holes != hole_count) {
        printf("  bad hole count at line %d: got %d, expected %d\n", NR, holes, hole_count)
        bad = 1
      }
    }
    END {
      if (checked == 0) {
        print "  no numeric configuration lines found"
        exit 2
      }
      if (bad) exit 1
    }
  ' "$file"
}

main() {
  local base_dir="."
  local npos=""
  local max_hole_order=4
  local dir ge_count hole_count file overall_ok=0

  while (($# > 0)); do
    case "$1" in
      --npos)
        npos="${2:-}"
        shift 2
        ;;
      --max-hole-order)
        max_hole_order="${2:-}"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        base_dir="$1"
        shift
        ;;
    esac
  done

  if [[ -z "$npos" ]]; then
    if ! npos=$(infer_npos "$base_dir"); then
      echo "Error: could not infer npos from $base_dir" >&2
      exit 1
    fi
  fi

  echo "Checking high-side OUTSOD files in $base_dir"
  echo "  inferred/selected npos = $npos"
  echo "  max hole order         = $max_hole_order"

  for ((hole_count = 1; hole_count <= max_hole_order; hole_count++)); do
    ge_count=$((npos - hole_count))
    if (( ge_count < 0 )); then
      break
    fi
    if ! dir=$(find_cluster_dir "$base_dir" "$ge_count"); then
      echo "[$hole_count hole(s)] missing directory for level $ge_count"
      continue
    fi
    file="$base_dir/$dir/OUTSOD"
    if [[ ! -f "$file" ]]; then
      echo "[$hole_count hole(s)] missing $file"
      continue
    fi

    echo "[$hole_count hole(s)] checking $file"
    if check_one_outsod "$file" "$npos" "$ge_count" "$hole_count"; then
      echo "  OK"
    else
      overall_ok=1
    fi
  done

  exit "$overall_ok"
}

main "$@"
