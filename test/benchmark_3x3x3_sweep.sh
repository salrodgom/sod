#!/bin/bash
# Benchmark sweep: modern SOD vs legacy combsod on quartz 3×3×3
# Measures enumeration-only time (FILER=-1) for each level N=0..6
set -e

MODERN_SOD="/home/salvador/coding/sod/sod/bin/sod"
LEGACY_COMBSOD="/home/salvador/software/sod/bin/combsod"
EXAMPLE_DIR="/home/salvador/coding/sod/sod/examples/spbemc_zeolites/quartz"
CSV_OUT="/home/salvador/coding/sod/sod/test/benchmark_3x3x3_sweep.csv"

MAX_N=6

echo "level,total_configs,inequivalent,modern_s,legacy_s" > "$CSV_OUT"

for N in $(seq 0 $MAX_N); do
    echo "--- Level N=$N ---"

    # Modern
    MDIR=$(mktemp -d /tmp/sod_sweep_m_XXXXXX)
    cp "$EXAMPLE_DIR"/{INSOD,SGO} "$MDIR/"
    sed -i 's/^2 2 2$/3 3 3/' "$MDIR/INSOD"
    sed -i 's/^11$/-1/' "$MDIR/INSOD"
    cd "$MDIR"
    # warm-up on first iteration
    if [ "$N" -eq 0 ]; then
        $MODERN_SOD comb -N 0 > /dev/null 2>&1
        rm -rf n00
    fi
    T_MODERN=$( { time $MODERN_SOD comb -N "$N" > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    # parse time output (format: 0m0.123s)
    T_MODERN=$(echo "$T_MODERN" | sed 's/m/*60+/;s/s//' | bc -l)
    TOTAL=$(grep "Total number" "$MDIR/n$(printf '%02d' $N)/OUTSOD" 2>/dev/null | head -1 || echo "")
    # Get total and inequiv from OUTSOD header
    OUTSOD_FILE="$MDIR/n$(printf '%02d' $N)/OUTSOD"
    if [ ${#N} -gt 2 ]; then
        OUTSOD_FILE="$MDIR/n0$(printf '%02d' $N)/OUTSOD"
    fi
    # Read from sod comb output instead
    OUT=$($MODERN_SOD comb -N "$N" 2>&1)
    TOTAL_CONFIGS=$(echo "$OUT" | grep "Total number of configurations" | awk '{print $NF}')
    INEQUIV=$(echo "$OUT" | grep "Number of inequivalent" | awk '{print $NF}')
    rm -rf "$MDIR"

    # Legacy
    LDIR=$(mktemp -d /tmp/sod_sweep_l_XXXXXX)
    cp "$EXAMPLE_DIR"/{INSOD,SGO} "$LDIR/"
    sed -i 's/^2 2 2$/3 3 3/' "$LDIR/INSOD"
    sed -i "/^# nsubs:/,/^[^#]/{s/^0$/$N/}" "$LDIR/INSOD"
    cd "$LDIR"
    T_LEGACY=$( { time $LEGACY_COMBSOD > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    T_LEGACY=$(echo "$T_LEGACY" | sed 's/m/*60+/;s/s//' | bc -l)
    rm -rf "$LDIR"

    echo "  total=$TOTAL_CONFIGS inequiv=$INEQUIV modern=${T_MODERN}s legacy=${T_LEGACY}s"
    echo "$N,$TOTAL_CONFIGS,$INEQUIV,$T_MODERN,$T_LEGACY" >> "$CSV_OUT"
done

echo ""
echo "Results written to: $CSV_OUT"
cat "$CSV_OUT"
