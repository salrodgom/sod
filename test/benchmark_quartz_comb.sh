#!/bin/bash
# Benchmark: modern SOD vs legacy SOD 0.62 — comb mode on quartz 2x2x2
# Compares wall-clock time for levels 0..5

set -e

MODERN_SOD="/home/salvador/coding/sod/sod/bin/sod"
LEGACY_COMBSOD="/home/salvador/software/sod/bin/combsod"
LEGACY_GENERSOD="/home/salvador/software/sod/bin/genersod"

EXAMPLE_DIR="/home/salvador/coding/sod/sod/examples/spbemc_zeolites/quartz"
LEVELS="0,1,2,3,4,5"

WORK_DIR=$(mktemp -d /tmp/sod_bench_XXXXXX)
echo "=== SOD comb benchmark: quartz 2x2x2, levels $LEVELS ==="
echo "Work directory: $WORK_DIR"
echo ""

# --- Prepare INSOD for legacy (nsubs=0 means full range in legacy) ---
# We need a modified INSOD with nsubs_min=0 nsubs_max=5 for a fair comparison.
# Legacy combsod reads "nsubs" and if 0, does all; otherwise it does that single level.
# Actually legacy supports nsubs_min nsubs_max on the nsubs line.
# Let's check and create appropriate INSODs.

# -----------------------------------------------------------------------
#  MODERN SOD
# -----------------------------------------------------------------------
echo "=============================="
echo "  MODERN SOD (sod comb)"
echo "=============================="

MODERN_DIR="$WORK_DIR/modern"
mkdir -p "$MODERN_DIR"
cp "$EXAMPLE_DIR/INSOD" "$MODERN_DIR/"
cp "$EXAMPLE_DIR/SGO" "$MODERN_DIR/"

cd "$MODERN_DIR"

echo ""
echo "Running: $MODERN_SOD comb -N $LEVELS"
echo ""

# Warm-up run (compile caches, etc.)
$MODERN_SOD comb -N 0 > /dev/null 2>&1
rm -rf n00

# Timed run
TIMEFORMAT='  Wall time: %3R s (user: %3U s, sys: %3S s)'
time $MODERN_SOD comb -N "$LEVELS" 2>&1

echo ""
echo "Output directories:"
ls -d n*/ 2>/dev/null | head -20
echo ""

# Count total configurations
MODERN_TOTAL=0
for d in n*/; do
    if [ -f "${d}OUTSOD" ]; then
        n=$(grep -c "^[[:space:]]*[0-9]" "${d}OUTSOD" 2>/dev/null || echo 0)
        level=$(basename "$d")
        echo "  $level: $n configurations"
        MODERN_TOTAL=$((MODERN_TOTAL + n))
    fi
done
echo "  TOTAL: $MODERN_TOTAL configurations"

# -----------------------------------------------------------------------
#  LEGACY SOD 0.62
# -----------------------------------------------------------------------
echo ""
echo "=============================="
echo "  LEGACY SOD 0.62 (combsod)"
echo "=============================="

LEGACY_DIR="$WORK_DIR/legacy"
mkdir -p "$LEGACY_DIR"
cp "$EXAMPLE_DIR/INSOD" "$LEGACY_DIR/"
cp "$EXAMPLE_DIR/SGO" "$LEGACY_DIR/"

cd "$LEGACY_DIR"

# Modify INSOD nsubs line to "0 5" (min=0, max=5) for legacy
# The nsubs line is the one after the "# nsubs:" comment
sed -i '/^# nsubs:/,/^[^#]/{s/^0$/0 5/}' INSOD

echo ""
echo "Running: combsod + genersod"
echo ""

# Warm-up
$LEGACY_COMBSOD > /dev/null 2>&1
rm -rf n* OUTSOD SUPERCELL EQMATRIX OPERATORS cSGO filer 2>/dev/null

# Timed run
time {
    $LEGACY_COMBSOD 2>&1
    # genersod generates the structure files (equivalent to what modern comb does)
    if [ -f filer ]; then
        $LEGACY_GENERSOD 2>&1
    fi
}

echo ""
echo "Output directories:"
ls -d n*/ 2>/dev/null | head -20
echo ""

# Count total configurations
LEGACY_TOTAL=0
for d in n*/; do
    if [ -f "${d}OUTSOD" ]; then
        n=$(grep -c "^[[:space:]]*[0-9]" "${d}OUTSOD" 2>/dev/null || echo 0)
        level=$(basename "$d")
        echo "  $level: $n configurations"
        LEGACY_TOTAL=$((LEGACY_TOTAL + n))
    fi
done
echo "  TOTAL: $LEGACY_TOTAL configurations"

# -----------------------------------------------------------------------
#  Summary
# -----------------------------------------------------------------------
echo ""
echo "=============================="
echo "  Cleanup"
echo "=============================="
echo "Benchmark data in: $WORK_DIR"
echo "To remove: rm -rf $WORK_DIR"
