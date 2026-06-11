#!/bin/bash
# Quick CheckV smoke test using the patched utility.py (bind-mounted).
# Confirms the completeness DIAMOND step clears [3/8] after the --db fix.
set -euo pipefail

REPO="${REPO:-/fastscratch/tylerdoe/PHINDER}"
DB="${CHECKV_DB:-/fastscratch/tylerdoe/databases/checkv-db-v1.5}"
PATCH="$REPO/assets/checkv_utility.py"
TARGET="/usr/local/lib/python3.10/site-packages/checkv/utility.py"
IMG="docker://quay.io/biocontainers/checkv:1.0.2--pyhdfd78af_0"

# Find an assembly produced by an earlier run
ASM=$(find "$REPO/work_3phages" -name 'SRR5131134_assembly.fasta' 2>/dev/null | head -1)
[ -z "$ASM" ] && ASM=$(find "$REPO/work_3phages" -name '*_assembly.fasta' 2>/dev/null | head -1)
[ -z "$ASM" ] && { echo "ERROR: no *_assembly.fasta found under $REPO/work_3phages"; exit 1; }

[ -f "$PATCH" ] || { echo "ERROR: patched utility not found at $PATCH"; exit 1; }

OUT="/tmp/checkv_smoke_$$"
rm -rf "$OUT"
echo "Assembly : $ASM"
echo "DB       : $DB"
echo "Patch    : $PATCH"
echo "Output   : $OUT"
echo "----------------------------------------"

apptainer exec --bind "$PATCH:$TARGET" "$IMG" \
    checkv end_to_end "$ASM" "$OUT" -t 4 -d "$DB"

echo "========================================"
echo "=== CHECKV OK (completeness passed) -> $OUT ==="
