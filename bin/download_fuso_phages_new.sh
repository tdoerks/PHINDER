#!/usr/bin/env bash
# Download the 9 new Fusobacterium phages not in the original OR492271-276 test set
# Output: /fastscratch/tylerdoe/PHINDER/test_data/
# Run from Beocat

set -euo pipefail

OUTDIR="/fastscratch/tylerdoe/PHINDER/test_data"
mkdir -p "$OUTDIR"

declare -A PHAGES=(
    ["Fnu1"]="NC_055035.1"
    ["FNU2"]="OQ808963.1"
    ["FNU3"]="OQ808965.1"
    ["TCUFN6"]="OP004060.1"
    ["TCUFN3"]="PQ464787.1"
    ["JDFnp1"]="ON464759.1"
    ["JDFnp6"]="ON464763.1"
    ["JDFnp4"]="ON464764.1"
    ["JDFnp7"]="PV175345.1"
)

echo "Downloading 9 new Fusobacterium phages..."
echo ""

for NAME in "${!PHAGES[@]}"; do
    ACC="${PHAGES[$NAME]}"
    OUTFILE="$OUTDIR/${NAME}.fasta"

    if [[ -f "$OUTFILE" ]] && [[ $(grep -c '^>' "$OUTFILE" 2>/dev/null) -gt 0 ]]; then
        echo "  SKIP: $NAME ($ACC) — already exists"
        continue
    fi

    echo -n "  Downloading $NAME ($ACC)... "
    efetch -db nuccore -id "$ACC" -format fasta > "$OUTFILE"
    sleep 1  # be polite to NCBI

    CONTIGS=$(grep -c '^>' "$OUTFILE" 2>/dev/null || echo 0)
    SIZE=$(awk '/^[^>]/{s+=length($0)} END{print s}' "$OUTFILE")
    if [[ "$CONTIGS" -gt 0 ]]; then
        echo "OK (${CONTIGS} contig, ~${SIZE} bp)"
    else
        echo "ERROR — 0 contigs, check $OUTFILE"
    fi
done

echo ""
echo "Done. Verifying all files:"
for NAME in "${!PHAGES[@]}"; do
    OUTFILE="$OUTDIR/${NAME}.fasta"
    CONTIGS=$(grep -c '^>' "$OUTFILE" 2>/dev/null || echo 0)
    echo "  $NAME: $CONTIGS contig(s)"
done
