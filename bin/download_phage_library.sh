#!/usr/bin/env bash
# Download all non-Fusobacterium phages in the PHINDER phage library
# Fusobacterium phages (OR492271-276, Fnu1, FNU2/3, TCUFN3/6, JDFnp1/4/6/7)
# are handled by bin/download_fuso_phages_new.sh and bin/download_fn_phages.sh
#
# Run from: /fastscratch/tylerdoe/PHINDER/
# Output:   test_data/<sample_id>.fasta

set -euo pipefail

OUTDIR="/fastscratch/tylerdoe/PHINDER/test_data"
mkdir -p "$OUTDIR"

# Format: ["sample_id"]="accession"
declare -A PHAGES

# Porphyromonas (oral anaerobe)
PHAGES["Pg_phage019b"]="PP754931.1"
PHAGES["Pg_phage010a"]="PP754930.1"
PHAGES["Pg_phage007a"]="PP754929.1"
PHAGES["Pg_phage005b"]="PP754928.1"

# Clostridioides difficile (gut anaerobe / therapeutic interest)
PHAGES["Cd_phiCD27"]="NC_011398.1"
PHAGES["Cd_phiCD508"]="OR295560.1"
PHAGES["Cd_phiCDKH02"]="PP767789.1"
PHAGES["Cd_AR1086"]="OR397124.1"
PHAGES["Cd_AR1074"]="OR397123.1"
PHAGES["Cd_AR1075"]="OQ703261.1"
PHAGES["Cd_phiSemix9P1"]="KX905163.1"
PHAGES["Cd_phiCD418"]="MW512573.1"
PHAGES["Cd_phiCD08011"]="MW512572.1"
PHAGES["Cd_CD2301"]="MW512571.1"
PHAGES["Cd_CD1801"]="MW512570.1"
PHAGES["Cd_JD032"]="MK473382.1"
PHAGES["Cd_JD033"]="MT193276.1"
PHAGES["Cd_phiCDKH01"]="MN718463.1"
PHAGES["Cd_LIBA2945"]="MF547663.1"
PHAGES["Cd_LIBA6276"]="MF547662.1"
PHAGES["Cd_HainSaunders"]="CP103806.1"

# Klebsiella pneumoniae (nosocomial pathogen)
PHAGES["Kp_phiW14"]="NC_105446.1"
PHAGES["Kp_Solomon"]="NC_103145.1"
PHAGES["Kp_Penguinator"]="NC_100761.1"
PHAGES["Kp_Opt817"]="NC_105496.1"
PHAGES["Kp_P1"]="NC_129957.1"
PHAGES["Kp_NPat"]="NC_106517.1"
PHAGES["Kp_PWKp20"]="NC_105644.1"
PHAGES["Kp_B1"]="NC_104812.1"
PHAGES["Kp_IME268"]="NC_129956.1"
PHAGES["Kp_FairDinkum"]="NC_105416.1"

# Pseudomonas aeruginosa (ESKAPE pathogen)
PHAGES["Pa_VIPEUCMC01"]="NC_130827.1"
PHAGES["Pa_POR1"]="NC_093432.1"
PHAGES["Pa_UMP151"]="NC_074747.1"
PHAGES["Pa_BUCT566"]="NC_074664.1"
PHAGES["Pa_LC3I3"]="NC_074663.1"
PHAGES["Pa_PP9W2"]="NC_074662.1"
PHAGES["Pa_Epa33"]="NC_073676.1"
PHAGES["Pa_LKA5"]="NC_073675.1"
PHAGES["Pa_E220"]="NC_073673.1"
PHAGES["Pa_PA8"]="NC_073672.1"
PHAGES["Pa_PsCh"]="NC_073623.1"
PHAGES["Pa_PsIn"]="NC_073622.1"
PHAGES["Pa_Ps12"]="NC_073621.1"

# Staphylococcus aureus (MRSA / therapeutic interest)
PHAGES["Sa_PhageK"]="NC_005880.2"
PHAGES["Sa_IMESA119"]="NC_047732.1"
PHAGES["Sa_IMESA118"]="NC_047731.1"
PHAGES["Sa_IMESA2"]="NC_047730.1"
PHAGES["Sa_IMESA1"]="NC_047729.1"
PHAGES["Sa_A5W"]="NC_047728.1"
PHAGES["Sa_P4W"]="NC_047727.1"
PHAGES["Sa_MSA6"]="NC_047726.1"
PHAGES["Sa_Fi200W"]="NC_047725.1"
PHAGES["Sa_676Z"]="NC_047724.1"
PHAGES["Sa_A3R"]="NC_047723.1"
PHAGES["Sa_Staph1N"]="NC_047722.1"

# Salmonella enterica
PHAGES["Sm_UTK0007"]="NC_130781.1"
PHAGES["Sm_UTK0006"]="NC_111484.1"
PHAGES["Sm_GSP032"]="NC_112025.1"
PHAGES["Sm_F61"]="NC_108349.1"
PHAGES["Sm_falkor"]="NC_102106.1"
PHAGES["Sm_GEC_N5"]="NC_103213.1"
PHAGES["Sm_LmqsSP1"]="NC_104259.1"

# Classic E. coli (validation standards — known biology)
PHAGES["Ec_T4"]="NC_000866.4"
PHAGES["Ec_T7"]="NC_001604.1"
PHAGES["Ec_Lambda"]="NC_001416.1"
PHAGES["Ec_P1"]="NC_005856.1"
PHAGES["Ec_Mu"]="NC_000929.1"
PHAGES["Ec_K1E"]="NC_007637.1"
PHAGES["Ec_VIPECOOM03"]="NC_111760.1"

# Acinetobacter baumannii (jumbo myovirus)
PHAGES["Ab_ME3"]="NC_041884.1"

# ─────────────────────────────────────────────────────────────────────────────

TOTAL=${#PHAGES[@]}
COUNT=0
SKIPPED=0
FAILED=0

echo "========================================"
echo "  PHINDER Phage Library Download"
echo "  ${TOTAL} phages to download"
echo "  Output: ${OUTDIR}"
echo "========================================"
echo ""

for NAME in $(echo "${!PHAGES[@]}" | tr ' ' '\n' | sort); do
    ACC="${PHAGES[$NAME]}"
    OUTFILE="$OUTDIR/${NAME}.fasta"

    if [[ -f "$OUTFILE" ]] && [[ $(grep -c '^>' "$OUTFILE" 2>/dev/null || echo 0) -gt 0 ]]; then
        SKIPPED=$((SKIPPED+1))
        echo "  SKIP: $NAME ($ACC)"
        continue
    fi

    echo -n "  [$((COUNT+SKIPPED+1))/${TOTAL}] $NAME ($ACC)... "
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" \
        > "$OUTFILE"
    sleep 1  # be polite to NCBI

    CONTIGS=$(grep -c '^>' "$OUTFILE" 2>/dev/null || echo 0)
    if [[ "$CONTIGS" -gt 0 ]]; then
        SIZE=$(awk '/^[^>]/{s+=length($0)} END{print s}' "$OUTFILE")
        echo "OK (${SIZE} bp)"
        COUNT=$((COUNT+1))
    else
        echo "FAILED — 0 contigs"
        FAILED=$((FAILED+1))
        rm -f "$OUTFILE"
    fi
done

echo ""
echo "========================================"
echo "  Downloaded: ${COUNT}"
echo "  Skipped (already exist): ${SKIPPED}"
echo "  Failed: ${FAILED}"
echo "========================================"

if [[ $FAILED -gt 0 ]]; then
    echo "Re-run this script to retry failed downloads."
    exit 1
fi
