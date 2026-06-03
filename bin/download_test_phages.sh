#!/bin/bash

#==============================================================================
# PHINDER Test Data Download Script
#==============================================================================
# Downloads well-characterized phage genomes from NCBI SRA for pipeline testing
#
# Usage:
#   ./bin/download_test_phages.sh [quick|full]
#
# Options:
#   quick - Download only Lambda phage (~100 MB, fast test)
#   full  - Download Lambda, T4, and T7 (~400 MB, comprehensive test)
#
# Requirements:
#   - SRA Toolkit (fasterq-dump)
#   - ~1-2 GB free space
#==============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default mode
MODE="${1:-quick}"

# Output directory
OUTDIR="test_data"
mkdir -p ${OUTDIR}

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  PHINDER Test Data Download${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Check for fasterq-dump
if ! command -v fasterq-dump &> /dev/null; then
    echo -e "${RED}ERROR: fasterq-dump not found!${NC}"
    echo "Please install SRA Toolkit:"
    echo "  module load sra-toolkit  # On HPC"
    echo "  conda install -c bioconda sra-tools  # With conda"
    exit 1
fi

# Function to download and process SRA data
download_phage() {
    local SRR=$1
    local NAME=$2
    local SIZE=$3

    echo -e "${YELLOW}Downloading ${NAME} (${SRR})...${NC}"
    echo "  Expected genome size: ${SIZE}"

    # Download
    if [ ! -f "${OUTDIR}/${SRR}_1.fastq" ]; then
        fasterq-dump ${SRR} -O ${OUTDIR}/ -e 8 -p
    else
        echo "  Reads already downloaded, skipping..."
    fi

    # Compress
    if [ ! -f "${OUTDIR}/${SRR}_1.fastq.gz" ]; then
        echo "  Compressing reads..."
        gzip -f ${OUTDIR}/${SRR}_*.fastq
    else
        echo "  Reads already compressed, skipping..."
    fi

    echo -e "${GREEN}  ✓ ${NAME} ready${NC}"
    echo ""
}

# Download based on mode
if [ "${MODE}" == "quick" ]; then
    echo "Mode: Quick test (Lambda phage only)"
    echo ""
    download_phage "SRR5131134" "Lambda phage" "~48.5 kb"

elif [ "${MODE}" == "full" ]; then
    echo "Mode: Full test (Lambda, T4, T7)"
    echo ""
    download_phage "SRR5131134" "Lambda phage" "~48.5 kb"
    download_phage "SRR5131135" "T4 phage" "~169 kb"
    download_phage "SRR5131136" "T7 phage" "~40 kb"

else
    echo -e "${RED}ERROR: Invalid mode '${MODE}'${NC}"
    echo "Usage: $0 [quick|full]"
    exit 1
fi

# Create samplesheet
echo -e "${YELLOW}Creating samplesheet...${NC}"

if [ "${MODE}" == "quick" ]; then
    cat > test_samplesheet.csv << EOF
sample,read1,read2
lambda,${OUTDIR}/SRR5131134_1.fastq.gz,${OUTDIR}/SRR5131134_2.fastq.gz
EOF
    SHEET="test_samplesheet.csv"

elif [ "${MODE}" == "full" ]; then
    cat > test_samplesheet_full.csv << EOF
sample,read1,read2
lambda,${OUTDIR}/SRR5131134_1.fastq.gz,${OUTDIR}/SRR5131134_2.fastq.gz
T4,${OUTDIR}/SRR5131135_1.fastq.gz,${OUTDIR}/SRR5131135_2.fastq.gz
T7,${OUTDIR}/SRR5131136_1.fastq.gz,${OUTDIR}/SRR5131136_2.fastq.gz
EOF
    SHEET="test_samplesheet_full.csv"
fi

echo -e "${GREEN}  ✓ Samplesheet created: ${SHEET}${NC}"
echo ""

# Summary
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Download Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Test data ready in: ${OUTDIR}/"
echo "Samplesheet: ${SHEET}"
echo ""
echo "Next steps:"
echo ""
echo "  1. Verify database paths in nextflow.config"
echo "  2. Run PHINDER:"
echo ""
if [ "${MODE}" == "quick" ]; then
    echo "     nextflow run main.nf \\"
    echo "         --input ${SHEET} \\"
    echo "         --input_mode reads \\"
    echo "         --outdir test_results \\"
    echo "         -profile slurm"
else
    echo "     nextflow run main.nf \\"
    echo "         --input ${SHEET} \\"
    echo "         --input_mode reads \\"
    echo "         --outdir test_results_full \\"
    echo "         -profile slurm"
fi
echo ""
echo "  3. Check results in test_results/multiqc/multiqc_report.html"
echo ""
