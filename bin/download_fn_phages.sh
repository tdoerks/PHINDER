#!/bin/bash
set -e
cd /fastscratch/tylerdoe/PHINDER
mkdir -p test_data

for acc in OR492271 OR492272 OR492273 OR492274 OR492275 OR492276; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=fasta&retmode=text" > test_data/${acc}.fasta
    count=$(grep -c '^>' test_data/${acc}.fasta || echo 0)
    if [[ "$count" -eq 0 ]]; then
        echo "ERROR: ${acc} got no sequence (rate-limited?) — check test_data/${acc}.fasta"
    else
        echo "Downloaded ${acc}: ${count} contig(s)"
    fi
    sleep 1
done

echo "Done. Files:"
ls -lh test_data/OR49*.fasta
