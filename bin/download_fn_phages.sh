#!/bin/bash
set -e
cd /fastscratch/tylerdoe/PHINDER
mkdir -p test_data

for acc in OR492271 OR492272 OR492273 OR492274 OR492275 OR492276; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=fasta&retmode=text" > test_data/${acc}.fasta
    echo "Downloaded ${acc}: $(grep -c '^>' test_data/${acc}.fasta) contig(s)"
done

echo "Done. Files:"
ls -lh test_data/OR49*.fasta
