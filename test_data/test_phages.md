# Test Phage Datasets for PHINDER

## Well-Characterized Model Phages (SRP093616)

These are classic, well-studied phages with known genomes - perfect for pipeline validation:

### 1. Lambda Phage (Escherichia coli)
- **SRR Accession:** SRR5131134
- **Genome Size:** ~48.5 kb
- **Type:** Temperate (lysogenic)
- **Reference:** NC_001416
- **Platform:** Illumina MiSeq
- **Layout:** Paired-end

### 2. T4 Phage (Escherichia coli)
- **SRR Accession:** SRR5131135
- **Genome Size:** ~169 kb
- **Type:** Lytic
- **Reference:** NC_000866
- **Platform:** Illumina MiSeq
- **Layout:** Paired-end

### 3. T7 Phage (Escherichia coli)
- **SRR Accession:** SRR5131136
- **Genome Size:** ~40 kb
- **Type:** Lytic
- **Reference:** NC_001604
- **Platform:** Illumina MiSeq
- **Layout:** Paired-end

## Recent Phage Isolates (Microbacterium phages)

Modern phage sequencing from SEA-PHAGES program:

### 4. Microbacterium phage Swervy
- **SRR Accession:** SRR17754903
- **Genome Size:** ~50-60 kb (estimated)
- **Type:** Novel isolate
- **Platform:** Illumina MiSeq
- **Year:** 2021

### 5. Microbacterium phage Fullmetal
- **SRR Accession:** SRR23711508
- **Genome Size:** ~50-70 kb (estimated)
- **Type:** Novel isolate
- **Platform:** Illumina MiSeq
- **Year:** 2023

## Download Commands

### Using SRA Toolkit (fasterq-dump)
```bash
# Download Lambda phage
fasterq-dump SRR5131134 -O test_data/ -e 8

# Download T4 phage
fasterq-dump SRR5131135 -O test_data/ -e 8

# Download T7 phage
fasterq-dump SRR5131136 -O test_data/ -e 8
```

### Using Nextflow (future feature)
```bash
# Will add SRA download module similar to COMPASS
nextflow run main.nf \
    --input_mode sra \
    --sra_accessions "SRR5131134,SRR5131135,SRR5131136" \
    --outdir results_test
```

## Recommended Test Strategy

### Quick Test (1 sample, ~10-15 min)
```bash
# Download Lambda (smallest)
fasterq-dump SRR5131134 -O test_data/

# Compress
gzip test_data/SRR5131134_*.fastq

# Create samplesheet
cat > test_samplesheet.csv << EOF
sample,read1,read2
lambda,test_data/SRR5131134_1.fastq.gz,test_data/SRR5131134_2.fastq.gz
EOF

# Run PHINDER
nextflow run main.nf \
    --input test_samplesheet.csv \
    --input_mode reads \
    --outdir test_results \
    -profile slurm
```

### Full Test (3 samples, ~30-45 min)
```bash
# Download all three classic phages
for srr in SRR5131134 SRR5131135 SRR5131136; do
    fasterq-dump $srr -O test_data/ -e 8
    gzip test_data/${srr}_*.fastq
done

# Create samplesheet
cat > test_samplesheet_full.csv << EOF
sample,read1,read2
lambda,test_data/SRR5131134_1.fastq.gz,test_data/SRR5131134_2.fastq.gz
T4,test_data/SRR5131135_1.fastq.gz,test_data/SRR5131135_2.fastq.gz
T7,test_data/SRR5131136_1.fastq.gz,test_data/SRR5131136_2.fastq.gz
EOF

# Run PHINDER
nextflow run main.nf \
    --input test_samplesheet_full.csv \
    --input_mode reads \
    --outdir test_results_full \
    -profile slurm
```

## Expected Results

### Lambda Phage
- **Assembly:** Single circular contig, ~48.5 kb
- **CheckV Quality:** Complete
- **VIBRANT:** Temperate/lysogenic
- **Genes:** ~70 CDSs
- **Notable:** Integrase gene present

### T4 Phage
- **Assembly:** Single circular contig, ~169 kb
- **CheckV Quality:** Complete
- **VIBRANT:** Lytic
- **Genes:** ~280 CDSs
- **Notable:** Large genome, many structural genes

### T7 Phage
- **Assembly:** Single circular contig, ~40 kb
- **CheckV Quality:** Complete
- **VIBRANT:** Lytic
- **Genes:** ~55 CDSs
- **Notable:** Well-annotated, compact genome

## Validation Metrics

Pipeline is working correctly if:
- ✅ All samples assemble to single circular contigs
- ✅ Assembly sizes match expected (~48kb, ~169kb, ~40kb)
- ✅ CheckV reports "Complete" quality
- ✅ Pharokka annotates 50-280 genes
- ✅ VIBRANT correctly predicts lifestyle
- ✅ MultiQC report generates without errors

## Notes

- **Data Size:** Each SRR ~100-500 MB compressed
- **Runtime:** ~10-15 min per sample on HPC
- **Storage:** Keep test data for future regression testing
- **References:** Compare final assemblies to RefSeq genomes for validation
