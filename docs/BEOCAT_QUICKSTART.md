# PHINDER Quick Start Guide for Beocat

Complete guide to running PHINDER on Kansas State's Beocat HPC cluster.

## Prerequisites

### 1. Clone PHINDER on Beocat

```bash
# SSH to Beocat
ssh <your-eid>@beocat.cis.ksu.edu

# Navigate to your workspace
cd /homes/$USER

# Clone repository
git clone https://github.com/tdoerks/PHINDER.git
cd PHINDER
```

### 2. Set Up Databases

PHINDER requires three databases. See [DATABASE_SETUP.md](DATABASE_SETUP.md) for details.

**Quick setup if databases already exist:**

```bash
# Edit nextflow.config to point to your database paths
nano nextflow.config
```

Update these lines:
```groovy
checkv_db    = '/homes/tylerdoe/databases/checkv-db-v1.5'
prophage_db  = '/homes/tylerdoe/databases/prophage_db.dmnd'
pharokka_db  = null  // Auto-downloads
```

### 3. Download Test Data

```bash
# Load SRA Toolkit
module load sra-toolkit

# Quick test (Lambda phage only, ~10 min)
./bin/download_test_phages.sh quick

# OR Full test (Lambda + T4 + T7, ~30 min)
./bin/download_test_phages.sh full
```

This downloads phage sequencing data and creates a samplesheet automatically.

## Running PHINDER

### Option 1: Using the Submission Script (Recommended)

```bash
# Submit to SLURM
sbatch bin/run_phinder_beocat.sh

# Check job status
squeue -u $USER

# Monitor progress
tail -f phinder_*.log
```

### Option 2: Interactive Session

```bash
# Start interactive session
salloc --partition=ksu-gen-bio.q --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=4:00:00

# Load Nextflow
module load Nextflow/24.04.4

# Run PHINDER
nextflow run main.nf \
    --input test_samplesheet.csv \
    --input_mode reads \
    --outdir results \
    -profile slurm \
    -resume

# Exit interactive session when done
exit
```

### Option 3: Direct Command Line

```bash
# From login node (Nextflow manages SLURM submissions)
module load Nextflow/24.04.4

nextflow run main.nf \
    --input test_samplesheet.csv \
    --input_mode reads \
    --outdir results \
    -profile slurm \
    -resume \
    -with-report results/report.html \
    -with-timeline results/timeline.html
```

## Understanding the Output

### Directory Structure

```
results/
├── fastqc/                      # Raw read quality
├── fastp/                       # Trimmed reads & QC
├── assemblies/                  # Assembled genomes
│   ├── lambda_assembly.fasta   # Single phage genome
│   ├── T4_assembly.fasta
│   └── T7_assembly.fasta
├── quast/                       # Assembly metrics
├── checkv/                      # Quality assessment
│   └── *_checkv/
│       └── quality_summary.tsv  # Completeness scores
├── pharokka/                    # Annotations
│   └── *_pharokka/
│       ├── *.gbk               # GenBank format
│       ├── *.gff               # Gene features
│       └── *_cds_functions.tsv # Gene functions
├── vibrant/                     # Lifestyle prediction
├── diamond_prophage/            # Taxonomy
├── phanotate/                   # Gene predictions
└── multiqc/
    └── multiqc_report.html      # ⭐ START HERE
```

### Key Output Files

1. **MultiQC Report** (`multiqc/multiqc_report.html`)
   - Overall quality metrics
   - Assembly statistics
   - Download to your computer to view

2. **Assembly** (`assemblies/*_assembly.fasta`)
   - Your phage genome(s)
   - Should be circular for complete genomes

3. **CheckV Quality** (`checkv/*/quality_summary.tsv`)
   - Completeness percentage
   - Contamination detection
   - Quality tier (Complete/High/Medium/Low)

4. **Pharokka Annotation** (`pharokka/*/`)
   - GenBank file (.gbk) for genome browsers
   - GFF file (.gff) for annotations
   - Functions table (.tsv) for genes

## Expected Results for Test Data

### Lambda Phage
- **Assembly size:** ~48.5 kb
- **Contigs:** 1 circular
- **Quality:** Complete
- **Lifestyle:** Temperate (lysogenic)
- **Genes:** ~70 CDSs
- **Runtime:** ~10-15 minutes

### T4 Phage
- **Assembly size:** ~169 kb
- **Contigs:** 1 circular
- **Quality:** Complete
- **Lifestyle:** Lytic
- **Genes:** ~280 CDSs
- **Runtime:** ~10-15 minutes

### T7 Phage
- **Assembly size:** ~40 kb
- **Contigs:** 1 circular
- **Quality:** Complete
- **Lifestyle:** Lytic
- **Genes:** ~55 CDSs
- **Runtime:** ~10-15 minutes

## Analyzing Your Own Data

### 1. Prepare Samplesheet

Create a CSV file with your sample information:

**For paired-end reads:**
```csv
sample,read1,read2
myPhage1,/path/to/phage1_R1.fastq.gz,/path/to/phage1_R2.fastq.gz
myPhage2,/path/to/phage2_R1.fastq.gz,/path/to/phage2_R2.fastq.gz
```

**For pre-assembled genomes:**
```csv
sample,assembly
myPhage1,/path/to/phage1.fasta
myPhage2,/path/to/phage2.fasta
```

### 2. Run Pipeline

```bash
# Edit submission script
nano bin/run_phinder_beocat.sh

# Change SAMPLESHEET line to your file:
SAMPLESHEET="my_samples.csv"

# Submit
sbatch bin/run_phinder_beocat.sh
```

### 3. Monitor Progress

```bash
# Check SLURM queue
squeue -u $USER

# Watch log in real-time
tail -f phinder_*.log

# Check Nextflow progress
tail -f .nextflow.log
```

## Advanced Options

### Skip Steps

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --skip_checkv \
    --skip_vibrant \
    --outdir results \
    -profile slurm
```

### Use Pre-assembled Genomes

```bash
nextflow run main.nf \
    --input assemblies.csv \
    --input_mode assembly \
    --outdir results \
    -profile slurm
```

### Adjust Resources

Edit `nextflow.config` to increase memory/CPUs for large phages:

```groovy
withName: UNICYCLER {
    cpus = 16
    memory = 64.GB
    time = 24.h
}
```

## Troubleshooting

### Job Fails Immediately

**Check:**
1. Database paths in `nextflow.config`
2. Input files exist: `ls -lh test_data/`
3. Samplesheet format correct

### Out of Memory Errors

**Solution:** Increase memory in SLURM script:
```bash
#SBATCH --mem=16G  # Instead of 8G
```

### Container Pull Failures

**Solution:** Pre-pull containers:
```bash
apptainer pull docker://staphb/fastqc:0.12.1
apptainer pull docker://staphb/fastp:0.23.4
# etc...
```

### Databases Not Found

**Verify paths:**
```bash
ls -lh /homes/tylerdoe/databases/checkv-db-v1.5/
ls -lh /homes/tylerdoe/databases/prophage_db.dmnd
```

### Resume Failed Run

PHINDER supports resume! Just re-run with `-resume`:
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile slurm \
    -resume  # ← Skips completed steps
```

## Getting Help

1. **Check logs:**
   - `phinder_*.log` - SLURM output
   - `.nextflow.log` - Nextflow details
   - `work/` - Individual process logs

2. **Common issues:**
   - Database setup: See [DATABASE_SETUP.md](DATABASE_SETUP.md)
   - General usage: See [README.md](../README.md)

3. **Open an issue:**
   - https://github.com/tdoerks/PHINDER/issues

## Quick Reference

```bash
# Download test data
./bin/download_test_phages.sh quick

# Run test
sbatch bin/run_phinder_beocat.sh

# Check status
squeue -u $USER

# View results
firefox results/multiqc/multiqc_report.html
```

---

**Ready to analyze phages on Beocat!** 🧬🚀
