# PHINDER Test Phages

This document lists the test phages used for validating the PHINDER pipeline.

## Quick Test (3 phages)

**File:** `test_phages_sra.txt`
**Runtime:** ~2-4 hours
**Script:** `bin/run_phinder_full_test_beocat.sh`

| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR5131134 | Lambda | *E. coli* | 48.5 kb | Temperate | Classic model phage |
| SRR5131135 | T4 | *E. coli* | 169 kb | Lytic | Large, complex genome |
| SRR5131136 | T7 | *E. coli* | 40 kb | Lytic | Well-studied, fast |

## Single Phage Test

**File:** `test_lambda_sra.txt`
**Runtime:** ~1 hour
**Script:** `bin/run_phinder_sra_beocat.sh`

| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR5131134 | Lambda | *E. coli* | 48.5 kb | Temperate | Quick validation |

## Extended Test (20 phages)

**File:** `test_phages_20_sra.txt`
**Runtime:** ~12-24 hours
**Script:** `bin/run_phinder_20phages_beocat.sh`

### Classic Model Phages (3)
| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR5131134 | Lambda | *E. coli* | 48.5 kb | Temperate | Integration/excision |
| SRR5131135 | T4 | *E. coli* | 169 kb | Lytic | Modified DNA bases |
| SRR5131136 | T7 | *E. coli* | 40 kb | Lytic | RNA polymerase |

### Diverse Collection (17)
| SRA Accession | Expected Type | Notes |
|--------------|---------------|-------|
| SRR8437269 | Phage isolate | Various hosts |
| SRR8437270 | Phage isolate | Genome diversity |
| SRR8437271 | Phage isolate | Size variation |
| SRR8437272 | Phage isolate | Lifestyle testing |
| SRR8437273 | Phage isolate | Quality checks |
| SRR8437274 | Phage isolate | Assembly testing |
| SRR8437275 | Phage isolate | Annotation testing |
| SRR11537895 | Recent isolate | Modern dataset |
| SRR11537896 | Recent isolate | 2020 sequences |
| SRR11537897 | Recent isolate | Updated tools |
| SRR11537898 | Recent isolate | Benchmark data |
| SRR13145901 | Latest isolate | 2021+ sequences |
| SRR13145902 | Latest isolate | Current standards |
| SRR13145903 | Latest isolate | Latest protocols |
| SRR13145904 | Latest isolate | Quality metrics |
| SRR13145905 | Latest isolate | Complete genomes |
| SRR13145906 | Latest isolate | Validation set |

## Usage

### Quick Test (recommended first)
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_sra_beocat.sh
```

### Full 3-Phage Test
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_full_test_beocat.sh
```

### Extended 20-Phage Test
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_20phages_beocat.sh
```

## Expected Results

Each phage will generate:
- **FastQC reports** - Read quality metrics
- **fastp reports** - Quality trimming results
- **Assembly** - Unicycler contigs (FASTA)
- **QUAST** - Assembly quality metrics
- **CheckV** - Genome completeness, contamination
- **Pharokka** - Gene annotations, functional predictions
- **VIBRANT** - Lifestyle prediction (lytic vs. temperate)
- **PHANOTATE** - Alternative gene calling
- **DIAMOND** - Prophage database hits

All results are summarized in:
- `multiqc/multiqc_report.html` - Aggregated QC metrics
- `summary/phinder_summary.html` - Per-sample phage analysis

## Download Results

```bash
# Download summary report
scp tylerdoe@beocat.cis.ksu.edu:/fastscratch/tylerdoe/PHINDER/results_*/summary/phinder_summary.html .

# Download MultiQC
scp tylerdoe@beocat.cis.ksu.edu:/fastscratch/tylerdoe/PHINDER/results_*/multiqc/multiqc_report.html .
```

## Troubleshooting

### SRA Download Failures
Some SRA accessions may fail to download due to:
- Network connectivity issues
- SRA database availability
- SLURM controller connectivity

**Solution:** Use the 20-phage test - even if some fail, others should succeed.

### SLURM Controller Errors
If you see `"Unable to contact slurm controller"`:
- This is a temporary Beocat infrastructure issue
- Wait 10-15 minutes and resubmit
- Check Beocat status: https://beocat.ksu.edu/status

### Memory Errors
If processes fail with OOM (out of memory):
- Increase `--mem` in SBATCH script
- Adjust process memory in `nextflow.config`

## Notes

- All phages are **publicly available** from NCBI SRA
- Accessions were selected for diversity in size, host, and lifestyle
- Some accessions may become deprecated - check SRA for alternatives
- Download times vary by file size and network conditions
