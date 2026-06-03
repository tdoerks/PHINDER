# PHINDER

**PHINDER**: PHage Isolate characterizatioN, Discovery & Evaluation Resource

A comprehensive Nextflow pipeline for end-to-end analysis of purified phage sequencing data.

## Overview

PHINDER is designed specifically for **purified phage isolates** from laboratory settings, providing a complete workflow from raw sequencing reads to publication-quality genome characterization.

**Key Features:**
- 🧬 **Complete Workflow**: Raw reads → Assembly → Quality assessment → Annotation → Classification
- 🎯 **Phage-Optimized**: Tools and parameters specifically tuned for phage biology
- 🔄 **Reproducible**: Containerized tools with Nextflow for consistent results
- 📊 **Comprehensive Output**: Quality metrics, annotations, classifications, and interactive reports
- 🚀 **HPC Ready**: Optimized for SLURM clusters with efficient resource management

## Pipeline Components

### Quality Control & Assembly
- **FastQC** - Raw read quality assessment
- **fastp** - Adapter trimming and quality filtering
- **Unicycler** - Phage genome assembly with circular genome detection
- **QUAST** - Assembly quality metrics

### Quality Assessment
- **CheckV** - Completeness estimation and contamination detection

### Annotation & Characterization
- **Pharokka** - Comprehensive phage annotation (CDS, tRNA, tmRNA, CRISPR)
- **PHANOTATE** - Phage-specific gene prediction
- **VIBRANT** - Lifestyle prediction (lytic vs lysogenic)
- **DIAMOND** - Prophage database comparison

### Reporting
- **MultiQC** - Aggregated quality control reports
- **Custom Reports** - Interactive HTML summaries with key metrics

## Quick Start

### Prerequisites
- Nextflow >= 24.04
- Apptainer/Singularity
- SLURM scheduler (optional)

### Installation

```bash
git clone https://github.com/tdoerks/PHINDER.git
cd PHINDER
```

### Basic Usage

#### Option 1: From Raw Reads (Recommended)
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --input_mode reads
```

**Samplesheet format (reads mode):**
```csv
sample,read1,read2
Phage1,/path/to/phage1_R1.fastq.gz,/path/to/phage1_R2.fastq.gz
Phage2,/path/to/phage2_R1.fastq.gz,/path/to/phage2_R2.fastq.gz
```

#### Option 2: From Assemblies
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --input_mode assembly
```

**Samplesheet format (assembly mode):**
```csv
sample,assembly
Phage1,/path/to/phage1.fasta
Phage2,/path/to/phage2.fasta
```

## Parameters

### Core Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Input samplesheet (CSV) | Required |
| `--input_mode` | Input type: `reads` or `assembly` | `reads` |
| `--outdir` | Output directory | `results` |

### Assembly Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--assembler` | Assembler to use: `unicycler` or `spades` | `unicycler` |
| `--skip_assembly` | Skip assembly (use with assembly mode) | `false` |

### Annotation Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_pharokka` | Skip Pharokka annotation | `false` |
| `--skip_vibrant` | Skip VIBRANT lifestyle prediction | `false` |

### Database Paths
| Parameter | Description | Setup Required |
|-----------|-------------|----------------|
| `--checkv_db` | CheckV database path | Yes |
| `--pharokka_db` | Pharokka database path | Auto-downloaded |
| `--prophage_db` | Prophage DIAMOND database (.dmnd) | Yes |

## Output Structure

```
results/
├── fastqc/                      # Raw read quality reports
├── fastp/                       # Trimming reports and trimmed reads
├── assemblies/                  # Assembled genomes
├── quast/                       # Assembly quality metrics
├── checkv/                      # Quality and completeness assessment
├── pharokka/                    # Comprehensive phage annotations
├── vibrant/                     # Lifestyle predictions
├── diamond_prophage/            # Prophage database comparisons
├── multiqc/                     # Aggregated QC report
└── summary/                     # Final integrated reports
    ├── phinder_summary.tsv
    └── phinder_report.html
```

## Tools & Versions

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Read quality assessment |
| fastp | 0.23.4 | Read trimming |
| Unicycler | 0.5.0 | Phage assembly |
| QUAST | 5.2.0 | Assembly metrics |
| CheckV | 1.0.2 | Quality assessment |
| Pharokka | 1.7.0 | Phage annotation |
| PHANOTATE | 1.6.7 | Gene prediction |
| VIBRANT | 4.0 | Lifestyle prediction |
| DIAMOND | 2.0 | Database search |
| MultiQC | 1.25.1 | Report aggregation |

## Database Setup

### CheckV Database
```bash
# Download CheckV database (~1.3 GB)
wget https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz
```

### Pharokka Database
Pharokka automatically downloads its database on first run, or you can pre-download:
```bash
pharokka.py --install-databases
```

### Prophage Database
See [docs/DATABASE_SETUP.md](docs/DATABASE_SETUP.md) for creating the prophage DIAMOND database.

## Citation

If you use PHINDER, please cite the tools it uses. See [CITATIONS.md](CITATIONS.md) for full citation list.

## License

MIT License - see [LICENSE](LICENSE) for details

## Contributing

Issues and pull requests welcome at: https://github.com/tdoerks/PHINDER

## Contact

**Author**: Tyler Doerksen (@tdoerks)
**Issues**: https://github.com/tdoerks/PHINDER/issues

---

**Version**: 1.0.0-dev
**Status**: Active Development
