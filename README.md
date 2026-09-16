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
- **SPAdes** - Phage genome assembly
- **QUAST** - Assembly quality metrics

### Quality Assessment
- **CheckV** - Completeness estimation and contamination detection

### Annotation & Characterization
- **Pharokka** - Comprehensive phage annotation (CDS, tRNA, tmRNA, CRISPR)
- **PHANOTATE** - Phage-specific gene prediction
- **VIBRANT** - Lifestyle prediction (lytic vs lysogenic)
- **BacPhlip** - ML-based lifestyle prediction (virulent vs temperate)
- **DIAMOND** - Prophage database comparison
- **geNomad** - Virus scoring and ICTV taxonomy assignment
- **PhageTerm** - Genome termini and packaging mechanism detection

### Safety Screening
- **AMRFinder Plus** - Antimicrobial resistance, virulence, and stress genes
- **Pharokka CARD/VFDB** - AMR and virulence factor screening

### Comparative Genomics
- **fastANI** - All-vs-all average nucleotide identity + species-level clustering (95% ANI)
- **vConTACT2** - Protein-sharing network clustering into viral clusters (custom container `ghcr.io/tdoerks/phinder-vcontact2`)
- **iPHoP** - Host genus prediction (optional; requires ~400 GB database)

### Reporting
- **MultiQC** - Aggregated quality control reports
- **Interactive dashboard** - Single-file HTML report: overview, lifestyle, annotation, quality, taxonomy, ANI heatmap, phylogenomics, AMR/VF screening, and full run provenance (parameters + software versions)

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

#### Option 1: From NCBI SRA Accessions (Easiest!)
```bash
nextflow run main.nf \
    --input test_lambda_sra.txt \
    --input_mode sra \
    --outdir results \
    -profile slurm
```

**SRA list format (one SRR per line):**
```
SRR5131134
SRR5131135
SRR5131136
```

#### Option 2: From Raw Reads
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

#### Option 3: From Assemblies
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
| `--input` | Input file (CSV samplesheet or TXT SRR list) | Required |
| `--input_mode` | Input type: `sra`, `reads`, or `assembly` | `reads` |
| `--outdir` | Output directory | `results` |

### Assembly Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--assembler` | Assembler to use | `spades` |
| `--skip_assembly` | Skip assembly (use with assembly mode) | `false` |

### Module Toggles
Every analysis module has a `--skip_<module>` flag (default `false` unless noted):
`skip_fastqc`, `skip_fastp`, `skip_pharokka`, `skip_vibrant`, `skip_phanotate`, `skip_diamond`, `skip_bacphlip`, `skip_checkv`, `skip_amrfinderplus`, `skip_genomad`, `skip_fastani`, `skip_phageterm`, `skip_vcontact2`, `skip_iphop` (default `true` — requires large database)

### Database Paths
| Parameter | Description | Setup Required |
|-----------|-------------|----------------|
| `--checkv_db` | CheckV database path | Yes |
| `--pharokka_db` | Pharokka database path | Yes — pre-install with `install_databases.py` (no runtime auto-download) |
| `--prophage_db` | Prophage DIAMOND database (.dmnd) | Yes |
| `--amrfinder_db` | AMRFinder database | No — bundled in container |
| `--genomad_db` | geNomad database | No — bundled in container |
| `--vcontact2_db` | vConTACT2 reference DB | No — defaults to clustering input phages only |
| `--iphop_db` | iPHoP database (~400 GB) | Only if `--skip_iphop false` |

Site-specific paths can be kept in a profile — see `conf/beocat.config` for an example, loaded with `-profile slurm,beocat`.

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
├── bacphlip/                    # ML lifestyle predictions
├── phanotate/                   # Phage gene predictions
├── diamond_prophage/            # Prophage database comparisons
├── amrfinderplus/               # AMR / virulence / stress genes
├── genomad/                     # Virus scores + ICTV taxonomy
├── phageterm/                   # Packaging mechanism / termini
├── fastani/                     # All-vs-all ANI matrix
├── vcontact2/                   # Viral cluster network results
├── multiqc/                     # Aggregated QC report
└── summary/                     # Final integrated reports
    ├── phinder_summary.tsv
    └── phinder_summary.html     # Interactive dashboard
```

## Tools & Versions

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Read quality assessment |
| fastp | 0.23.4 | Read trimming |
| SPAdes | 3.15.5 | Phage assembly |
| QUAST | 5.2.0 | Assembly metrics |
| CheckV | 1.0.2 | Quality assessment |
| Pharokka | 1.7.5 | Phage annotation |
| PHANOTATE | 1.6.7 | Gene prediction |
| VIBRANT | 4.0 | Lifestyle prediction |
| BacPhlip | 0.9.6 | ML lifestyle prediction |
| DIAMOND | 2.1.8 | Database search |
| AMRFinder Plus | 4.2.7 | AMR / virulence screening |
| geNomad | 1.12.0 | Virus scoring + taxonomy |
| PhageTerm | 1.0.12 | Packaging mechanism |
| fastANI | 1.34 | Average nucleotide identity |
| vConTACT2 | 0.11.3 | Viral cluster networks |
| iPHoP | 1.3.3 | Host prediction (optional) |
| MultiQC | 1.25.1 | Report aggregation |

Exact container URIs are pinned in each module under `modules/`; versions used in a given run are recorded in the dashboard's Run Info tab.

## Database Setup

### CheckV Database
```bash
# Download CheckV database (~1.3 GB)
wget https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz
```

### Pharokka Database
Pharokka does **not** auto-download its database at runtime — pre-install once:
```bash
install_databases.py -o /path/to/pharokka_db
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

**Version**: 1.0.0
**Status**: Stable
