# PHINDER Database Setup

PHINDER requires several databases for comprehensive phage analysis. This guide covers downloading and setting up each database.

## Required Databases

### 1. CheckV Database (~1.3 GB)

CheckV requires a database for phage quality assessment and completeness estimation.

**Download and setup:**
```bash
# Create databases directory
mkdir -p ~/databases
cd ~/databases

# Download CheckV database
wget https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz

# Verify
ls -lh checkv-db-v1.5/
```

**REQUIRED: rebuild the DIAMOND database with the container's own DIAMOND.**

The `checkv_reps.dmnd` shipped in the CheckV download is built with a *newer*
DIAMOND than the one bundled in the `checkv:1.0.2` container (db build 167 vs
container build 162). `diamond dbinfo` reads it fine (header only), but the
completeness step's `diamond blastp` fails with the misleading message
`DIAMOND task failed. Program should be rerun.` — while contamination (HMMER)
succeeds. Rebuild the `.dmnd` from the bundled `checkv_reps.faa` using the
container's DIAMOND so the database build matches the binary:

```bash
cd checkv-db-v1.5/genome_db
mv checkv_reps.dmnd checkv_reps.dmnd.orig   # keep the original
apptainer exec docker://quay.io/biocontainers/checkv:1.0.2--pyhdfd78af_0 \
    diamond makedb --in checkv_reps.faa --db checkv_reps
```

This is a one-time step per database install. Confirm it worked:
```bash
apptainer exec docker://quay.io/biocontainers/checkv:1.0.2--pyhdfd78af_0 \
    diamond blastp --query <any_proteins.faa> \
    --db checkv-db-v1.5/genome_db/checkv_reps.dmnd \
    --out /tmp/test.tsv --outfmt 6 --threads 4 --tmpdir /tmp
# Exit 0 + alignments reported = good. See bin/diagnose_checkv_diamond.sh.
```

**Update `nextflow.config`:**
```groovy
params {
    checkv_db = '/path/to/checkv-db-v1.5'
}
```

### 2. Pharokka Database (REQUIRED — manual, one-time)

> **Important:** Pharokka does **not** download its database at runtime. The
> `pharokka.py` annotation run needs the PHROGs database to already exist, and
> the BioContainers image does **not** bundle it. You must install it once with
> `install_databases.py` and point `pharokka_db` at it, or every run fails at
> the PHAROKKA step.

**Install once (via the pipeline's own container, no conda needed):**
```bash
apptainer exec docker://quay.io/biocontainers/pharokka:1.7.5--pyhdfd78af_2 \
    install_databases.py -o /fastscratch/tylerdoe/databases/pharokka_db
```

> Note: use build `1.7.5--pyhdfd78af_2`. The earlier `1.7.0--pyhdfd78af_0`
> build ships an incompatible numpy/pandas pairing (`numpy.dtype size changed`
> ABI error) that breaks both `install_databases.py` and `pharokka.py`.

**Set in `nextflow.config`** (already defaulted to this path):
```groovy
params {
    pharokka_db = '/fastscratch/tylerdoe/databases/pharokka_db'
}
```

### 3. Prophage DIAMOND Database (~500 MB)

DIAMOND prophage database for taxonomic classification.

**Option A: Use existing database**
If you have a prophage protein FASTA file:

```bash
# Create DIAMOND database from FASTA
apptainer exec docker://staphb/diamond diamond makedb \
    --in prophage_proteins.faa \
    --db prophage_db

# This creates prophage_db.dmnd
```

**Option B: Download from custom source**
Contact your lab for the prophage database or use a public database.

**Update `nextflow.config`:**
```groovy
params {
    prophage_db = '/path/to/prophage_db.dmnd'
}
```

## Optional Databases

### VIBRANT Databases

VIBRANT downloads its own databases automatically on first run. No manual setup required.

### PHANOTATE

PHANOTATE has built-in models. No external database needed.

## Verify Setup

After setting up databases, verify paths:

```bash
# Check CheckV
ls -lh ~/databases/checkv-db-v1.5/

# Check Prophage DB
ls -lh ~/databases/prophage_db.dmnd

# Check Pharokka (if pre-downloaded)
ls -lh ~/databases/pharokka_db/
```

## Storage Requirements

| Database | Size | Auto-download |
|----------|------|---------------|
| CheckV | ~1.3 GB | No - manual |
| Pharokka | ~800 MB | No - manual (install_databases.py) |
| Prophage DIAMOND | ~500 MB | No - manual |
| VIBRANT | ~800 MB | Yes (first run) |

**Total:** ~3.4 GB

## Troubleshooting

### CheckV errors

**Error:** `Database not found`
- Verify the path in `nextflow.config` matches your actual database location
- Ensure database was fully extracted (check for subdirectories)

### Pharokka database issues

**Error:** `Cannot download database`
- Try manual download with `pharokka.py --install-databases`
- Check internet connectivity
- Verify sufficient disk space

### DIAMOND database errors

**Error:** `Error opening database file`
- Ensure `.dmnd` extension is present
- Verify file permissions: `chmod 644 prophage_db.dmnd`
- Re-create database if corrupted

## HPC-Specific Notes

For HPC environments (e.g., Beocat), store databases in a shared location:

```bash
# Recommended structure
/homes/$USER/databases/
├── checkv-db-v1.5/
├── pharokka_db/
└── prophage_db.dmnd
```

Update paths in your configuration profile accordingly.

## Database Updates

Databases should be updated periodically:

- **CheckV**: Check https://portal.nersc.gov/CheckV/ for updates
- **Pharokka**: Run `pharokka.py --update-databases`
- **Prophage DB**: Re-download source FASTA and rebuild DIAMOND database

## Contact

For database access issues, contact:
- Tyler Doerksen (@tdoerks)
- Open an issue: https://github.com/tdoerks/PHINDER/issues
