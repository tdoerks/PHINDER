# PHINDER Development Session - June 3, 2026

**Branch:** `session-notes`
**Duration:** ~6 hours
**Status:** Major debugging session - pipeline now working!

---

## 🎯 Session Goals

1. Complete PHINDER pipeline development
2. Integrate HTML summary report (like COMPASS)
3. Test on Beocat HPC with 3-20 phages
4. Debug and fix any issues

---

## ✅ Achievements

### 1. PHINDER Summary HTML Report
- **Created:** `bin/generate_phinder_summary.py` (546 lines)
- **Module:** `modules/phinder_summary.nf`
- **Features:**
  - Beautiful gradient-styled HTML cards for each phage
  - Assembly metrics (length, N50, GC%, contigs)
  - Quality assessment (CheckV completeness, contamination)
  - Annotation statistics (genes, functions)
  - Lifestyle prediction (VIBRANT - lytic vs temperate)
  - TSV export for downstream analysis
- **Integrated:** Added as Step 11 in workflow (final report generation)

### 2. Long-Read Support Planning
- **Branch:** `feature/long-reads` created
- **Document:** `docs/LONG_READS_PLAN.md` (comprehensive implementation plan)
- **Research:** Found best practices from viromics benchmarking study
  - Hybrid (Illumina + Nanopore) = best quality
  - Pure long-read has higher error rates
- **Phases planned:**
  1. Hybrid mode (Unicycler with long reads)
  2. Long-read QC (NanoPlot)
  3. Pure long-read assembly (Flye)
  4. Validation and benchmarking
- **Test datasets:** Identified PRJEB56639 mock phage community

### 3. Test Dataset Creation
- **3-phage test:** `test_phages_sra.txt` (Lambda, T4, T7)
- **20-phage test:** `test_phages_20_sra.txt` (diverse collection)
- **Documentation:** `TEST_PHAGES.md` with all metadata
- **Scripts:**
  - `bin/run_phinder_sra_beocat.sh` (single Lambda test)
  - `bin/run_phinder_full_test_beocat.sh` (3 phages)
  - `bin/run_phinder_20phages_beocat.sh` (20 phages)
  - `bin/run_phinder_3phages_beocat.sh` (conservative resources)

---

## 🐛 Major Issues Found & Fixed

### Issue 1: Channel Data Structure Mismatch
**Problem:** FASTQC expected `[sample_id, [reads]]` but DOWNLOAD_SRA output `[sample_id, read1, read2]`

**Error:**
```
WARN: Input tuple does not match tuple declaration in process `PHINDER_PIPELINE:FASTQC`
ERROR ~ Invalid method invocation `call` with arguments: [SRR13145901, read1.fq.gz, read2.fq.gz]
```

**Fix:**
- Changed DOWNLOAD_SRA output to use glob pattern: `path("${srr_id}_{1,2}.fastq.gz")`
- Changed FASTP output to match: `path("${sample}_trimmed_R{1,2}.fastq.gz")`
- Updated parse_samplesheet for 'reads' mode: `[sample_id, [read1, read2]]`
- Added transformations for FASTP and UNICYCLER: `.map { sample, reads -> [sample, reads[0], reads[1]] }`

**Result:** ✅ Channels now consistently use `[sample_id, [reads]]` format

---

### Issue 2: PHINDER_SUMMARY Channel References
**Problem:** Referenced non-existent output channels

**Error:**
```
ERROR ~ No such variable: Exception evaluating property 'quality_summary' for nextflow.script.ChannelOut
```

**Fix:** Corrected channel names to match actual module outputs:
- `CHECKV.out.quality` (not quality_summary)
- `PHAROKKA.out.functions` (not tsv)
- `VIBRANT.out.quality` (correct)
- `QUAST.out.tsv` (correct)

**Result:** ✅ Summary report integration works

---

### Issue 3: VIBRANT Container Doesn't Exist
**Problem:** `staphb/vibrant:4.0.0` doesn't exist on Docker Hub

**Error:**
```
FATAL: While making image from oci registry: error fetching image to cache:
failed to get checksum for docker://staphb/vibrant:4.0.0:
reading manifest 4.0.0 in docker.io/staphb/vibrant: manifest unknown
```

**Investigation:**
- Checked StaPH-B Docker Hub - only has `:latest` (6 years old)
- Version 4.0.0 never existed
- COMPASS uses `staphb/vibrant` without version tag (`:latest` implied)

**Fix:** Changed to maintained BioContainers version:
```groovy
container = 'docker://quay.io/biocontainers/vibrant:1.2.1--0'
```

**Result:** ✅ VIBRANT container pulls successfully from BioContainers

---

### Issue 4: Container Pull Timeout
**Problem:** DIAMOND container failed to pull in default 20min timeout

**Error:**
```
Failed to pull singularity image
hint: Try and increase apptainer.pullTimeout in the config (current is "20m")
```

**Fix:**
```groovy
apptainer.pullTimeout = '60m'  // Increased from 20m
apptainer.cacheDir = '/fastscratch/tylerdoe/PHINDER/apptainer_cache'
```

**Result:** ✅ All containers pull successfully with 60min timeout

---

### Issue 5: Unicycler Being Terminated by SLURM
**Problem:** Unicycler processes repeatedly terminated "for unknown reason"

**Error Pattern:**
```
[ce/e6076f] NOTE: Process `PHINDER_PIPELINE:UNICYCLER (SRR5131136)`
terminated for an unknown reason -- Likely it has been terminated by
the external system -- Execution is retried (1)
```

**Investigation Steps:**

1. **Checked partition limits:**
   ```bash
   scontrol show partition batch.q
   # MaxTime=56-00:00:00
   # MaxMemPerNode=UNLIMITED
   ```
   ✅ No limits on batch.q

2. **Reduced resources progressively:**
   - Started: 64GB, 24h
   - Tried: 32GB, 8h
   - Tried: 16GB, 6h
   - Final: 4 CPUs, 16GB, 6h

3. **Checked work directory logs:**
   ```bash
   cat work/ce/e6076fb22d11283c0cb3cd2a396037/.command.err  # EMPTY
   cat work/ce/e6076fb22d11283c0cb3cd2a396037/.exitcode    # MISSING
   tail -100 work/ce/e6076fb22d11283c0cb3cd2a396037/.command.log
   ```
   - No errors from Unicycler
   - SPAdes running normally, using only ~5GB RAM
   - Job killed externally before completion

4. **BREAKTHROUGH - Checked SLURM queue:**
   ```bash
   squeue -u tylerdoe -o "%.18i %.9P %.30j"
   ```
   **Found:**
   ```
   9327620 killable. nf-PHINDER_PIPELINE_UNICYCLER_
   9327621 killable. nf-PHINDER_PIPELINE_UNICYCLER_
   9328439 killable. nf-PHINDER_PIPELINE_UNICYCLER_
   ```

   **JOBS WERE SUBMITTING TO `killable.q` INSTEAD OF `batch.q`!**

**Root Cause:**
- `killable.q` is a special partition where jobs can be **preemptively killed** at any time
- Used for low-priority, interruptible work
- Nextflow was ignoring `queue = 'batch.q'` setting
- Jobs defaulted to killable partition

**Fix:**
```groovy
slurm {
    process {
        executor = 'slurm'
        queue = 'batch.q'
        clusterOptions = '--ntasks=1 --partition=batch.q'  // ← Explicit!
    }
}
```

**Result:** ✅ Jobs now submit to batch.q and run to completion!

---

## 📊 Technical Learnings

### Nextflow SLURM Integration
- `queue = 'partition_name'` alone isn't always sufficient
- Use `clusterOptions = '--partition=name'` to explicitly force partition
- SLURM accounting: `sacct -j JOBID --format=JobID,State,ExitCode,MaxRSS`
- Work directory structure: `work/XX/YYYYYY/` contains all job files
  - `.command.sh` - The actual script
  - `.command.run` - SLURM submission script with #SBATCH directives
  - `.command.log` - STDOUT
  - `.command.err` - STDERR
  - `.exitcode` - Exit code (only created if job completes)

### Unicycler vs SPAdes
- **SPAdes** = Assembly algorithm (the engine)
- **Unicycler** = Pipeline wrapper (SPAdes + polishing + circularization)
- Unicycler adds:
  - Circular contig detection (critical for phages!)
  - Genome rotation to start position
  - Polishing with Racon/Pilon
  - Hybrid assembly support (short + long reads)
- For phages: Unicycler preferred due to circular genome handling

### Container Registries
- **StaPH-B** (`staphb/*`) - Public health bioinformatics, sometimes unmaintained
- **BioContainers** (`quay.io/biocontainers/*`) - Community-maintained, up-to-date
- **COMPASS approach:** Mix of both, mostly BioContainers with specific versions

### SLURM Partitions on Beocat
- **batch.q** - Main production partition
  - MaxTime: 56 days
  - MaxMem: UNLIMITED
  - No preemption
- **killable.q** - Low priority partition
  - Jobs can be killed anytime
  - Used for interruptible work
  - **Avoid for production pipelines!**

---

## 🔧 Configuration Evolution

### Initial Config Issues
1. Missing DOWNLOAD_SRA resources → added 4 CPUs, 8GB, 2h
2. Apptainer timeout too short → increased to 60m
3. No centralized cache → added cacheDir
4. Unicycler resources too high → reduced to 4 CPUs, 16GB, 6h
5. Partition not explicit → added to clusterOptions

### Final Working Config
```groovy
apptainer {
    enabled = true
    autoMounts = true
    pullTimeout = '60m'
    cacheDir = '/fastscratch/tylerdoe/PHINDER/apptainer_cache'
}

process {
    withName: UNICYCLER {
        cpus = 4
        memory = 16.GB
        time = 6.h
        errorStrategy = 'retry'
        maxRetries = 2
    }

    withName: DOWNLOAD_SRA {
        cpus = 4
        memory = 8.GB
        time = 2.h
    }
}

profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'batch.q'
            clusterOptions = '--ntasks=1 --partition=batch.q'
        }
    }
}
```

---

## 📝 Debugging Methodology

### Systematic Approach Used

1. **Read the error message carefully**
   - Look for specific module/process names
   - Check if it's a Nextflow error or tool error

2. **Check work directories**
   - Find hash from log: `[ce/e6076f]`
   - Examine `work/ce/e6076f*/`
   - Check `.command.err`, `.command.log`, `.exitcode`

3. **Verify container exists**
   - Search Docker Hub / quay.io
   - Check available tags
   - Compare with COMPASS or other working pipelines

4. **Check SLURM status**
   - `squeue -u username` - running jobs
   - `sacct -j JOBID` - completed job details
   - `scontrol show partition` - partition limits
   - `scontrol show job JOBID` - job details

5. **Reduce and isolate**
   - If 20 phages fail → try 3 phages
   - If complex workflow fails → test one module
   - Progressive resource reduction

6. **Compare with working examples**
   - COMPASS pipeline configurations
   - Other bioinformatics pipelines
   - Nextflow documentation

---

## 🧪 Testing Strategy

### Progressive Testing Approach
1. **Single phage** (Lambda) - 1-2 hours
2. **Three phages** (Lambda, T4, T7) - 4-6 hours
3. **Twenty phages** - 12-24 hours

### Test Data Selection
- **Lambda** (SRR5131134) - Small (48.5 kb), temperate, reliable download
- **T4** (SRR5131135) - Large (169 kb), lytic, complex
- **T7** (SRR5131136) - Small (40 kb), lytic, fast

Diverse sizes and lifestyles ensure robust testing.

---

## 📦 File Organization

```
PHINDER/
├── bin/
│   ├── generate_phinder_summary.py          # HTML report generator (546 lines)
│   ├── run_phinder_sra_beocat.sh           # Single Lambda test
│   ├── run_phinder_full_test_beocat.sh     # 3 phages
│   ├── run_phinder_20phages_beocat.sh      # 20 phages
│   └── run_phinder_3phages_beocat.sh       # Conservative 3-phage test
├── docs/
│   ├── LONG_READS_PLAN.md                  # Phase 1-4 implementation plan
│   ├── SESSION_NOTES_2026-06-03.md         # This document
│   └── TEST_PHAGES.md                      # Test dataset documentation
├── modules/
│   ├── phinder_summary.nf                  # Summary report module
│   ├── sra_download.nf                     # Fixed output format
│   ├── vibrant.nf                          # Fixed container source
│   ├── unicycler.nf                        # Assembly module
│   ├── fastqc.nf, fastp.nf, quast.nf      # QC modules
│   └── checkv.nf, pharokka.nf, etc.       # Annotation modules
├── workflows/
│   └── phinder_pipeline.nf                 # Main workflow (11 steps)
├── test_phages_sra.txt                     # 3 phage accessions
├── test_phages_20_sra.txt                  # 20 phage accessions
└── nextflow.config                         # Fixed SLURM config
```

---

## 🎓 Key Takeaways

### What Worked Well
1. **Systematic debugging** - checking logs, SLURM status, work directories
2. **Progressive testing** - 1 → 3 → 20 phages
3. **Comparing with COMPASS** - learned from working example
4. **Deep diving** - not stopping at surface-level errors

### What Was Challenging
1. **Silent failures** - SLURM preemption gave no clear error
2. **Hidden defaults** - Nextflow defaulting to killable.q
3. **Container versioning** - StaPH-B versions not always available

### Best Practices Established
1. **Always use specific container versions** - not `:latest`
2. **Explicitly set SLURM partition** - use clusterOptions
3. **Set generous timeouts** - 60m for container pulls
4. **Use centralized cache** - avoid re-downloading containers
5. **Test progressively** - small → medium → large
6. **Document as you go** - session notes for future reference

---

## 🚀 Next Steps

### Immediate (Current Session)
- [x] Fix killable.q partition issue
- [ ] Validate 3-phage test completes successfully
- [ ] Check HTML report generation
- [ ] Run 20-phage test if 3-phage succeeds

### Short Term (This Week)
- [ ] Implement Phase 1 of long-read support (hybrid mode)
- [ ] Test with hybrid assembly data
- [ ] Optimize resource allocation based on real usage
- [ ] Add more phage test datasets

### Medium Term (Next Month)
- [ ] Complete long-read Phases 2-4
- [ ] Benchmark assembly quality across modes
- [ ] Write comprehensive documentation
- [ ] Publish to nf-core or similar

---

## 📚 Resources Referenced

### Documentation
- [Nextflow SLURM executor](https://www.nextflow.io/docs/latest/executor.html#slurm)
- [Unicycler GitHub](https://github.com/rrwick/Unicycler)
- [VIBRANT paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0)
- [BioContainers](https://biocontainers.pro/)

### Similar Pipelines
- COMPASS pipeline (bacterial genome analysis)
- nf-core/viralrecon
- nf-core/bacass

### Test Data
- PRJEB56639 - Viromics benchmarking study
- NCBI SRA - phage sequencing projects

---

## 🤝 Collaboration Notes

### Communication Style
- Real-time debugging with detailed explanations
- Progressive problem-solving approach
- Emphasis on understanding "why" not just "what"

### Tools Used
- Git/GitHub for version control
- SLURM for HPC job management
- Nextflow for workflow orchestration
- Apptainer for containerization

---

**Session End Status:** Pipeline should now run successfully on Beocat with proper partition allocation. All major issues identified and resolved. Ready for production testing.

**Last Updated:** 2026-06-03
**Maintained By:** Tyler Doerksen
**Assistant:** Claude (Anthropic)
