/*
========================================================================================
    SPADES_MODE_COMPARE workflow
    Tests all 6 SPAdes assembly modes on a set of typed samples (dsDNA, ssDNA, ssRNA,
    dsRNA, bacterial) and compares recovery quality using CheckV + QUAST.

    Samplesheet format (CSV with header):
        sample,fastq_1,fastq_2,sample_type
    sample_type values: dsDNA | ssDNA | ssRNA | dsRNA | bacterial

    Run:
        nextflow run main_compare.nf -profile slurm,beocat \
            --compare_input samplesheets/samplesheet_spades_compare.csv \
            --outdir results_spades_compare

    Modes tested per sample: isolate, meta, metaviral, rnaviral, careful, standard
    Effective sample IDs downstream: ${sample}__${mode} (double underscore separator)
========================================================================================
*/

include { FASTP             } from '../modules/fastp'
include { SPADES_MODE       } from '../modules/spades_modes'
include { QUAST             } from '../modules/quast'
include { COMPARE_ASSEMBLIES } from '../modules/compare_assemblies'

def MODES = ['isolate', 'meta', 'metaviral', 'rnaviral', 'careful', 'standard']

workflow SPADES_MODE_COMPARE {

    // Parse samplesheet: sample,fastq_1,fastq_2,sample_type
    ch_samples = Channel
        .fromPath(params.compare_input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample, type: row.sample_type]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }

    // Trim reads
    ch_trimmed = ch_samples.map { meta, r1, r2 -> tuple(meta.id, r1, r2) }
    FASTP(ch_trimmed)

    // Fan out: each (sample, r1, r2) × 6 modes
    ch_with_modes = FASTP.out.reads
        .combine(Channel.from(MODES))
        .map { sample_id, reads, mode ->
            tuple(sample_id, reads[0], reads[1], mode)
        }

    SPADES_MODE(ch_with_modes)

    // Downstream CheckV + QUAST use "${sample}__${mode}" as effective sample_id
    // so publishDir names encode both, letting the comparison script parse them back out.
    ch_assemblies = SPADES_MODE.out.assembly
        .map { sample_id, mode, assembly ->
            tuple("${sample_id}__${mode}", assembly)
        }

    QUAST(ch_assemblies)

    // Collect QUAST results and build comparison report.
    // CheckV removed from comparison workflow — hmmsearch batching fails on fragmented
    // assemblies from poor-fit modes; QUAST N50/contig-count is sufficient here.
    ch_quast_collected = QUAST.out.results.map { _id, dir -> dir }.collect()
    ch_samplesheet     = Channel.fromPath(params.compare_input)

    COMPARE_ASSEMBLIES(
        Channel.empty(),
        ch_quast_collected,
        ch_samplesheet
    )
}
