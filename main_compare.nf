#!/usr/bin/env nextflow

/*
========================================================================================
    PHINDER — SPAdes mode comparison entry point
    Tests all 6 SPAdes modes (isolate, meta, metaviral, rnaviral, careful, standard)
    on typed phage/bacterial samples and produces a comparison report.

    Usage:
        nextflow run main_compare.nf \
            -profile slurm,beocat \
            --compare_input samplesheets/samplesheet_spades_compare.csv \
            --outdir results_spades_compare

    See samplesheets/samplesheet_spades_compare_template.csv for samplesheet format.
========================================================================================
*/

nextflow.enable.dsl = 2

include { SPADES_MODE_COMPARE } from './workflows/spades_mode_compare'

params.compare_input = 'samplesheets/samplesheet_spades_compare.csv'
params.skip_checkv   = true

workflow {
    SPADES_MODE_COMPARE()
}
