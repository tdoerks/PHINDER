/*
========================================================================================
    COMPARE_ASSEMBLIES — aggregate CheckV + QUAST results across all SPAdes modes
    Produces a per-sample × per-mode comparison TSV and HTML report.
========================================================================================
*/

process COMPARE_ASSEMBLIES {
    publishDir "${params.outdir}/comparison", mode: 'copy'
    container = 'quay.io/biocontainers/pandas:1.5.2'

    input:
    path(checkv_dirs, stageAs: 'checkv_results/*')    // all *_checkv dirs
    path(quast_dirs,  stageAs: 'quast_results/*')     // all *_quast dirs
    path(samplesheet)                                  // original CSV with sample_type col

    output:
    path "spades_mode_comparison.tsv",  emit: tsv
    path "spades_mode_comparison.html", emit: html

    script:
    """
    python3 ${projectDir}/bin/compare_spades_modes.py \\
        --checkv-dir checkv_results \\
        --quast-dir  quast_results \\
        --samplesheet ${samplesheet} \\
        --outdir .
    """
}
