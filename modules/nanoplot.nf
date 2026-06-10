process NANOPLOT {
    tag "$sample_id"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'
    // 1.42.0 shipped a numpy 1.x/2.x ABI break (pandas/pyarrow _ARRAY_API not found);
    // 1.47.0 has a coherent numpy stack. Verified `NanoPlot --version` runs clean
    // in the container on Beocat (2026-06-10). Same lesson as the Pharokka ABI fix.
    container = 'quay.io/biocontainers/nanoplot:1.47.0--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(long_reads)

    output:
    tuple val(sample_id), path("${sample_id}_nanoplot"), emit: results
    path "${sample_id}_nanoplot/NanoStats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    """
    NanoPlot \\
        --fastq ${long_reads} \\
        --outdir ${sample_id}_nanoplot \\
        --threads ${task.cpus} \\
        --tsv_stats \\
        --N50

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version | sed 's/NanoPlot //g')
    END_VERSIONS
    """
}
