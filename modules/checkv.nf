process CHECKV {
    tag "$sample_id"
    publishDir "${params.outdir}/checkv", mode: 'copy'
    container = 'docker://quay.io/biocontainers/checkv:1.0.2--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_checkv"), emit: results
    path "${sample_id}_checkv/quality_summary.tsv", emit: quality
    path "${sample_id}_checkv/completeness.tsv", emit: completeness, optional: true
    path "versions.yml", emit: versions

    script:
    """
    checkv end_to_end \\
        ${assembly} \\
        ${sample_id}_checkv \\
        -t ${task.cpus} \\
        -d ${params.checkv_db}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h | grep 'version' | sed 's/.*version //g')
    END_VERSIONS
    """
}
