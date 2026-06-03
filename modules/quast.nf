process QUAST {
    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'
    container = 'docker://staphb/quast:5.2.0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_quast"), emit: results
    path "${sample_id}_quast/report.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    quast.py \\
        ${assembly} \\
        -o ${sample_id}_quast \\
        -t ${task.cpus} \\
        --min-contig 500

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version | sed 's/QUAST v//g')
    END_VERSIONS
    """
}
