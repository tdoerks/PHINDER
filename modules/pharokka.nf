process PHAROKKA {
    tag "$sample_id"
    publishDir "${params.outdir}/pharokka", mode: 'copy'
    container = 'docker://quay.io/biocontainers/pharokka:1.7.0--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_pharokka"), emit: results
    path "${sample_id}_pharokka/${sample_id}.gbk", emit: genbank
    path "${sample_id}_pharokka/${sample_id}.gff", emit: gff
    path "${sample_id}_pharokka/${sample_id}_cds_functions.tsv", emit: functions
    path "versions.yml", emit: versions

    script:
    def db_arg = params.pharokka_db ? "-d ${params.pharokka_db}" : ""
    """
    pharokka.py \\
        -i ${assembly} \\
        -o ${sample_id}_pharokka \\
        -t ${task.cpus} \\
        -p ${sample_id} \\
        ${db_arg} \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version | sed 's/pharokka //g')
    END_VERSIONS
    """
}
