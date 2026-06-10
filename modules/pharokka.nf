process PHAROKKA {
    tag "$sample_id"
    publishDir "${params.outdir}/pharokka", mode: 'copy'
    // Use build _0, NOT _2: the _2 rebuild drifted to MMseqs2 v18 (pharokka 1.7.5
    // hard-requires v13.45111) and a numpy-2 stack whose subnormal UserWarning
    // polluted phanotate's version string -> check_dependencies crashed with
    // "ValueError: invalid literal for int()". The original _0 build has the
    // correct MMseqs2 v13.45111 + numpy-1.x (clean phanotate parse) and passes
    // the full dependency check end-to-end (verified on a Beocat node, 2026-06-10).
    // Same 1.7.5 version, so the installed pharokka_db is unchanged.
    container = 'docker://quay.io/biocontainers/pharokka:1.7.5--pyhdfd78af_0'

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
