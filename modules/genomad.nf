process GENOMAD {
    tag "$sample_id"
    publishDir "${params.outdir}/genomad/${sample_id}", mode: 'copy'
    container = 'staphb/genomad:1.12.0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_genomad/"), emit: results
    path "${sample_id}_genomad/${sample_id}_summary/${sample_id}_virus_summary.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    // Rename assembly to sample_id so geNomad output files use sample_id as prefix
    // geNomad names all output files after the input filename stem
    def db = params.genomad_db ?: '/genomad_db'
    """
    cp ${assembly} ${sample_id}.fasta

    genomad end-to-end \\
        ${sample_id}.fasta \\
        ${sample_id}_genomad \\
        ${db} \\
        --threads ${task.cpus} \\
        --splits 8 \\
        --enable-score-calibration

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(genomad --version 2>&1 | head -1)
    END_VERSIONS
    """
}
