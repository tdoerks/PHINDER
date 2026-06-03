process PHANOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/phanotate", mode: 'copy'
    container = 'docker://quay.io/biocontainers/phanotate:1.6.7--py311he264feb_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_phanotate.gff"), emit: gff
    tuple val(sample_id), path("${sample_id}_phanotate.faa"), emit: proteins
    path "versions.yml", emit: versions

    script:
    """
    # Run PHANOTATE for phage-specific gene prediction
    phanotate.py ${assembly} -o ${sample_id}_phanotate.gff -f gff3

    # Also generate protein sequences
    phanotate.py ${assembly} -o ${sample_id}_phanotate.faa -f fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phanotate: "1.6.7"
    END_VERSIONS
    """
}
