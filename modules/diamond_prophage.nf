process DIAMOND_PROPHAGE {
    tag "$sample_id"
    publishDir "${params.outdir}/diamond_prophage", mode: 'copy'
    container = 'docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_diamond_results.tsv"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    # Run DIAMOND blastx against prophage database
    diamond blastx \\
        --query ${assembly} \\
        --db ${params.prophage_db} \\
        --out ${sample_id}_diamond_results.tsv \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \\
        --evalue 1e-5 \\
        --max-target-seqs 10 \\
        --threads ${task.cpus} \\
        --sensitive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version | sed 's/diamond version //g')
    END_VERSIONS
    """
}
