process SPADES {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    container = 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    path "${sample_id}_spades.log", emit: log
    path "versions.yml", emit: versions

    script:
    """
    spades.py \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${sample_id}_spades \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        --isolate \\
        --only-assembler

    # Use contigs.fasta (best for isolate mode)
    if [ -f ${sample_id}_spades/contigs.fasta ] && [ -s ${sample_id}_spades/contigs.fasta ]; then
        cp ${sample_id}_spades/contigs.fasta ${sample_id}_assembly.fasta
        cp ${sample_id}_spades/spades.log ${sample_id}_spades.log
    else
        echo "ERROR: SPAdes assembly failed for ${sample_id}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: "3.15.5"
    END_VERSIONS
    """
}
