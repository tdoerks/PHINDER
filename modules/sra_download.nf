process DOWNLOAD_SRA {
    tag "$srr_id"
    publishDir "${params.outdir}/fastq", mode: 'copy', pattern: '*_*.fastq.gz'
    container = 'quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0'
    errorStrategy = 'retry'
    maxRetries = 3

    input:
    val srr_id

    output:
    tuple val(srr_id), path("${srr_id}_1.fastq.gz"), path("${srr_id}_2.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    """
    # Download and convert to FASTQ with fasterq-dump
    fasterq-dump ${srr_id} \\
        --threads ${task.cpus} \\
        --split-files \\
        --skip-technical \\
        --progress

    # Compress the FASTQ files
    gzip ${srr_id}_*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-tools: \$(fasterq-dump --version 2>&1 | grep fasterq-dump | sed 's/.*: //g')
    END_VERSIONS
    """
}
