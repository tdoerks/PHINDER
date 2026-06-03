process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{json,html}"
    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy', pattern: "*_trimmed*.fastq.gz"
    container = 'docker://staphb/fastp:0.23.4'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: reads
    path "${sample_id}_fastp.json", emit: json
    path "${sample_id}_fastp.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${sample_id}_trimmed_R1.fastq.gz \\
        --out2 ${sample_id}_trimmed_R2.fastq.gz \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} \\
        --qualified_quality_phred 20 \\
        --length_required 50 \\
        --detect_adapter_for_pe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //g')
    END_VERSIONS
    """
}
