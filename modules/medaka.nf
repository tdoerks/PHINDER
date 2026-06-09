process MEDAKA {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}_assembly.fasta"
    // VERIFY container tag resolves in registry before first run
    container = 'quay.io/biocontainers/medaka:1.11.3--py39h05d5c5e_2'

    input:
    tuple val(sample_id), path(draft), path(long_reads), val(platform)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    path "versions.yml", emit: versions

    script:
    // Nanopore consensus polishing. Model auto-selected by medaka unless
    // params.medaka_model is set (basecaller-specific, e.g. r1041_e82_400bps_sup_v4.2.0)
    def model_arg = params.medaka_model ? "-m ${params.medaka_model}" : ""
    """
    medaka_consensus \\
        -i ${long_reads} \\
        -d ${draft} \\
        -o ${sample_id}_medaka \\
        -t ${task.cpus} \\
        ${model_arg}

    if [ -f ${sample_id}_medaka/consensus.fasta ] && [ -s ${sample_id}_medaka/consensus.fasta ]; then
        cp ${sample_id}_medaka/consensus.fasta ${sample_id}_assembly.fasta
    else
        echo "ERROR: Medaka polishing failed for ${sample_id}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version | sed 's/medaka //g')
    END_VERSIONS
    """
}
