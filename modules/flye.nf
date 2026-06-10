process FLYE {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}_assembly.fasta"
    publishDir "${params.outdir}/flye", mode: 'copy', pattern: "${sample_id}_flye.log"
    // Verified on quay.io/biocontainers (2026-06-10)
    container = 'quay.io/biocontainers/flye:2.9.5--py39hdf45acc_0'

    input:
    tuple val(sample_id), path(long_reads), val(platform)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    path "${sample_id}_flye.log", emit: log
    path "versions.yml", emit: versions

    script:
    // Read-type preset by platform. Modern defaults: --nano-hq (Guppy5+/R10),
    // --pacbio-hifi. Override via params.flye_read_type for older chemistries
    // (e.g. --nano-raw, --pacbio-raw).
    def preset = params.flye_read_type ?: (platform == 'pacbio' ? '--pacbio-hifi' : '--nano-hq')
    """
    flye \\
        ${preset} ${long_reads} \\
        --out-dir ${sample_id}_flye \\
        --threads ${task.cpus}

    if [ -f ${sample_id}_flye/assembly.fasta ] && [ -s ${sample_id}_flye/assembly.fasta ]; then
        cp ${sample_id}_flye/assembly.fasta ${sample_id}_assembly.fasta
        cp ${sample_id}_flye/flye.log ${sample_id}_flye.log
    else
        echo "ERROR: Flye assembly failed for ${sample_id}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """
}
