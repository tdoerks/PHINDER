process VCONTACT2 {
    publishDir "${params.outdir}/vcontact2", mode: 'copy'
    container = 'quay.io/biocontainers/vcontact2:0.11.3--pyhdfd78af_1'

    input:
    val(sample_ids)
    path(faa_files)

    output:
    path "vcontact2_results/", emit: results
    path "vcontact2_results/genome_by_genome_overview.csv", emit: clusters, optional: true
    path "versions.yml", emit: versions

    script:
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]
    def fns = (faa_files instanceof List ? faa_files : [faa_files]).collect { it.name }
    // Rename each .faa to sample_id.faa so protein-to-genome mapping is unambiguous
    def renames = [ids, fns].transpose().collect { id, fn ->
        fn != "${id}.faa" ? "cp '${fn}' '${id}.faa'" : "true"
    }.join(' && ')
    def db_arg = params.vcontact2_db ? "--db '${params.vcontact2_db}'" : "--db 'None'"
    """
    ${renames}

    # Build gene-to-genome CSV: every protein header maps to its sample_id
    python3 -c "
import os, sys
print('protein_id,contig_id,keywords')
for fname in sorted(f for f in os.listdir('.') if f.endswith('.faa')):
    genome_id = fname[:-4]
    with open(fname) as fh:
        for line in fh:
            if line.startswith('>'):
                prot_id = line[1:].split()[0]
                print(f'{prot_id},{genome_id},')
" > gene_to_genome.csv

    cat *.faa > all_proteins.faa

    vcontact2 \\
        --raw-proteins all_proteins.faa \\
        --rel-mode Diamond \\
        --proteins-fp gene_to_genome.csv \\
        ${db_arg} \\
        --pcs-mode MCL \\
        --vcs-mode ClusterONE \\
        --output-dir vcontact2_results \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact2: \$(vcontact2 --version 2>&1 | head -1)
    END_VERSIONS
    """
}
