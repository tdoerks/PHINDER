process VCONTACT2 {
    publishDir "${params.outdir}/vcontact2", mode: 'copy'
    container = 'ghcr.io/tdoerks/phinder-vcontact2:0.11.3'

    input:
    val(sample_ids)
    path(faa_files)

    output:
    path "vcontact2_results/", emit: results
    path "vcontact2_results/genome_by_genome_overview.csv", emit: clusters, optional: true
    path "versions.yml", emit: versions

    script:
    def db_arg = params.vcontact2_db ? "--db '${params.vcontact2_db}'" : "--db 'None'"
    // Staged .faa files are already named <sample_id>.faa by the Pharokka module.
    // Do NOT rename/pair against sample_ids here: the two collected channels
    // arrive in different orders, and pairwise cp swaps clobber file contents.
    """
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

    set +e
    vcontact2 \\
        --raw-proteins all_proteins.faa \\
        --rel-mode Diamond \\
        --proteins-fp gene_to_genome.csv \\
        ${db_arg} \\
        --pcs-mode MCL \\
        --vcs-mode ClusterONE \\
        --c1-bin /opt/conda/bin/cluster_one-1.0.jar \\
        --output-dir vcontact2_results \\
        --threads ${task.cpus} 2> vcontact2.stderr
    status=\$?
    set -e
    cat vcontact2.stderr >&2

    if [ \$status -ne 0 ]; then
        if grep -q "No edge in the similarity network" vcontact2.stderr; then
            # Legitimate outcome for small/dissimilar genome sets with --db None:
            # no genome pair shares enough protein clusters to form an edge.
            echo "WARNING: vConTACT2 found no edges — genomes too dissimilar to cluster" >&2
            mkdir -p vcontact2_results
        else
            exit \$status
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact2: \$(vcontact2 --version 2>&1 | head -1)
    END_VERSIONS
    """
}
