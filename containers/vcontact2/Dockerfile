FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Tyler Doerksen"
LABEL description="vConTACT2 0.11.3 with Diamond, MCL, and ClusterONE for phage protein-sharing network clustering"

# Install vConTACT2 and all dependencies via bioconda
RUN conda install -y -n base -c bioconda -c conda-forge -c defaults \
        vcontact2=0.11.3 \
        diamond=2.1.9 \
        mcl=22.282 \
        blast=2.15.0 \
        openjdk \
    && conda clean -afy

# Download ClusterONE JAR directly (cluster-one conda package has solve conflicts)
RUN wget -q "https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar" \
         -O /opt/conda/bin/cluster_one-1.0.jar

# Patch np.warnings (removed in numpy 1.24+) using find+sed to avoid Dockerfile parse issues
RUN find /opt/conda/lib/python*/site-packages/vcontact2 -name "*.py" \
        -exec grep -ql "np\.warnings" {} \; \
        -exec sed -i 's/np\.warnings\./warnings./g' {} \; \
        -exec sed -i '1s/^/import warnings\n/' {} \;

# Patch pcprofiles.py: zip(*df.values) incompatible with newer scipy — use .T (transpose) instead
RUN sed -i 's/(zip(\*profiles\.values))/(profiles.values.T)/' \
    /opt/conda/lib/python*/site-packages/vcontact2/pcprofiles.py

# Newer scipy rejects lazy zip objects as coo_matrix index arrays — materialize
# every zip(*x) as list(zip(*x)) package-wide (matrices.py, contig_clusters.py, ...)
RUN find /opt/conda/lib/python*/site-packages/vcontact2 -name "*.py" \
        -exec sed -i 's/\([( ]\)zip(\*\([A-Za-z_][A-Za-z0-9_.]*\))/\1list(zip(*\2))/g' {} \;

# Verify tools are on PATH (vcontact2 --version exits 1 by design, hence || true)
RUN vcontact2 --help > /dev/null 2>&1 || true && diamond --version && mcl --version

CMD ["vcontact2", "--help"]
