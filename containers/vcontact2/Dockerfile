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

# modules.py: groupby("module") yields float keys (NaN in the column) and float pos arrays;
# newer scipy rejects inexact sparse indices ("Inexact indices into sparse matrices are not
# allowed"). Hit only when enough related genomes reach module detection (PHINDER stress test).
RUN sed -i 's/^\(\s*\)N\[:, m\] = matrix\[:, pos\]\.sum(1)$/\1N[:, int(m)] = matrix[:, pos.astype(int)].sum(1)/' \
        /opt/conda/lib/python*/site-packages/vcontact2/modules.py \
    && grep -q "N\[:, int(m)\]" /opt/conda/lib/python*/site-packages/vcontact2/modules.py

# modules.py link_modules_and_clusters_df: the first merge keeps modules' "pos" column, so the
# second merge's suffixes would rename it to "pos_module", which already exists. pandas >= 2
# raises MergeError ("Passing 'suffixes' which cause duplicate columns {'pos_module'}");
# old pandas silently created a duplicate column. Drop the redundant "pos" (== pos_module).
# Patch the ONE file Python imports: /opt/conda/lib/python* can match several names that point at the
# same directory, which inserted this line three times in an earlier build (KeyError on the 2nd drop).
# The grep guard makes the insert idempotent; the final check requires exactly one copy.
RUN F="$(python -c 'import os, vcontact2; print(os.path.join(os.path.dirname(vcontact2.__file__), "modules.py"))')" \
    && (grep -q 'data = data.drop(columns="pos")' "$F" \
        || sed -i 's/^\(\s*\)data = pandas\.merge(data, modules, left_on="pos_module", right_on="pos",  how="left")$/&\n\1data = data.drop(columns="pos")/' "$F") \
    && [ "$(grep -c 'data = data.drop(columns="pos")' "$F")" -eq 1 ]

# Verify tools are on PATH (vcontact2 --version exits 1 by design, hence || true)
RUN vcontact2 --help > /dev/null 2>&1 || true && diamond --version && mcl --version

CMD ["vcontact2", "--help"]
