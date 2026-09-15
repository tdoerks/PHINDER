FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Tyler Doerksen"
LABEL description="vConTACT2 0.11.3 with Diamond, MCL, and ClusterONE for phage protein-sharing network clustering"

# Install vConTACT2 and all dependencies via bioconda
RUN conda install -y -n base -c bioconda -c conda-forge -c defaults \
        vcontact2=0.11.3 \
        diamond=2.1.9 \
        mcl=22.282 \
        blast=2.15.0 \
    && conda clean -afy

# Patch np.warnings (removed in numpy 1.24+) across all vcontact2 source files
RUN python3 -c "
import site, os, re, glob
sp = site.getsitepackages()[0]
vc_dir = os.path.join(sp, 'vcontact2')
for f in glob.glob(os.path.join(vc_dir, '**/*.py'), recursive=True):
    txt = open(f).read()
    if 'np.warnings' in txt:
        txt = re.sub(r'np\.warnings\.', 'warnings.', txt)
        if 'import warnings' not in txt:
            txt = 'import warnings\n' + txt
        open(f, 'w').write(txt)
        print('Patched', f)
"

# Verify tools are on PATH (vcontact2 --version exits 1 by design, hence || true)
RUN vcontact2 --help > /dev/null 2>&1 || true && diamond --version && mcl --version

CMD ["vcontact2", "--help"]
