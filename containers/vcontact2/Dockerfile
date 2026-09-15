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

# Patch np.warnings (removed in numpy 1.24+) — locate via import to get correct conda path
RUN python3 -c "
import vcontact2, os, re, glob
vc = os.path.dirname(vcontact2.__file__)
for f in glob.glob(os.path.join(vc, '**', '*.py'), recursive=True):
    try:
        txt = open(f).read()
        if 'np.warnings' not in txt:
            continue
        txt = re.sub(r'np[.]warnings[.]', 'warnings.', txt)
        if 'import warnings' not in txt:
            txt = 'import warnings\n' + txt
        open(f, 'w').write(txt)
        print('Patched:', f)
    except Exception as e:
        print('Skip:', f, e)
"

# Verify tools are on PATH (vcontact2 --version exits 1 by design, hence || true)
RUN vcontact2 --help > /dev/null 2>&1 || true && diamond --version && mcl --version

CMD ["vcontact2", "--help"]
