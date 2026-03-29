FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="nf-meta-gwas"
LABEL description="gwaslab environment for GWAS summary statistics processing"

COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV NUMBA_CACHE_DIR=/tmp/numba_cache
ENV PATH="/opt/conda/bin:$PATH"

RUN python -c "import gwaslab as gl; gl.download_ref('ensembl_hg38_gtf')"