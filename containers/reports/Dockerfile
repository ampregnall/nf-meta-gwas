FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="nf-meta-gwas"
LABEL description="Quarto reporting environment for GWAS meta-analysis HTML reports"

COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y --no-install-recommends curl && \
    curl -fsSL https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.57/quarto-1.5.57-linux-amd64.deb \
        -o /tmp/quarto.deb && \
    dpkg -i /tmp/quarto.deb && \
    rm /tmp/quarto.deb && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV JUPYTER_PLATFORM_DIRS=1
ENV PATH="/opt/conda/bin:$PATH"
