FROM python:3.11-slim

LABEL maintainer="nf-meta-gwas"
LABEL description="gwaslab environment for GWAS summary statistics processing"

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    bcftools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir "gwaslab>=4.0"

# Verify installations
RUN python -c "from importlib.metadata import version; print('gwaslab', version('gwaslab'), 'installed successfully')"
RUN bcftools --version | head -1
RUN tabix --version | head -1

ENV PYTHONUNBUFFERED=1
