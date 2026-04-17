FROM python:3.11-slim

LABEL maintainer="nf-meta-gwas"
LABEL description="MAGMA gene-based analysis environment"

# Place magma_v1.10_static.zip in this directory before building.
# The static binary requires no additional shared libraries.
COPY magma_v1.10_static.zip /tmp/magma.zip

RUN apt-get update && \
    apt-get install -y --no-install-recommends unzip && \
    cd /tmp && unzip magma.zip && \
    mv magma /usr/local/bin/magma && \
    chmod +x /usr/local/bin/magma && \
    rm -f /tmp/magma.zip && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pandas numpy

ENV PYTHONUNBUFFERED=1
