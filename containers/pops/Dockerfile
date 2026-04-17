FROM python:3.8-slim

LABEL maintainer="nf-meta-gwas"
LABEL description="PoPS (Polygenic Priority Score) environment, pinned to FinucaneLab/pops@76eb86c"

RUN apt-get update && \
    apt-get install -y --no-install-recommends git && \
    git clone https://github.com/FinucaneLab/pops.git /opt/pops && \
    cd /opt/pops && git checkout 76eb86c && \
    pip install --no-cache-dir \
        scipy==1.5.2 \
        pandas==1.0.5 \
        numpy==1.19.5 \
        scikit-learn==0.24.2 && \
    sed -i '1s|^|#!/usr/bin/env python3\n|' /opt/pops/pops.py /opt/pops/munge_feature_directory.py && \
    chmod +x /opt/pops/pops.py /opt/pops/munge_feature_directory.py && \
    rm -rf /opt/pops/.git && \
    apt-get remove -y git && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/pops:$PATH"
ENV PYTHONUNBUFFERED=1
