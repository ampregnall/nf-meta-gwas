FROM python:3.8-slim

LABEL maintainer="nf-meta-gwas"
LABEL description="FLAMES gene prioritization environment, pinned to ampregnall/FLAMES@2008304"

RUN apt-get update && \
    apt-get install -y --no-install-recommends git && \
    git clone https://github.com/ampregnall/FLAMES.git /opt/flames && \
    cd /opt/flames && git checkout 7fc943c && \
    pip install --no-cache-dir \
        "Brotli==1.1.0" \
        "numpy==1.24.4" \
        "pandas==2.0.3" \
        "fastparquet==2023.7.0" \
        "pyarrow==8.0.0" \
        "pyliftover==0.4" \
        "scikit-learn==1.2.2" \
        "scipy==1.10.1" \
        "xgboost==1.7.5" \
        "requests==2.31.0" && \
    chmod +x /opt/flames/FLAMES.py && \
    rm -rf /opt/flames/.git && \
    apt-get remove -y git && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/flames:$PATH"
ENV PYTHONUNBUFFERED=1
