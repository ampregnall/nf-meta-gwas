FROM registry.access.redhat.com/ubi9/ubi:9.4

LABEL maintainer="nf-meta-gwas"
LABEL description="gwaslab environment for GWAS summary statistics processing"

RUN dnf install -y \
        python3.11 \
        python3.11-devel \
        gcc \
        gcc-c++ \
        make \
        bzip2 \
        xz \
        zlib-devel \
        bzip2-devel \
        xz-devel \
        libcurl-devel \
        openssl-devel \
    && dnf clean all

# Build htslib (provides tabix and bgzip)
ARG HTSLIB_VERSION=1.21
RUN curl -fsSL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    | tar -xj && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION}

# Build bcftools against installed htslib
ARG BCFTOOLS_VERSION=1.21
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    | tar -xj && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf bcftools-${BCFTOOLS_VERSION}

RUN python3.11 -m ensurepip && \
    python3.11 -m pip install --no-cache-dir --upgrade pip && \
    python3.11 -m pip install --no-cache-dir "gwaslab>=4.0"

# Make python3 resolve to python3.11 for script shebangs
RUN ln -sf /usr/bin/python3.11 /usr/local/bin/python3

# Verify installations
RUN python3 -c "from importlib.metadata import version; print('gwaslab', version('gwaslab'), 'installed successfully')"
RUN bcftools --version | head -1
RUN tabix --version | head -1

ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp/matplotlib
