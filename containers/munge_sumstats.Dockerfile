FROM bioconductor/bioconductor_docker:RELEASE_3_21

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R --no-echo -e "install.packages( \
    c('vroom', 'logger', 'argparse', 'box', 'data.table'), \
    repos = 'https://cloud.r-project.org' \
)"

# Install Bioconductor packages
RUN R --no-echo -e "BiocManager::install(c( \
    'MungeSumstats', \
    'SNPlocs.Hsapiens.dbSNP155.GRCh38', \
    'BSgenome.Hsapiens.NCBI.GRCh38', \
    'BSgenome.Hsapiens.1000genomes.hs37d5' \
), ask = FALSE, update = FALSE)"

ENV LC_ALL=C
ENV LANG=C
ENV OPENSSL_CONF=/dev/null

LABEL org.opencontainers.image.description="MungeSumstats with dbSNP155 for GWAS summary statistic harmonization"
LABEL org.opencontainers.image.source="https://github.com/ampregnall/nf-meta-gwas"
