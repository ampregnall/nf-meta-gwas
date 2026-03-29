FROM rocker/tidyverse:4.4

LABEL maintainer="nf-meta-gwas"
LABEL description="R environment for GWAS meta-analysis and eQTL/pQTL colocalization with HyPrColoc"

# Install system dependencies for package compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install arrow using Posit Package Manager binary for faster install
RUN Rscript -e "install.packages('arrow', repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')"

# Install CRAN packages (vroom, glue, dplyr already in rocker/tidyverse)
RUN Rscript -e "install.packages(c('logger', 'argparse', 'box'), repos = 'https://cloud.r-project.org')"

# Install specific RcppEigen version required by HyPrColoc
RUN Rscript -e "remotes::install_version('RcppEigen', version = '0.3.3.9.3', repos = 'https://cloud.r-project.org')"

# Install HyPrColoc from GitHub
RUN Rscript -e "remotes::install_github('jrs95/hyprcoloc', build_vignettes = FALSE)"
