# processor-cytonorm: CytoNorm batch normalisation of arcsinh-transformed FCS files.
# R 4.6.1 (rocker/r-ver, Ubuntu 24.04), binary packages from Posit Public Package Manager
# pinned to one snapshot date, Bioconductor 3.23, CytoNorm pinned to a GitHub commit.
FROM --platform=linux/amd64 rocker/r-ver:4.6.1@sha256:bafb3b9a72791eeeb810f0bc0472e74c7f4939b20b46c96cb7ebc6838b6a25ef

ARG P3M_SNAPSHOT=2026-10-26
ARG BIOC_VERSION=3.23
# saeyslab/CytoNorm main as of 2026-07-07 (DESCRIPTION Version 2.0.12)
ARG CYTONORM_SHA=278f352d8edffd7a5d3e6b0b76c402d991a1e7a9

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    OMP_NUM_THREADS=1

# Runtime libraries the binary packages link against (igraph: glpk, gmp; XML: libxml2).
RUN apt-get update \
 && apt-get install -y --no-install-recommends libglpk40 libgmp10 libxml2 \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "\
  options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/noble/${P3M_SNAPSHOT}'), \
          BioC_mirror = 'https://p3m.dev/bioconductor/__linux__/noble/${P3M_SNAPSHOT}', \
          Ncpus = 4); \
  install.packages(c('BiocManager', 'remotes')); \
  BiocManager::install(c('flowCore', 'FlowSOM'), version = '${BIOC_VERSION}', update = FALSE, ask = FALSE); \
  install.packages(c('emdist', 'dplyr', 'stringr', 'pheatmap', 'ggplot2', 'ggpubr', 'gridExtra', \
                     'ggridges', 'tidyr', 'jsonlite')); \
  remotes::install_github('saeyslab/CytoNorm@${CYTONORM_SHA}', dependencies = FALSE, upgrade = 'never'); \
  for (p in c('flowCore', 'FlowSOM', 'CytoNorm', 'jsonlite', 'ggplot2')) \
    if (!requireNamespace(p, quietly = TRUE)) stop('package missing: ', p)" \
 && rm -rf /tmp/downloaded_packages /tmp/Rtmp* /root/.cache

WORKDIR /app
COPY VERSION /app/VERSION
COPY R/ /app/R/

# No network at run time; everything the run needs is in the image.
ENV INPUT_DIR=/data/input \
    OUTPUT_DIR=/data/output

ENTRYPOINT ["Rscript", "/app/R/main.R"]
