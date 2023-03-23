FROM r-base:4.2.2
LABEL maintainer="Andre Fonseca" \
    description="scPackages - Container with single-cell R dependencies"

# Install ps, for Nextflow. https://www.nextflow.io/docs/latest/tracing.html
RUN apt-get update && \
    apt-get install -y procps \
    pandoc \
    libcurl4-openssl-dev \
    r-cran-curl \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \ 
    libfribidi-dev \
    libfreetype6-dev \ 
    libpng-dev \ 
    libtiff5-dev \
    libjpeg-dev

# Custom cache invalidation
ARG CACHEBUST=1

# Install required R packages
ARG R_DEPS="c( \
    'pkgdown', \
    'patchwork', \
    'devtools',\
    'rmarkdown', \
    'tidyverse', \
    'bookdown', \
    'BiocManager', \
    'Seurat', \
    'harmony',\
    'immunarch', \
    'ggpubr', \
    'scater', \
    'viridis', \
    'RcmdrMisc' \
    )"

ARG R_BIOC_DEPS="c( \
    'Biobase', \
    'singleCellTK', \
    'SingleCellExperiment', \
    'DropletUtils' \ 
    )"

ARG DEV_DEPS="c(\
    'PaulingLiu/ROGUE',\
    'immunogenomics/lisi',\
    'theislab/kBET'\
    )"

# Caching R-lib on the building process
RUN --mount=type=cache,target=/usr/local/lib/R Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
RUN --mount=type=cache,target=/usr/local/lib/R Rscript -e "BiocManager::install(${R_BIOC_DEPS})"
RUN --mount=type=cache,target=/usr/local/lib/R Rscript -e "devtools::install_github(${DEV_DEPS})"
    
CMD ["R"]
