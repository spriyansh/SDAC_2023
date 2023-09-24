# Use rocker/rstudio as the base image
FROM rocker/rstudio:4

# Maintainer
LABEL maintainer = Priyansh Srivastava <spriyansh29@gmail.com>

# Install Required OS libs
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libpng-dev \
    gdal-bin \
    libgdal-dev \
    libfontconfig1-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk40 \
    texlive \
    texlive-latex-extra \
    texlive-fonts-recommended \
    libxt6 \
    libcairo2-dev \
    patch \
    libgeos-dev \
    libxt-dev \
    cmake \
    libgsl-dev \
    build-essential \
    libfftw3-dev \
    libmagick++-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Installation in batches to avoid error-code 137
RUN Rscript -e 'install.packages(c("tidyverse", "BiocManager", "metap", "intergraph", "spdep", "eHOF", "devtools", "doParallel", "MatrixExtra", "Seurat", "ggpubr", "ggrastr", "Cairo", "cellranger", "mutoss", "rvest", "xml2", "lava", "qqconf", "revealjs", "paletteer", "ggpubr", "markdown", "magick"), repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'

# BioConductor
RUN Rscript -e "BiocManager::install('qvalue', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('multtest', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('BiocGenerics', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DelayedArray', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DelayedMatrixStats', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('limma', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('S4Vectors', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SingleCellExperiment', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SummarizedExperiment', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('batchelor', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('Matrix.utils', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('rhdf5', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('HDF5Array', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('slingshot', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('miloR', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('tradeSeq', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    #Rscript -e "BiocManager::install('celldex', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('scater', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('glmGamPoi', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SingleR', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    #Rscript -e "BiocManager::install('dyngen', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    #Rscript -e "BiocManager::install('splatter', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DESeq2', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('scran', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('limma', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('org.Hs.eg.db', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('GO.db', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('BiocStyle', update = TRUE, ask = FALSE, Ncpus = 24)"

# Github
RUN Rscript -e 'remotes::install_github(c("rspatial/terra", "satijalab/seurat-data", "satijalab/azimuth")); \
    devtools::install_github(c("cole-trapnell-lab/leidenbase", "cole-trapnell-lab/monocle3", "mojaveazure/seurat-disk", "ConesaLab/RColorConesa"))'


# Monocle3 Compatible and Azimuth Reference
#RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); remotes::install_version("igraph", version = "1.4.3", force = TRUE)'
RUN Rscript -e 'install.packages("bonemarrowref.SeuratData_1.0.0.tar.gz", repos = NULL, Ncpus = 4, quiet = TRUE, type = "source")'


# Start the lighter stage
#FROM rocker/rstudio:4

# Copy the R library directory from the builder to the final image
#COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

RUN Rscript -e 'install.packages(c("wesanderson"), repos = "https://cloud.r-project.org/", dependencies = TRUE, Ncpus = 24, quiet = TRUE)'

# Create ESR data
RUN mkdir /home/esr

# Expose port 8787 (RStudio runs on this port)
EXPOSE 8787

CMD ["/init"]


