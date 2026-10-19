FROM bioconductor/bioconductor_docker:RELEASE_3_23-R-4.6.1

LABEL software="scRepertoire"
LABEL software.version="2.9.1"
LABEL about.summary="scRepertoire - a toolkit for single-cell immune receptor profiling"
LABEL about.home="https://github.com/BorchLab/scRepertoire"
LABEL about.documentation="https://www.borch.dev/uploads/scRepertoire/"

# Bioconductor 3.23 ships scRepertoire 2.8.0 / immApex 1.6.0, but scRepertoire 2.9.1
# needs immApex >= 1.7.1, so both are installed from pinned GitHub tags.
# Bioconductor/CRAN dependencies come from the 3.23 release repos.
RUN Rscript -e 'BiocManager::install(c("remotes", "BiocParallel", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "scater"), Ncpus=8, ask=FALSE, update=FALSE)'

RUN Rscript -e 'remotes::install_github("BorchLab/immApex@1.7.1", repos=BiocManager::repositories(), Ncpus=8, upgrade="never")'

RUN Rscript -e 'remotes::install_github("BorchLab/scRepertoire@v2.9.1", repos=BiocManager::repositories(), Ncpus=8, upgrade="never")'

# Commonly used Suggests (Seurat workflows, dowser export, plotting extras)
RUN Rscript -e 'BiocManager::install(c("Seurat", "dowser", "circlize", "scales", "Peptides"), Ncpus=8, ask=FALSE, update=FALSE)'

RUN Rscript -e 'library(scRepertoire); stopifnot(packageVersion("scRepertoire") == "2.9.1"); stopifnot(packageVersion("immApex") >= "1.7.1"); library(Seurat); library(dowser)'

CMD ["/bin/bash"]
