FROM ghcr.io/tron-bioinformatics/splice2neo:v0.6.13

RUN Rscript -e 'BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10"), Ncpus=8)'
RUN Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38"), Ncpus=8)'
RUN Rscript -e 'BiocManager::install(c("EnsDb.Mmusculus.v95"), Ncpus=8)'
