# =============================================================================
# genomics — Docker image for de novo peptide sequencing & proteomics analysis
#
# Includes:
#   - BLAST+ (blastp / makeblastdb)
#   - Python 3.11 + BioPython, pandas, numpy
#   - Exact L/I-aware peptide matching scripts
#   - Jupyter Lab for interactive analysis
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="myhoang04"
LABEL description="Genomics, Proteomics / de novo peptide analysis: BLAST+, BioPython, pandas"

# ---------------------------------------------------------------------------
# 0. Non-interactive installs
# ---------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        curl \
        git \
        unzip \
        bzip2 \
        ca-certificates \
        build-essential \
        gcc \
        g++ \
        make \
        ncbi-blast+ \
        python3.11 \
        python3.11-dev \
        python3-pip \
        python3.11-venv \
        libxml2-dev \
        libxslt1-dev \
        libhdf5-dev \
        less \
        vim \
        parallel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Make python3.11 the default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
 && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
 && python -m pip install --upgrade pip setuptools wheel

# ---------------------------------------------------------------------------
# 2. Python packages
# ---------------------------------------------------------------------------
RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.1 \
    biopython==1.83 \
    "pyteomics[xml,unimod]==4.7.3" \
    lxml==5.2.2 \
    sqlalchemy==2.0.30 \
    pynumpress==0.0.8 \
    openpyxl==3.1.2 \
    xlrd==2.0.1 \
    matplotlib==3.9.0 \
    seaborn==0.13.2 \
    plotly==5.22.0 \
    kaleido==0.2.1 \
    jupyterlab==4.2.1 \
    ipywidgets==8.1.3 \
    tqdm==4.66.4 \
    rich==13.7.1 \
    click==8.1.7 \
    requests==2.32.3 \
    pyahocorasick==2.1.0 \
    scikit-learn==1.5.0 \
    statsmodels==0.14.2

# ---------------------------------------------------------------------------
# 3. Copy analysis scripts into image
#
#   scripts/
#   |-- blast_peptides.sh
#   |     Shell wrapper that runs the full BLAST pipeline end-to-end:
#   |     list_to_fasta -> makeblastdb -> blastp -> parse_blast_results.
#   |     Usage: bash blast_peptides.sh <peptide_list.txt> <ref_fasta> <out_dir>
#   |
#   |-- exact_li_match.py
#   |     L/I-aware exact string matching: treats Leucine and Isoleucine as
#   |     equivalent (indistinguishable by mass) when comparing peptide lists.
#   |     Usage: python exact_li_match.py <query_list.txt> <ref_list.txt>
#   |
#   |-- kmer_pile_assign.py
#   |     Assigns k-mer "piles" (groups of overlapping k-mers) to peptide
#   |     sequences; used for coverage and redundancy analysis.
#   |     Usage: python kmer_pile_assign.py <peptides.txt> --k <int>
#   |
#   |-- kmer_ref_annotate.py
#   |     Annotates peptides against a k-mer reference database, reporting
#   |     which reference proteins share k-mers with each query peptide.
#   |     Usage: python kmer_ref_annotate.py <peptides.txt> <ref_fasta>
#   |
#   |-- list_to_fasta.py
#   |     Converts a plain text peptide list (one sequence per line) into a
#   |     FASTA file suitable for BLAST or other sequence tools.
#   |     Usage: python list_to_fasta.py <peptide_list.txt> <output.fasta>
#   |
#   |-- parse_blast_results.py
#   |     Parses BLAST tabular output (fmt 6), filters by e-value / identity,
#   |     and writes a tidy TSV summary with hit annotations.
#   |     Usage: python parse_blast_results.py <blast_out.txt> <output.tsv>
#   |
#   +-- plot_instanovo_match_vs_nomatch.py
#         Joins InstaNovo+ predictions against a peptide reference file
#         (columns: peptide, n_proteins, protein_ids), classifies each
#         prediction as MATCH (n_proteins > 0) or NO_MATCH, then produces
#         side-by-side box plots of log probability and delta_mass_ppm.
#         Usage: python plot_instanovo_match_vs_nomatch.py \
#                    --instanovo <instanovo.csv> \
#                    --reference <ref.csv> \
#                    --out       <plot.png>
# ---------------------------------------------------------------------------
WORKDIR /opt/genomics

COPY scripts/ /opt/genomics/
RUN chmod +x /opt/genomics/blast_peptides.sh

ENV PATH="/opt/genomics:${PATH}"

# ---------------------------------------------------------------------------
# 4. Working directories
# ---------------------------------------------------------------------------
RUN mkdir -p /data/input \
             /data/output \
             /data/fasta \
             /data/blast_db \
             /data/blast_results

WORKDIR /data

# ---------------------------------------------------------------------------
# 5. Verify installations
# ---------------------------------------------------------------------------
RUN python -c "import Bio; import pandas; import numpy; import scipy; import pyteomics; import lxml; import sqlalchemy; import matplotlib; import sklearn; import ahocorasick; print('All Python packages OK')" \
 && blastp -version \
 && makeblastdb -version

# ---------------------------------------------------------------------------
# 6. Default command — open Jupyter Lab (override for batch runs)
# ---------------------------------------------------------------------------
EXPOSE 8888

CMD ["jupyter", "lab", \
     "--ip=0.0.0.0", \
     "--port=8888", \
     "--no-browser", \
     "--allow-root", \
     "--NotebookApp.token=''", \
     "--NotebookApp.password=''", \
     "--notebook-dir=/data"]
