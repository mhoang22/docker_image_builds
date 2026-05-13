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
        # Core tools
        wget \
        curl \
        git \
        unzip \
        bzip2 \
        ca-certificates \
        # Build tools (needed by some Python packages)
        build-essential \
        gcc \
        g++ \
        make \
        # NCBI BLAST+
        ncbi-blast+ \
        # Python
        python3.11 \
        python3.11-dev \
        python3-pip \
        python3.11-venv \
        # Misc utilities
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
    # Core data science
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.1 \
    # Bioinformatics
    biopython==1.83 \
    pyteomics==4.7.3 \
    # File I/O
    openpyxl==3.1.2 \
    xlrd==2.0.1 \
    # Visualization
    matplotlib==3.9.0 \
    seaborn==0.13.2 \
    plotly==5.22.0 \
    # Interactive analysis
    jupyterlab==4.2.1 \
    ipywidgets==8.1.3 \
    # Utilities
    tqdm==4.66.4 \
    rich==13.7.1 \
    click==8.1.7 \
    # String matching
    pyahocorasick==2.1.0

# ---------------------------------------------------------------------------
# 3. Copy analysis scripts into image
# ---------------------------------------------------------------------------
WORKDIR /opt/genomics

# --- Script 1: exact L/I-aware matching ------------------------------------
COPY scripts/exact_li_match.py       /opt/genomics/exact_li_match.py

# --- Script 2: BLAST result parser -----------------------------------------
COPY scripts/parse_blast_results.py  /opt/genomics/parse_blast_results.py

# --- Script 3: FASTA helper (peptide list → FASTA) -------------------------
COPY scripts/list_to_fasta.py        /opt/genomics/list_to_fasta.py

# --- Shell wrapper for full BLAST pipeline ---------------------------------
COPY scripts/blast_peptides.sh       /opt/genomics/blast_peptides.sh
RUN chmod +x /opt/genomics/blast_peptides.sh

# Add scripts to PATH
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
RUN python -c "import Bio; import pandas; import numpy; print('Python packages OK')" \
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
