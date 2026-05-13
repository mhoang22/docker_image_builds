#!/usr/bin/env bash
# =============================================================================
# blast_peptides.sh
# Full BLAST pipeline: peptide list → makeblastdb → blastp → raw TSV
#
# Usage:
#   bash blast_peptides.sh <tool_name> <peptides.txt> <reference.fasta>
#
# Example:
#   bash blast_peptides.sh casanovo /data/input/casanovo_peptides.txt \
#                                   /data/fasta/reference.fasta
# =============================================================================

set -euo pipefail

TOOL="${1:?Usage: $0 <tool_name> <peptides.txt> <reference.fasta>}"
PEPTIDES="${2:?Missing peptides file}"
FASTA="${3:?Missing reference FASTA}"

BLAST_DIR="/data/blast_results"
DB_DIR="/data/blast_db"
mkdir -p "$BLAST_DIR" "$DB_DIR"

echo "=========================================="
echo " Tool     : $TOOL"
echo " Peptides : $PEPTIDES"
echo " FASTA    : $FASTA"
echo "=========================================="

# ---------------------------------------------------------------------------
# Step 1 — Convert peptide list to query FASTA
# ---------------------------------------------------------------------------
echo "[1/3] Converting peptide list to FASTA..."
python /opt/genomics/list_to_fasta.py \
    --peptides "$PEPTIDES" \
    --tool     "$TOOL" \
    --out      "${BLAST_DIR}/${TOOL}_queries.fasta"

# ---------------------------------------------------------------------------
# Step 2 — Build BLAST protein database (idempotent)
# ---------------------------------------------------------------------------
DB="${DB_DIR}/fullpep_db"
if [ ! -f "${DB}.pin" ]; then
    echo "[2/3] Building BLAST database..."
    makeblastdb \
        -in          "$FASTA" \
        -dbtype      prot \
        -out         "$DB" \
        -parse_seqids
else
    echo "[2/3] BLAST database already exists, skipping."
fi

# ---------------------------------------------------------------------------
# Step 3 — Run blastp (short peptide settings)
# ---------------------------------------------------------------------------
echo "[3/3] Running blastp..."
blastp \
    -task            blastp-short \
    -query           "${BLAST_DIR}/${TOOL}_queries.fasta" \
    -db              "$DB" \
    -out             "${BLAST_DIR}/${TOOL}_blast_raw.tsv" \
    -outfmt          "6 qseqid sseqid pident length mismatch qseq sseq evalue bitscore" \
    -evalue          1 \
    -word_size       2 \
    -matrix          BLOSUM62 \
    -num_threads     "$(nproc)" \
    -max_target_seqs 10

echo ""
echo "Done. Raw BLAST output → ${BLAST_DIR}/${TOOL}_blast_raw.tsv"
echo ""
echo "Next step — filter and apply L/I equivalence:"
echo "  python /opt/genomics/parse_blast_results.py \\"
echo "      --tool $TOOL \\"
echo "      --blast_tsv ${BLAST_DIR}/${TOOL}_blast_raw.tsv"
