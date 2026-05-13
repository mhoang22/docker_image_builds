#!/usr/bin/env python3
"""
list_to_fasta.py
----------------
Convert a plain peptide list (one per line) to a FASTA file
suitable for use as a BLAST query.

Usage:
    python list_to_fasta.py \\
        --peptides /data/input/casanovo_peptides.txt \\
        --tool     casanovo \\
        --out      /data/blast_results/casanovo_queries.fasta
"""

import argparse


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True, help="One peptide per line")
    ap.add_argument("--tool",     required=True, help="Tool name (used as FASTA ID prefix)")
    ap.add_argument("--out",      default=None)
    args = ap.parse_args()

    with open(args.peptides) as fh:
        seqs = [ln.strip() for ln in fh if ln.strip()]

    outfile = args.out or f"/data/blast_results/{args.tool}_queries.fasta"
    with open(outfile, "w") as fh:
        for i, seq in enumerate(seqs, 1):
            fh.write(f">{args.tool}_peptide_{i}\n{seq}\n")

    print(f"Written {len(seqs):,} sequences → {outfile}")
