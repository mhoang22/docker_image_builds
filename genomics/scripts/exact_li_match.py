#!/usr/bin/env python3
"""
exact_li_match.py
-----------------
Exact substring matching with L=I equivalence.
Every query peptide is searched as a substring of each full-length protein.

Usage:
    python exact_li_match.py \\
        --peptides  /data/input/casanovo_peptides.txt \\
        --fasta     /data/fasta/reference.fasta \\
        --tool      casanovo \\
        --out       /data/output/casanovo_exact_match.tsv
"""

import argparse
from collections import defaultdict
from Bio import SeqIO
import pandas as pd


def normalize(seq: str) -> str:
    """Collapse I→L so both are treated as identical."""
    return seq.upper().replace("I", "L")


def build_protein_index(fasta_path: str) -> dict:
    """Load FASTA and return {prot_id: (orig_seq, norm_seq)}."""
    print(f"Loading reference FASTA: {fasta_path}")
    proteins = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        orig = str(rec.seq).upper()
        proteins[rec.id] = (orig, normalize(orig))
    print(f"  Loaded {len(proteins):,} proteins")
    return proteins


def match_peptides(peptides: list, proteins: dict, tool: str) -> pd.DataFrame:
    results = []
    for pep in peptides:
        norm_pep = normalize(pep)
        hits = [pid for pid, (_, norm_seq) in proteins.items()
                if norm_pep in norm_seq]
        results.append({
            "tool":       tool,
            "peptide":    pep,
            "n_proteins": len(hits),
            "protein_ids": ";".join(hits) if hits else "NO_MATCH",
        })

    df = pd.DataFrame(results)
    matched = (df["protein_ids"] != "NO_MATCH").sum()
    print(f"[{tool}] {matched:,}/{len(peptides):,} peptides matched "
          f"({matched / len(peptides) * 100:.1f}%)")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True, help="One peptide per line")
    ap.add_argument("--fasta",    required=True, help="Reference FASTA")
    ap.add_argument("--tool",     required=True,
                    help="Tool name (casanovo | deepnovo | instanovo | cascadia)")
    ap.add_argument("--out",      default=None)
    args = ap.parse_args()

    with open(args.peptides) as fh:
        peptides = [ln.strip() for ln in fh if ln.strip()]
    print(f"Loaded {len(peptides):,} query peptides")

    proteins = build_protein_index(args.fasta)
    df = match_peptides(peptides, proteins, args.tool)

    outfile = args.out or f"/data/output/{args.tool}_exact_match.tsv"
    df.to_csv(outfile, sep="\t", index=False)
    print(f"Saved → {outfile}")
