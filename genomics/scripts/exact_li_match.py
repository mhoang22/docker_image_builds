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
        --out       /data/output/exact_match.tsv \\
        --log       /data/output/exact_match.log
"""
import argparse
import sys
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


def match_peptides(peptides: list, proteins: dict) -> pd.DataFrame:
    results = []
    for pep in peptides:
        norm_pep = normalize(pep)
        hits = [pid for pid, (_, norm_seq) in proteins.items()
                if norm_pep in norm_seq]
        results.append({
            "peptide":     pep,
            "n_proteins":  len(hits),
            "protein_ids": ";".join(hits) if hits else "NO_MATCH",
        })
    df = pd.DataFrame(results)
    matched = (df["protein_ids"] != "NO_MATCH").sum()
    print(f"{matched:,}/{len(peptides):,} peptides matched "
          f"({matched / len(peptides) * 100:.1f}%)")
    return df


class Tee:
    """Write to both stdout and a log file simultaneously."""
    def __init__(self, log_path: str):
        self.terminal = sys.stdout
        self.log = open(log_path, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True, help="One peptide per line")
    ap.add_argument("--fasta",    required=True, help="Reference FASTA")
    ap.add_argument("--out",      default="/data/output/exact_match.tsv")
    ap.add_argument("--log",      default=None,  help="Optional path to save screen output")
    args = ap.parse_args()

    if args.log:
        sys.stdout = Tee(args.log)

    try:
        with open(args.peptides) as fh:
            peptides = [ln.strip() for ln in fh if ln.strip()]
        print(f"Loading input peptides: {args.peptides}")
        print(f"Loaded {len(peptides):,} query peptides")

        proteins = build_protein_index(args.fasta)
        df = match_peptides(peptides, proteins)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"Saved → {args.out}")

    finally:
        if args.log and hasattr(sys.stdout, "close"):
            log_path = args.log
            sys.stdout.close()
            sys.stdout = sys.__stdout__
            print(f"Log saved → {log_path}")
