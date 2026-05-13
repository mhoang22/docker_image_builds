#!/usr/bin/env python3
"""
parse_blast_results.py
----------------------
Filter raw BLAST tabular output and apply L/I equivalence as a final check.

Usage:
    python parse_blast_results.py \\
        --tool      casanovo \\
        --blast_tsv /data/blast_results/casanovo_blast_raw.tsv \\
        --min_pident 90 \\
        --max_mismatch 1 \\
        --out       /data/output/casanovo_blast_mapped.tsv
"""

import argparse
from collections import defaultdict
import pandas as pd


COLS = ["query_id", "subject_id", "pident", "length",
        "mismatch", "qseq", "sseq", "evalue", "bitscore"]


def normalize(seq: str) -> str:
    return seq.upper().replace("I", "L")


def is_li_only_mismatch(row) -> bool:
    """Return True when all mismatches are explained by L↔I swaps."""
    if row["mismatch"] == 0:
        return True
    return normalize(row["qseq"]) == normalize(row["sseq"])


def parse_blast(tsv: str, tool: str,
                min_pident: float, max_mismatch: int) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t", header=None, names=COLS)
    print(f"[{tool}] Raw BLAST hits : {len(df):,}")

    # Alignment must cover the full query (no partial hits)
    df["query_len"] = df["qseq"].str.len()
    df = df[df["length"] == df["query_len"]]
    print(f"[{tool}] Full-length aln : {len(df):,}")

    # Identity / mismatch hard filters
    df = df[(df["pident"] >= min_pident) & (df["mismatch"] <= max_mismatch)]
    print(f"[{tool}] After pident/mm  : {len(df):,}")

    # L/I equivalence
    df = df[df.apply(is_li_only_mismatch, axis=1)].copy()
    print(f"[{tool}] After L/I filter : {len(df):,}")

    df["tool"] = tool
    return df[["tool", "query_id", "subject_id",
               "pident", "mismatch", "qseq", "sseq", "evalue"]]


def summarize(df: pd.DataFrame, tool: str) -> pd.DataFrame:
    mapping = defaultdict(set)
    for _, row in df.iterrows():
        mapping[row["qseq"]].add(row["subject_id"])
    rows = [{"tool": tool, "peptide": pep,
             "n_proteins": len(prots),
             "protein_ids": ";".join(sorted(prots))}
            for pep, prots in mapping.items()]
    return pd.DataFrame(rows).sort_values("peptide")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool",       required=True)
    ap.add_argument("--blast_tsv",  required=True)
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--max_mismatch", type=int, default=1)
    ap.add_argument("--out",        default=None)
    args = ap.parse_args()

    df = parse_blast(args.blast_tsv, args.tool,
                     args.min_pident, args.max_mismatch)
    summary = summarize(df, args.tool)

    outfile = args.out or f"/data/output/{args.tool}_blast_mapped.tsv"
    summary.to_csv(outfile, sep="\t", index=False)
    print(f"\nSaved {len(summary):,} peptide mappings → {outfile}")
    print(summary.head(10).to_string(index=False))
