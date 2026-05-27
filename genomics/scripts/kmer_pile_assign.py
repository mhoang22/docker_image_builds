#!/usr/bin/env python3
"""
kmer_pile_assign.py
-------------------
Split protein sequences from two FASTA files into k-mers and assign each
k-mer to one of three piles:
  - shared      : present in both file1 and file2
  - only_file1  : present only in file1
  - only_file2  : present only in file2

Each k-mer records which transcript(s) it came from per source file.

Usage:
    python kmer_pile_assign.py \\
        --file1  proteins_A.fa \\
        --file2  proteins_B.fa \\
        --k      8 9 10 11 \\
        --outdir ./kmer_output

    # Override the pile label names in the output:
    python kmer_pile_assign.py \\
        --file1  short_read.fa  --label1 short \\
        --file2  long_read.fa   --label2 long  \\
        --k      9 \\
        --outdir ./kmer_output

Output (one TSV per k size):
    kmer_output/kmers_k{k}.tsv
    Columns:
        kmer | pile | file1_transcripts | file2_transcripts | file1_count | file2_count
    (or using custom labels:)
        kmer | pile | short_transcripts | long_transcripts  | short_count | long_count
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AMBIGUOUS_AA = set("XBZUO*")   # skip k-mers containing these by default


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_fasta(filepath: str, label: str) -> dict[str, str]:
    """
    Parse a FASTA file → {record_id: sequence}.
    Sequences are uppercased; trailing stop-codon asterisks are stripped.
    Empty sequences are skipped with a warning.
    """
    records: dict[str, str] = {}
    duplicates = 0

    print(f"[{label}] Reading {filepath} ...")
    for rec in SeqIO.parse(filepath, "fasta"):
        seq = str(rec.seq).upper().strip("*")
        if not seq:
            print(f"  [{label}] WARNING: '{rec.id}' has an empty sequence — skipped.")
            continue
        if rec.id in records:
            duplicates += 1
        records[rec.id] = seq   # last one wins on duplicate IDs

    print(f"  [{label}] Loaded {len(records):,} sequences.", end="")
    if duplicates:
        print(f"  WARNING: {duplicates:,} duplicate IDs found (last sequence kept).", end="")
    print()
    return records


def check_id_overlap(seqs1: dict, seqs2: dict, label1: str, label2: str) -> None:
    """
    Warn if the two files share transcript IDs.
    This won't break anything, but overlapping IDs are ambiguous in the output.
    """
    overlap = set(seqs1.keys()) & set(seqs2.keys())
    if overlap:
        print(
            f"\nWARNING: {len(overlap):,} transcript ID(s) appear in both files "
            f"(e.g. '{next(iter(overlap))}').\n"
            f"  The '{label1}_transcripts' and '{label2}_transcripts' columns in the\n"
            f"  output will both show the same ID — they cannot be distinguished.\n"
            f"  Consider prefixing your FASTA headers before running this script.\n"
        )


def extract_kmers(
    sequences: dict[str, str],
    k: int,
    label: str,
    skip_ambiguous: bool = True,
) -> dict[str, set[str]]:
    """
    Sliding-window k-mer extraction over all sequences.

    Returns:
        kmer_to_transcripts: {kmer_string: {transcript_id, ...}}
    """
    kmer_to_transcripts: dict[str, set[str]] = defaultdict(set)
    skipped_short = 0
    skipped_ambig = 0

    for tid, seq in tqdm(sequences.items(), desc=f"  [{label}] k={k}", unit="seq"):
        if len(seq) < k:
            skipped_short += 1
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if skip_ambiguous and AMBIGUOUS_AA.intersection(kmer):
                skipped_ambig += 1
                continue
            kmer_to_transcripts[kmer].add(tid)

    notes = []
    if skipped_short:
        notes.append(f"{skipped_short:,} seqs shorter than k skipped")
    if skipped_ambig:
        notes.append(f"{skipped_ambig:,} ambiguous k-mers skipped")
    if notes:
        print(f"    ({'; '.join(notes)})")

    return kmer_to_transcripts


def assign_piles(
    kmers1: dict[str, set[str]],
    kmers2: dict[str, set[str]],
    label1: str,
    label2: str,
) -> pd.DataFrame:
    """
    Categorise every unique k-mer into shared / only_<label1> / only_<label2>.

    Returns a DataFrame with columns:
        kmer | pile | {label1}_transcripts | {label2}_transcripts
             | {label1}_count | {label2}_count
    """
    all_kmers = set(kmers1.keys()) | set(kmers2.keys())

    rows = []
    for kmer in all_kmers:
        in1 = kmer in kmers1
        in2 = kmer in kmers2

        if in1 and in2:
            pile = "shared"
        elif in1:
            pile = f"only_{label1}"
        else:
            pile = f"only_{label2}"

        t1 = ";".join(sorted(kmers1[kmer])) if in1 else ""
        t2 = ";".join(sorted(kmers2[kmer])) if in2 else ""
        c1 = len(kmers1[kmer]) if in1 else 0
        c2 = len(kmers2[kmer]) if in2 else 0

        rows.append((kmer, pile, t1, t2, c1, c2))

    df = pd.DataFrame(
        rows,
        columns=[
            "kmer", "pile",
            f"{label1}_transcripts", f"{label2}_transcripts",
            f"{label1}_count",       f"{label2}_count",
        ],
    )

    pile_order = {"shared": 0, f"only_{label1}": 1, f"only_{label2}": 2}
    df["_order"] = df["pile"].map(pile_order)
    df = df.sort_values(["_order", "kmer"]).drop(columns="_order").reset_index(drop=True)
    return df


def print_summary(df: pd.DataFrame, k: int, label1: str, label2: str) -> None:
    piles = ["shared", f"only_{label1}", f"only_{label2}"]
    counts = df["pile"].value_counts()
    total  = len(df)
    print(f"\n  ── k={k} summary {'─'*30}")
    for pile in piles:
        n   = counts.get(pile, 0)
        pct = 100 * n / total if total else 0
        print(f"  {pile:<20}: {n:>10,}  ({pct:.1f}%)")
    print(f"  {'TOTAL':<20}: {total:>10,}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Assign protein k-mers to shared / only_file1 / only_file2 piles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--file1",  required=True, help="First input FASTA file")
    parser.add_argument("--file2",  required=True, help="Second input FASTA file")
    parser.add_argument(
        "--label1", default="file1",
        help="Short label for file1 used in column names and pile names (default: file1)"
    )
    parser.add_argument(
        "--label2", default="file2",
        help="Short label for file2 used in column names and pile names (default: file2)"
    )
    parser.add_argument(
        "--k", nargs="+", type=int, default=[9], metavar="K",
        help="One or more k-mer sizes (e.g. --k 8 9 10 11). Default: 9"
    )
    parser.add_argument("--outdir", default=".", help="Output directory (created if absent)")
    parser.add_argument(
        "--keep-ambiguous", action="store_true",
        help="Include k-mers containing ambiguous amino acids (X/B/Z/U/O/*). Default: skip."
    )
    args = parser.parse_args()

    # --- Validate ---
    for k in args.k:
        if k < 1:
            sys.exit(f"ERROR: k must be >= 1, got {k}.")
        if k > 30:
            print(f"WARNING: k={k} is very large — expect very sparse results.")

    if args.label1 == args.label2:
        sys.exit("ERROR: --label1 and --label2 must be different strings.")

    os.makedirs(args.outdir, exist_ok=True)
    skip_amb = not args.keep_ambiguous

    # --- Parse FASTAs once (reused across all k) ---
    seqs1 = parse_fasta(args.file1, args.label1)
    seqs2 = parse_fasta(args.file2, args.label2)

    check_id_overlap(seqs1, seqs2, args.label1, args.label2)

    # --- Per-k processing ---
    for k in args.k:
        print(f"\n{'='*55}")
        print(f"  Processing k = {k}")
        print(f"{'='*55}")

        kmers1 = extract_kmers(seqs1, k, args.label1, skip_ambiguous=skip_amb)
        kmers2 = extract_kmers(seqs2, k, args.label2, skip_ambiguous=skip_amb)

        print(
            f"  Unique k-mers — {args.label1}: {len(kmers1):,}"
            f"  |  {args.label2}: {len(kmers2):,}"
        )

        df = assign_piles(kmers1, kmers2, args.label1, args.label2)
        print_summary(df, k, args.label1, args.label2)

        outpath = os.path.join(args.outdir, f"kmers_k{k}.tsv")
        df.to_csv(outpath, sep="\t", index=False)
        print(f"  Saved → {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
