#!/usr/bin/env python3
"""
kmer_ref_annotate.py
--------------------
Steps 1–3 of the peptidome novelty analysis pipeline.

Takes the output TSV(s) from kmer_pile_assign.py (which compared short-read
vs long-read peptidomes) and annotates each k-mer with whether it has an
exact match in a reference proteome FASTA.

This produces a final 6-category classification per k-mer:

    pile (short/long/shared) × ref_match (True/False)
    ─────────────────────────────────────────────────
    shared       + in_reference   → canonical, high-confidence
    shared       + novel          → novel, robustly detected by both methods
    only_short   + in_reference   → known peptide, short-read only
    only_short   + novel          → novel + short-read only (lower confidence)
    only_long    + in_reference   → known peptide, long-read only
    only_long    + novel          → novel + long-read only (interesting)

Usage:
    # Annotate a single TSV (one k value):
    python kmer_ref_annotate.py \\
        --input  kmer_output/kmers_k9.tsv \\
        --ref    reference_proteome.fa \\
        --outdir ./annotated

    # Annotate multiple TSVs at once:
    python kmer_ref_annotate.py \\
        --input  kmer_output/kmers_k9.tsv kmer_output/kmers_k10.tsv \\
        --ref    reference_proteome.fa \\
        --outdir ./annotated

    # Also write per-category split TSVs:
    python kmer_ref_annotate.py \\
        --input  kmer_output/kmers_k9.tsv \\
        --ref    reference_proteome.fa \\
        --outdir ./annotated \\
        --split-categories

    # Track which reference protein(s) each k-mer matched:
    python kmer_ref_annotate.py \\
        --input  kmer_output/kmers_k9.tsv \\
        --ref    reference_proteome.fa \\
        --outdir ./annotated \\
        --track-ref-proteins

Output:
    annotated/kmers_k{k}_annotated.tsv
    Columns (appended to original columns):
        ref_match          : True / False
        ref_proteins       : semicolon-joined reference protein IDs (if --track-ref-proteins)
        category           : one of the 6 categories above
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AMBIGUOUS_AA = set("XBZUO*")

CATEGORY_ORDER = [
    "shared|in_reference",
    "shared|novel",
    "only_short|in_reference",
    "only_short|novel",
    "only_long|in_reference",
    "only_long|novel",
]


# ---------------------------------------------------------------------------
# Step 1: Build reference k-mer index
# ---------------------------------------------------------------------------

def build_ref_index(
    ref_fasta: str,
    k: int,
    track_proteins: bool = False,
    skip_ambiguous: bool = True,
) -> dict:
    """
    Extract all k-mers from the reference proteome.

    Args:
        ref_fasta       : path to reference FASTA
        k               : k-mer size
        track_proteins  : if True, map kmer → set of protein IDs
                          if False, return a plain set (faster, less memory)
        skip_ambiguous  : skip k-mers containing ambiguous amino acids

    Returns:
        ref_kmers : set of k-mer strings  (track_proteins=False)
                  | dict {kmer: set(protein_ids)}  (track_proteins=True)
    """
    print(f"\n[Reference] Building k={k} index from {ref_fasta} ...")

    if track_proteins:
        ref_kmers: dict = defaultdict(set)
    else:
        ref_kmers: set = set()

    n_seqs = 0
    skipped_short = 0
    skipped_ambig = 0

    for rec in tqdm(SeqIO.parse(ref_fasta, "fasta"), desc="  [ref] scanning", unit="seq"):
        seq = str(rec.seq).upper().strip("*")
        if not seq:
            continue
        n_seqs += 1
        if len(seq) < k:
            skipped_short += 1
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if skip_ambiguous and AMBIGUOUS_AA.intersection(kmer):
                skipped_ambig += 1
                continue
            if track_proteins:
                ref_kmers[kmer].add(rec.id)
            else:
                ref_kmers.add(kmer)

    n_unique = len(ref_kmers)
    print(f"  [ref] {n_seqs:,} sequences → {n_unique:,} unique k-mers (k={k})")
    if skipped_short:
        print(f"  [ref] {skipped_short:,} sequences shorter than k skipped")
    if skipped_ambig:
        print(f"  [ref] {skipped_ambig:,} ambiguous k-mers skipped")

    return ref_kmers


# ---------------------------------------------------------------------------
# Step 2: Annotate pile TSV with reference match
# ---------------------------------------------------------------------------

def annotate_with_reference(
    df: pd.DataFrame,
    ref_kmers,
    track_proteins: bool = False,
) -> pd.DataFrame:
    """
    Add ref_match (bool), optionally ref_proteins (str), and category columns.

    Args:
        df            : DataFrame from kmer_pile_assign output
        ref_kmers     : set or dict from build_ref_index
        track_proteins: whether ref_kmers is a dict with protein sets

    Returns:
        Annotated DataFrame
    """
    print("  Annotating k-mers against reference ...")

    if track_proteins:
        ref_match   = df["kmer"].apply(lambda km: km in ref_kmers)
        ref_proteins = df["kmer"].apply(
            lambda km: ";".join(sorted(ref_kmers[km])) if km in ref_kmers else ""
        )
    else:
        ref_match = df["kmer"].isin(ref_kmers)

    df = df.copy()
    df["ref_match"] = ref_match

    if track_proteins:
        df["ref_proteins"] = ref_proteins

    # Build category label: "{pile}|{ref_status}"
    ref_label = df["ref_match"].map({True: "in_reference", False: "novel"})
    df["category"] = df["pile"] + "|" + ref_label

    return df


def print_summary(df: pd.DataFrame, k: int) -> None:
    """Print a formatted summary table of the 6-category breakdown."""
    total = len(df)
    counts = df["category"].value_counts()

    print(f"\n  ── k={k} novelty summary {'─'*35}")
    print(f"  {'Category':<35} {'Count':>10}  {'%':>6}")
    print(f"  {'─'*35} {'─'*10}  {'─'*6}")

    for cat in CATEGORY_ORDER:
        n   = counts.get(cat, 0)
        pct = 100 * n / total if total else 0
        marker = "  ◀ novel" if "novel" in cat else ""
        print(f"  {cat:<35} {n:>10,}  {pct:>5.1f}%{marker}")

    print(f"  {'─'*35} {'─'*10}")
    print(f"  {'TOTAL':<35} {total:>10,}")

    # Quick summary: how many kmers are novel overall
    n_novel = df["ref_match"].eq(False).sum()
    pct_novel = 100 * n_novel / total if total else 0
    print(f"\n  Overall novel (no reference match): {n_novel:,} / {total:,} ({pct_novel:.1f}%)")
    print()


# ---------------------------------------------------------------------------
# Step 3: Write output
# ---------------------------------------------------------------------------

def write_outputs(
    df: pd.DataFrame,
    outdir: str,
    k: int,
    split_categories: bool = False,
) -> None:
    """
    Write the annotated TSV and optionally per-category split TSVs.
    """
    os.makedirs(outdir, exist_ok=True)

    # Main annotated TSV
    main_path = os.path.join(outdir, f"kmers_k{k}_annotated.tsv")
    df.to_csv(main_path, sep="\t", index=False)
    print(f"  Saved → {main_path}")

    # Optional per-category split
    if split_categories:
        split_dir = os.path.join(outdir, f"k{k}_by_category")
        os.makedirs(split_dir, exist_ok=True)
        for cat in df["category"].unique():
            subset = df[df["category"] == cat]
            safe_name = cat.replace("|", "_")
            cat_path = os.path.join(split_dir, f"{safe_name}.tsv")
            subset.to_csv(cat_path, sep="\t", index=False)
            print(f"    [{cat}] {len(subset):,} rows → {cat_path}")


# ---------------------------------------------------------------------------
# Utility: infer k from TSV filename  (e.g. "kmers_k9.tsv" → 9)
# ---------------------------------------------------------------------------

def infer_k_from_filename(path: str) -> int | None:
    m = re.search(r"_k(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate kmer_pile_assign output TSVs with reference proteome match status.\n"
            "Produces a 6-category classification: pile × ref_match."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", nargs="+", required=True, metavar="TSV",
        help="One or more kmer_pile_assign output TSVs (e.g. kmers_k9.tsv kmers_k10.tsv)"
    )
    parser.add_argument(
        "--ref", required=True, metavar="FASTA",
        help="Reference proteome FASTA file"
    )
    parser.add_argument(
        "--outdir", default=".", metavar="DIR",
        help="Output directory (created if absent). Default: current directory"
    )
    parser.add_argument(
        "--track-ref-proteins", action="store_true",
        help=(
            "Add a 'ref_proteins' column listing which reference protein(s) "
            "each k-mer matched. Uses more memory. Default: off."
        )
    )
    parser.add_argument(
        "--split-categories", action="store_true",
        help="Also write one TSV per category into a subdirectory. Default: off."
    )
    parser.add_argument(
        "--keep-ambiguous", action="store_true",
        help="Include ambiguous amino acid k-mers when building the reference index."
    )
    args = parser.parse_args()

    # Validate inputs
    for path in args.input:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: Input TSV not found: {path}")
    if not os.path.isfile(args.ref):
        sys.exit(f"ERROR: Reference FASTA not found: {args.ref}")

    skip_amb = not args.keep_ambiguous

    # Collect unique k values needed (to avoid re-indexing reference per file)
    k_values = {}
    for path in args.input:
        k = infer_k_from_filename(path)
        if k is None:
            sys.exit(
                f"ERROR: Cannot infer k from filename '{os.path.basename(path)}'.\n"
                f"  Expected pattern: kmers_k{{k}}.tsv  (e.g. kmers_k9.tsv)"
            )
        k_values[path] = k

    # Group input paths by k so we only build each reference index once
    k_to_paths: dict[int, list[str]] = defaultdict(list)
    for path, k in k_values.items():
        k_to_paths[k].append(path)

    # Process each unique k
    for k, paths in sorted(k_to_paths.items()):
        print(f"\n{'='*60}")
        print(f"  k = {k}")
        print(f"{'='*60}")

        # Step 1: build reference index for this k
        ref_kmers = build_ref_index(
            args.ref, k,
            track_proteins=args.track_ref_proteins,
            skip_ambiguous=skip_amb,
        )

        # Step 2 & 3: annotate each TSV for this k
        for path in paths:
            print(f"\n  Input: {path}")
            df = pd.read_csv(path, sep="\t", dtype=str)

            # Coerce count columns back to int if present
            for col in df.columns:
                if col.endswith("_count"):
                    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

            df = annotate_with_reference(df, ref_kmers, track_proteins=args.track_ref_proteins)
            print_summary(df, k)

            # Step 3: write outputs
            write_outputs(df, args.outdir, k, split_categories=args.split_categories)

    print("\nDone.")


if __name__ == "__main__":
    main()
