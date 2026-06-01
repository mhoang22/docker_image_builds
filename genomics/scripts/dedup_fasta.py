#!/usr/bin/env python3
"""
dedup_fasta.py

Deduplicate a FASTA file by:
  1. Collapsing identical sequences      -> keep the first-encountered entry as representative
  2. Removing entries whose sequence is a substring of a longer entry

Uses Aho-Corasick multi-pattern search for speed (handles ~40k sequences in seconds).

Outputs:
  <prefix>_dedup.fa      -- deduplicated FASTA; header = representative ID
  <prefix>_clusters.tsv  -- provenance table:
        representative_id | related_id | relationship
        relationship is either 'identical' or 'substring'

Usage:
    python dedup_fasta.py input.fa
    python dedup_fasta.py input.fa --output_dir outdir_path --output_prefix my_output

Requirements:
    pip install pyahocorasick --break-system-packages
"""

import argparse
import sys
from pathlib import Path

try:
    import ahocorasick
except ImportError:
    sys.exit(
        "Error: pyahocorasick not found.\n"
        "Install with: pip install pyahocorasick --break-system-packages"
    )


# --------------------------------------------------------------------------- #
# I/O helpers
# --------------------------------------------------------------------------- #

def parse_fasta(filepath):
    """Return list of (name, seq) tuples; sequence is uppercased."""
    entries = []
    name, seq_parts = None, []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]   # ID = everything before first space
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            entries.append((name, "".join(seq_parts).upper()))
    return entries


def write_fasta(entries, filepath):
    """Write (name, seq) pairs to a FASTA file, sequences wrapped at 60 chars."""
    with open(filepath, "w") as f:
        for name, seq in entries:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def write_tsv(provenance, filepath):
    """Write provenance records to a TSV file."""
    with open(filepath, "w") as f:
        f.write("representative_id\trelated_id\trelationship\n")
        for rep, related, rel in provenance:
            f.write(f"{rep}\t{related}\t{rel}\n")


# --------------------------------------------------------------------------- #
# Core deduplication
# --------------------------------------------------------------------------- #

def collapse_exact_duplicates(entries):
    """
    Group entries with identical sequences.

    Returns:
        unique      : list of (rep_name, seq) -- first-encountered name per sequence,
                      preserving original entry order
        provenance  : list of (rep_name, duplicate_name, 'identical')
    """
    seq_to_names = {}
    order = []   # first-seen order of sequences

    for name, seq in entries:
        if seq not in seq_to_names:
            seq_to_names[seq] = []
            order.append(seq)
        seq_to_names[seq].append(name)

    unique = []
    provenance = []
    for seq in order:
        names = seq_to_names[seq]
        rep = names[0]
        unique.append((rep, seq))
        for other in names[1:]:
            provenance.append((rep, other, "identical"))

    return unique, provenance


def find_substrings(unique):
    """
    Given a list of (name, seq) with unique sequences, find all entries whose
    sequence is a substring of some longer entry's sequence.

    Uses Aho-Corasick: build an automaton of ALL sequences, then scan each
    sequence -- any pattern match by a *different, shorter* sequence flags
    a substring relationship.

    Returns:
        kept       : list of (name, seq) that are NOT substrings of anything longer,
                     in original entry order
        provenance : list of (absorber_name, absorbed_name, 'substring')
    """
    # Sort longest-first so that when we scan seq[i] and find seq[j] inside it,
    # j > i (j is shorter or equal length but later), meaning it is the substring.
    indexed = list(enumerate(unique))
    indexed.sort(key=lambda x: -len(x[1][1]))   # descending seq length

    sorted_seqs = [(name, seq) for _, (name, seq) in indexed]

    print("  Building Aho-Corasick automaton ...", flush=True)
    A = ahocorasick.Automaton()
    for sort_idx, (name, seq) in enumerate(sorted_seqs):
        A.add_word(seq, (sort_idx, name))
    A.make_automaton()

    print("  Scanning for substrings ...", flush=True)
    # absorbed_by[sort_idx of absorbed] = sort_idx of absorber (longer seq)
    absorbed_by = {}

    for sort_idx, (name, seq) in enumerate(sorted_seqs):
        for _, (match_sort_idx, match_name) in A.iter(seq):
            if match_sort_idx == sort_idx:
                continue                      # skip self-match
            if match_sort_idx in absorbed_by:
                continue                      # already recorded
            # match_sort_idx sequence was found inside current seq.
            # Since sorted longest-first, match_sort_idx > sort_idx means
            # the matched sequence is strictly shorter (or equal-length but later).
            if match_sort_idx > sort_idx:
                absorbed_by[match_sort_idx] = sort_idx

    kept = []
    provenance = []
    for sort_idx, (name, seq) in enumerate(sorted_seqs):
        if sort_idx in absorbed_by:
            absorber_name = sorted_seqs[absorbed_by[sort_idx]][0]
            provenance.append((absorber_name, name, "substring"))
        else:
            kept.append((name, seq))

    # Restore original entry order for stable, reproducible output
    kept_names_set = {name for name, _ in kept}
    kept_ordered = [(name, seq) for name, seq in unique if name in kept_names_set]

    return kept_ordered, provenance


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument(
        "--output_prefix",
        default=None,
        help="Prefix for output files (default: input file stem)",
    )
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Directory for output files (default: current working directory)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_fasta)
    if not input_path.exists():
        sys.exit(f"Error: file not found: {input_path}")

    prefix   = args.output_prefix or input_path.stem
    out_dir  = Path(args.output_dir) if args.output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fa   = out_dir / (prefix + "_dedup.fa")
    out_tsv  = out_dir / (prefix + "_clusters.tsv")

    # ---- Parse ----
    print(f"Reading {input_path} ...")
    entries = parse_fasta(input_path)
    print(f"  Loaded {len(entries):,} entries")

    # ---- Pass 1: exact duplicates ----
    print("Pass 1: collapsing exact duplicates ...")
    unique, prov_identical = collapse_exact_duplicates(entries)
    print(f"  {len(entries) - len(unique):,} duplicates collapsed -> {len(unique):,} unique sequences")

    # ---- Pass 2: substrings ----
    print("Pass 2: finding substring relationships ...")
    kept, prov_substring = find_substrings(unique)
    print(f"  {len(unique) - len(kept):,} substrings removed -> {len(kept):,} representatives kept")

    # ---- Write outputs ----
    all_provenance = prov_identical + prov_substring

    print(f"Writing {out_fa} ...")
    write_fasta(kept, out_fa)

    print(f"Writing {out_tsv} ...")
    write_tsv(all_provenance, out_tsv)

    # ---- Summary ----
    print("\n-- Summary -----------------------------------------------")
    print(f"  Input entries        : {len(entries):,}")
    print(f"  Identical collapsed  : {len(prov_identical):,}")
    print(f"  Substrings removed   : {len(prov_substring):,}")
    print(f"  Representatives kept : {len(kept):,}")
    print(f"  Output FASTA         : {out_fa}")
    print(f"  Cluster TSV          : {out_tsv}")
    print("----------------------------------------------------------")


if __name__ == "__main__":
    main()