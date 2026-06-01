"""
compare_peptide_lists_DatabaseVsDenovo.py
========================
Compares peptide lists from database search tools (e.g. FragPipe) against
de novo sequencing tools (e.g. InstaNovo, Cascadia) for immunopeptidomics data.

USAGE
-----
    python compare_peptide_lists_DatabaseVsDenovo.py \\
        --db      fragpipe_run3.txt fragpipe_run8.txt \\
        --denovo  instanovo.txt cascadia.txt \\
        --output  results/

    Each file should be a plain text file with one peptide sequence per line.
    File names (without extension) are used as labels throughout the output.

MATCHING STRATEGY
-----------------
Two levels of matching are performed:

  1. EXACT MATCH (after I/L normalization)
     All peptides are normalized by replacing isoleucine (I) with leucine (L)
     before comparison, because de novo sequencing tools cannot distinguish
     between I and L (they have identical masses). Matches are then reported
     in terms of the *original* sequences.

  2. SUBSTRING MATCH (de novo peptide as substring of database peptide only)

     *** IMPORTANT DISCLAIMER ***
     Substring matching is performed in ONE direction only:
       - We check if a de novo peptide is contained WITHIN a database peptide.
       - We do NOT check if a database peptide is contained within a de novo peptide.

     Rationale: De novo tools occasionally miss one or two residues at the
     N- or C-terminus of a peptide, producing a truncated version of the true
     sequence. Checking whether a de novo peptide is a substring of a database
     peptide catches these near-miss cases.

     The reverse direction (database peptide as substring of de novo peptide)
     is NOT checked because database search tools like FragPipe use strict
     enzymatic cleavage rules and length constraints — so a database peptide
     being "shorter than" a de novo prediction would not represent the same
     underlying peptide. Including that direction would inflate overlap counts
     with biologically meaningless matches.

     Substring matches that are also exact matches are excluded from the
     substring-only count to avoid double-counting.
     *** END DISCLAIMER ***

OUTPUT FILES
------------
  overlap_summary.tsv       — pairwise counts table (exact + substring)
  exact_overlaps/           — per-pair TSV of exact-matching peptides
  substring_overlaps/       — per-pair TSV of substring-only matching peptides
  all_files_overlap.txt     — peptides present in ALL files (exact, I/L-norm)
  venn_counts.tsv           — counts for every subset of the power set (for Venn)
"""

import argparse
import itertools
import os
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_peptides(filepath: str) -> list[str]:
    """Load peptides from a one-peptide-per-line file. Skips blank lines."""
    with open(filepath, "r") as fh:
        peptides = [line.strip() for line in fh if line.strip()]
    return peptides


def normalize_il(peptide: str) -> str:
    """Replace all I with L for I/L-agnostic comparison."""
    return peptide.replace("I", "L")


def build_normalized_map(peptides: list[str]) -> dict[str, list[str]]:
    """
    Returns a dict mapping normalized sequence -> list of original sequences.
    One-to-many because e.g. ALSGHLETL and ALSGHLETL (already L) both normalize
    to the same key. In practice there will usually be just one original per key.
    """
    mapping: dict[str, list[str]] = defaultdict(list)
    for pep in peptides:
        mapping[normalize_il(pep)].append(pep)
    return mapping


# ---------------------------------------------------------------------------
# Matching logic
# ---------------------------------------------------------------------------

def exact_overlap(
    norm_map_a: dict[str, list[str]],
    norm_map_b: dict[str, list[str]],
) -> list[tuple[str, str]]:
    """
    Return list of (original_from_a, original_from_b) pairs that match
    exactly after I/L normalization.
    """
    shared_keys = set(norm_map_a.keys()) & set(norm_map_b.keys())
    pairs = []
    for key in shared_keys:
        for orig_a in norm_map_a[key]:
            for orig_b in norm_map_b[key]:
                pairs.append((orig_a, orig_b))
    return pairs


def substring_overlap(
    denovo_peptides: list[str],
    db_norm_map: dict[str, list[str]],
    exact_norm_keys: set[str],
) -> list[tuple[str, str]]:
    """
    Check whether each de novo peptide (after I/L normalization) is a substring
    of any database peptide (also normalized). Excludes pairs already captured
    by exact matching.

    Returns list of (denovo_original, db_original) pairs.
    """
    # Build a list of (norm_db_peptide, orig_db_peptide) for substring search
    db_norm_list: list[tuple[str, str]] = []
    for norm_seq, originals in db_norm_map.items():
        for orig in originals:
            db_norm_list.append((norm_seq, orig))

    matches = []
    for dn_pep in denovo_peptides:
        dn_norm = normalize_il(dn_pep)
        # Skip if this is already an exact match
        if dn_norm in exact_norm_keys:
            continue
        for db_norm, db_orig in db_norm_list:
            # De novo (normalized) is a proper substring of a db peptide (normalized)
            if dn_norm in db_norm and dn_norm != db_norm:
                matches.append((dn_pep, db_orig))
                break  # one match per de novo peptide is enough
    return matches


# ---------------------------------------------------------------------------
# All-files overlap
# ---------------------------------------------------------------------------

def all_files_exact_overlap(all_norm_maps: dict[str, dict[str, list[str]]]) -> set[str]:
    """
    Return the set of I/L-normalized sequences present in every file.
    """
    norm_key_sets = [set(m.keys()) for m in all_norm_maps.values()]
    return set.intersection(*norm_key_sets)


# ---------------------------------------------------------------------------
# Power-set Venn counts
# ---------------------------------------------------------------------------

def venn_counts(
    all_norm_maps: dict[str, dict[str, list[str]]],
) -> dict[frozenset[str], int]:
    """
    For every non-empty subset of file labels, count how many I/L-normalized
    peptides appear in exactly that subset (exclusive membership).
    Used for Venn / UpSet plot construction.
    """
    labels = list(all_norm_maps.keys())
    # Build a dict: norm_seq -> set of labels it appears in
    seq_to_labels: dict[str, set[str]] = defaultdict(set)
    for label, norm_map in all_norm_maps.items():
        for norm_seq in norm_map:
            seq_to_labels[norm_seq].add(label)

    counts: dict[frozenset[str], int] = defaultdict(int)
    for norm_seq, label_set in seq_to_labels.items():
        counts[frozenset(label_set)] += 1
    return counts


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_pairwise_tsv(filepath: str, header: list[str], rows: list[tuple]) -> None:
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def write_summary(
    output_dir: str,
    db_labels: list[str],
    denovo_labels: list[str],
    exact_counts: dict[tuple[str, str], int],
    substring_counts: dict[tuple[str, str], int],
) -> None:
    filepath = os.path.join(output_dir, "overlap_summary.tsv")
    header = ["db_file", "denovo_file", "exact_overlap", "substring_only_overlap",
              "total_overlap"]
    rows = []
    for db_label in db_labels:
        for dn_label in denovo_labels:
            key = (db_label, dn_label)
            exact = exact_counts.get(key, 0)
            substr = substring_counts.get(key, 0)
            rows.append((db_label, dn_label, exact, substr, exact + substr))
    write_pairwise_tsv(filepath, header, rows)
    print(f"  Written: {filepath}")


def write_all_files_overlap(
    output_dir: str,
    shared_norm_keys: set[str],
    all_norm_maps: dict[str, dict[str, list[str]]],
) -> None:
    filepath = os.path.join(output_dir, "all_files_overlap.txt")
    with open(filepath, "w") as fh:
        fh.write("# Peptides (I/L-normalized) present in ALL input files\n")
        fh.write("# Original sequences from each file shown alongside\n")
        labels = list(all_norm_maps.keys())
        fh.write("normalized_sequence\t" + "\t".join(labels) + "\n")
        for norm_seq in sorted(shared_norm_keys):
            originals = [
                ";".join(all_norm_maps[label].get(norm_seq, ["-"]))
                for label in labels
            ]
            fh.write(norm_seq + "\t" + "\t".join(originals) + "\n")
    print(f"  Written: {filepath}")


def write_venn_counts(output_dir: str, counts: dict[frozenset[str], int]) -> None:
    filepath = os.path.join(output_dir, "venn_counts.tsv")
    with open(filepath, "w") as fh:
        fh.write("file_combination\tpeptide_count\n")
        for subset, count in sorted(counts.items(), key=lambda x: (-len(x[0]), sorted(x[0]))):
            combo = " & ".join(sorted(subset))
            fh.write(f"{combo}\t{count}\n")
    print(f"  Written: {filepath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare peptide lists from database search vs de novo tools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--db", nargs="+", required=True, metavar="FILE",
        help="One or more database-search peptide list files (e.g. FragPipe output).",
    )
    parser.add_argument(
        "--denovo", nargs="+", required=True, metavar="FILE",
        help="One or more de novo peptide list files (e.g. InstaNovo, Cascadia output).",
    )
    parser.add_argument(
        "--output", default="overlap_results", metavar="DIR",
        help="Output directory (default: overlap_results/).",
    )
    args = parser.parse_args()

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    exact_out_dir = os.path.join(output_dir, "exact_overlaps")
    substr_out_dir = os.path.join(output_dir, "substring_overlaps")

    # --- Load all files ---
    print("\n=== Loading peptide files ===")

    db_data: dict[str, list[str]] = {}
    db_norm_maps: dict[str, dict[str, list[str]]] = {}
    for f in args.db:
        label = Path(f).stem
        peps = load_peptides(f)
        db_data[label] = peps
        db_norm_maps[label] = build_normalized_map(peps)
        print(f"  [DB]     {label}: {len(peps):,} peptides")

    denovo_data: dict[str, list[str]] = {}
    denovo_norm_maps: dict[str, dict[str, list[str]]] = {}
    for f in args.denovo:
        label = Path(f).stem
        peps = load_peptides(f)
        denovo_data[label] = peps
        denovo_norm_maps[label] = build_normalized_map(peps)
        print(f"  [DeNovo] {label}: {len(peps):,} peptides")

    # --- Pairwise comparisons ---
    print("\n=== Pairwise comparisons ===")
    exact_counts: dict[tuple[str, str], int] = {}
    substring_counts: dict[tuple[str, str], int] = {}

    for db_label, dn_label in itertools.product(db_norm_maps.keys(), denovo_norm_maps.keys()):
        db_nmap = db_norm_maps[db_label]
        dn_nmap = denovo_norm_maps[dn_label]
        dn_peps = denovo_data[dn_label]
        pair_key = (db_label, dn_label)

        # Exact matches
        exact_pairs = exact_overlap(db_nmap, dn_nmap)
        exact_norm_keys = set(normalize_il(p[0]) for p in exact_pairs)
        exact_counts[pair_key] = len(exact_pairs)

        print(f"  {db_label} vs {dn_label}:")
        print(f"    Exact matches:     {len(exact_pairs):,}")

        # Write exact overlap file
        if exact_pairs:
            out_path = os.path.join(exact_out_dir, f"{db_label}__{dn_label}.tsv")
            write_pairwise_tsv(
                out_path,
                [f"{db_label}_original", f"{dn_label}_original", "il_normalized"],
                [(a, b, normalize_il(a)) for a, b in exact_pairs],
            )

        # Substring matches (de novo as substring of db only)
        substr_pairs = substring_overlap(dn_peps, db_nmap, exact_norm_keys)
        substring_counts[pair_key] = len(substr_pairs)
        print(f"    Substring-only:    {len(substr_pairs):,}")
        print(f"    Total overlap:     {len(exact_pairs) + len(substr_pairs):,}")

        # Write substring overlap file
        if substr_pairs:
            out_path = os.path.join(substr_out_dir, f"{db_label}__{dn_label}.tsv")
            write_pairwise_tsv(
                out_path,
                [f"{dn_label}_denovo", f"{db_label}_db_containing_peptide"],
                substr_pairs,
            )

    # --- All-files overlap ---
    print("\n=== All-files overlap ===")
    all_norm_maps = {**db_norm_maps, **denovo_norm_maps}
    shared = all_files_exact_overlap(all_norm_maps)
    print(f"  Peptides in ALL {len(all_norm_maps)} files: {len(shared):,}")
    write_all_files_overlap(output_dir, shared, all_norm_maps)

    # --- Venn / UpSet counts ---
    print("\n=== Venn subset counts ===")
    vcounts = venn_counts(all_norm_maps)
    write_venn_counts(output_dir, vcounts)
    for subset, count in sorted(vcounts.items(), key=lambda x: (-len(x[0]), sorted(x[0]))):
        combo = " & ".join(sorted(subset))
        print(f"  {combo}: {count:,}")

    # --- Summary table ---
    print("\n=== Writing summary ===")
    write_summary(output_dir, list(db_norm_maps.keys()), list(denovo_norm_maps.keys()),
                  exact_counts, substring_counts)

    print("\nDone. Results written to:", output_dir)


if __name__ == "__main__":
    main()
