#!/usr/bin/env python3
"""
extract_cascadia_peptide.py

Extract cleaned peptide sequences from a Cascadia .ssl output file.
Modifications (e.g. C[+57.02]) are stripped, leaving bare amino acid sequences.
One peptide per line in the output file.

Usage:
    python extract_cascadia_peptide.py <input_ssl> <output_file>
"""

import sys
import re
import pandas as pd


def clean_sequence(seq: str) -> str:
    """Remove bracket-enclosed modifications from a peptide sequence."""
    return re.sub(r'\[.*?\]', '', seq).strip()


def extract_peptides(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path, sep='\t')

    if 'sequence' not in df.columns:
        raise ValueError(
            f"Expected a 'sequence' column in {input_path}, "
            f"but found: {df.columns.tolist()}"
        )

    cleaned = (
        df['sequence']
        .dropna()
        .apply(clean_sequence)
        .loc[lambda s: s.str.len() > 0]  # drop any empty strings
    )

    with open(output_path, 'w') as f:
        for peptide in cleaned:
            f.write(peptide + '\n')

    print(f"Done. {len(cleaned)} peptides written to {output_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python extract_cascadia_peptide.py <input_ssl> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_peptides(input_file, output_file)