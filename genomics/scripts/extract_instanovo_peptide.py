#!/usr/bin/env python3
"""
extract_instanovo_peptide.py

Extract cleaned peptide sequences from an InstaNovo predictions CSV file.
Uses the 'instanovoplus_predictions' column (refined InstaNovo+ output).
Modifications (e.g. S[UNIMOD:21]) are stripped, leaving bare amino acid sequences.
One peptide per line in the output file.

Usage:
    python extract_instanovo_peptide.py <input_csv> <output_file>
"""

import sys
import re
import ast
import pandas as pd


def parse_token_list(value: str) -> list[str]:
    """Parse a stringified Python list of tokens, e.g. \"['G', 'S[UNIMOD:21]', 'R']\"."""
    return ast.literal_eval(value)


def clean_token(token: str) -> str:
    """Strip bracket-enclosed modifications from a single token, e.g. S[UNIMOD:21] -> S."""
    return re.sub(r'\[.*?\]', '', token)


def tokens_to_sequence(value: str) -> str:
    """Convert a stringified token list into a clean bare amino acid sequence."""
    tokens = parse_token_list(value)
    return ''.join(clean_token(t) for t in tokens)


def extract_peptides(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path)

    col = 'instanovoplus_predictions'
    if col not in df.columns:
        raise ValueError(
            f"Expected an '{col}' column in {input_path}, "
            f"but found: {df.columns.tolist()}"
        )

    cleaned = (
        df[col]
        .dropna()
        .apply(tokens_to_sequence)
        .loc[lambda s: s.str.len() > 0]  # drop any empty strings
    )

    with open(output_path, 'w') as f:
        for peptide in cleaned:
            f.write(peptide + '\n')

    print(f"Done. {len(cleaned)} peptides written to {output_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python extract_instanovo_peptide.py <input_csv> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_peptides(input_file, output_file)