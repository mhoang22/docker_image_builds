#!/usr/bin/env python3
"""
Convert SQANTI3 QC output to protein FASTA.
Usage: python sqanti3_to_fasta.py <input_txt> <output_fasta>
"""

import sys
import csv

def sqanti3_to_fasta(input_file, output_file):
    written = 0
    skipped = 0

    with open(input_file, "r", newline="") as infile, open(output_file, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            isoform = row.get("isoform", "").strip()
            orf_seq = row.get("ORF_seq", "").strip()

            # Skip if no protein sequence or sequence is NA/empty
            if not orf_seq or orf_seq == "NA":
                skipped += 1
                continue

            # Remove trailing stop codon asterisk if present
            orf_seq = orf_seq.rstrip("*")

            outfile.write(f">{isoform}\n{orf_seq}\n")
            written += 1

    print(f"Done. Wrote {written} entries, skipped {skipped} (no protein sequence).")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sqanti3_to_fasta.py <input_txt> <output_fasta>")
        sys.exit(1)

    sqanti3_to_fasta(sys.argv[1], sys.argv[2])