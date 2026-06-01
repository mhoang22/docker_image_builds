#!/usr/bin/env bash
# Usage: filter_peptides.sh <input_file> <output_file>
# Keeps unique peptide sequences with 8+ amino acids (length > 7)

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_file>" >&2
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input file '$INPUT' not found." >&2
    exit 1
fi

# sort | uniq | keep only lines where length > 7
sort -u "$INPUT" | awk 'length($0) > 7' > "$OUTPUT"

echo "Done: $(wc -l < "$OUTPUT") peptides written to $OUTPUT"