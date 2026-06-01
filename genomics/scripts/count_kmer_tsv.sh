#!/usr/bin/env bash
# Usage: bash count_kmer_tsv.sh <path_to_kN_by_category_dir>
# Example: bash count_kmer_tsv.sh result_ref_annotate/k8_by_category
#
# Counts data rows (total rows minus 1 header) for each of the 6 TSV files
# and prints a two-column table: kmer_type | count

DIR="${1:-.}"   # default to current directory if no argument given

TSV_FILES=(
    "only_long_in_reference.tsv"
    "only_long_novel.tsv"
    "only_short_in_reference.tsv"
    "only_short_novel.tsv"
    "shared_in_reference.tsv"
    "shared_novel.tsv"
)

printf "%-35s %s\n" "kmer type" "count"
printf "%-35s %s\n" "$(printf '%.0s-' {1..35})" "-----"

for tsv in "${TSV_FILES[@]}"; do
    filepath="${DIR}/${tsv}"
    kmer_type="${tsv%.tsv}"   # strip .tsv suffix

    if [[ -f "$filepath" ]]; then
        total_lines=$(wc -l < "$filepath")
        count=$(( total_lines - 1 ))
    else
        count="FILE NOT FOUND"
    fi

    printf "%-35s %s\n" "$kmer_type" "$count"
done
