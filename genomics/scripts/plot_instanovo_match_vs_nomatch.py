"""
Plot distribution of log_probs and delta_mass_ppm for InstaNovo predictions,
comparing peptides that MATCH vs NO_MATCH against a custom reference file.

Usage:
    python plot_instanovo_match_vs_nomatch.py \
        --instanovo  <instanovo_output.csv> \
        --reference  <custom_file.csv/tsv> \
        --out        <output_plot.png>          # optional, default: output_plot.png

File expectations:
    InstaNovo file : CSV with columns including 'instanovoplus_predictions',
                     'instanovoplus_prediction_log_probability', 'delta_mass_ppm'
    Reference file : CSV/TSV with columns 'peptide', 'n_proteins', 'protein_ids'
                     (NO_MATCH when n_proteins == 0)

Peptide matching logic:
    instanovoplus_predictions is a stringified token list e.g. ['G','S[UNIMOD:21]','R']
    → tokens are joined and mods stripped → 'GSR'
    This matches the 'peptide' column in the reference file (already stripped).
"""

import argparse
import sys
import re
import ast
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── helpers ──────────────────────────────────────────────────────────────────

def load_table(path: str) -> pd.DataFrame:
    """Load CSV or TSV, auto-detecting separator."""
    p = Path(path)
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    return pd.read_csv(p, sep=sep, low_memory=False)


def parse_instanovoplus(value) -> str:
    """
    Parse instanovoplus_predictions token list into a clean bare amino acid sequence.
    Handles:
      - stringified list:  "['G', 'S[UNIMOD:21]', 'R']"  →  'GSR'
      - plain string:      'GSR'                           →  'GSR'
      - already joined:    'GS[UNIMOD:21]R'               →  'GSR'
    """
    s = str(value).strip()
    if s.startswith("["):
        try:
            tokens = ast.literal_eval(s)
            joined = "".join(str(t) for t in tokens)
        except Exception:
            joined = s
    else:
        joined = s
    # strip any bracket-enclosed modifications
    cleaned = re.sub(r'\[.*?\]', '', joined).strip().upper()
    return cleaned


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Box-plot log probability & delta_mass_ppm: MATCH vs NO_MATCH"
    )
    parser.add_argument("--instanovo", required=True,
                        help="InstaNovo output CSV (must contain "
                             "'instanovoplus_predictions', "
                             "'instanovoplus_prediction_log_probability', "
                             "'delta_mass_ppm')")
    parser.add_argument("--reference", required=True,
                        help="Reference file with columns 'peptide', 'n_proteins'")
    parser.add_argument("--out", default="output_plot.png",
                        help="Output figure path (default: output_plot.png)")
    args = parser.parse_args()

    # ── 1. load ──────────────────────────────────────────────────────────────
    print(f"Loading InstaNovo file : {args.instanovo}")
    insta = load_table(args.instanovo)

    print(f"Loading reference file : {args.reference}")
    ref   = load_table(args.reference)

    # ── 2. validate required columns ─────────────────────────────────────────
    required_insta = {"instanovoplus_predictions",
                      "instanovoplus_prediction_log_probability",
                      "delta_mass_ppm"}
    required_ref   = {"peptide", "n_proteins"}

    missing_insta = required_insta - set(insta.columns)
    missing_ref   = required_ref   - set(ref.columns)

    if missing_insta:
        sys.exit(f"ERROR – InstaNovo file is missing columns: {missing_insta}\n"
                 f"Found columns: {list(insta.columns)}")
    if missing_ref:
        sys.exit(f"ERROR – Reference file is missing columns: {missing_ref}\n"
                 f"Found columns: {list(ref.columns)}")

    # ── 3. build match/no-match lookup from reference ─────────────────────────
    ref["_peptide_key"] = ref["peptide"].astype(str).str.strip().str.upper()
    ref["_is_match"]    = ref["n_proteins"].astype(int) > 0
    match_lookup = dict(zip(ref["_peptide_key"], ref["_is_match"]))

    # ── 4. parse instanovoplus sequences and annotate ─────────────────────────
    print("Parsing instanovoplus_predictions column…")
    insta["_pred_key"] = insta["instanovoplus_predictions"].apply(parse_instanovoplus)
    insta["_match_status"] = insta["_pred_key"].map(match_lookup)

    n_total   = len(insta)
    n_found   = insta["_match_status"].notna().sum()
    n_match   = (insta["_match_status"] == True).sum()
    n_nomatch = (insta["_match_status"] == False).sum()
    n_absent  = n_total - n_found

    print(f"\nAnnotation summary")
    print(f"  Total InstaNovo predictions : {n_total}")
    print(f"  Found in reference          : {n_found}")
    print(f"    → MATCH  (n_proteins > 0) : {n_match}")
    print(f"    → NO_MATCH                : {n_nomatch}")
    print(f"  Not found in reference      : {n_absent}  "
          f"(peptide not present in reference file)")

    if n_found == 0:
        # show a few parsed keys to help debug
        print("\nSample parsed keys from InstaNovo:")
        for k in insta["_pred_key"].head(5):
            print(f"  '{k}'")
        print("\nSample keys from reference:")
        for k in list(match_lookup.keys())[:5]:
            print(f"  '{k}'")
        sys.exit(
            "\nERROR – No predictions matched any peptide in the reference file.\n"
            "Check that both files use the same peptide format."
        )

    # ── 5. prepare data groups ───────────────────────────────────────────────
    match_df   = insta[insta["_match_status"] == True]
    nomatch_df = insta[insta["_match_status"] == False]

    colors = ["#4C9BE8", "#E87B4C"]   # blue = match, orange = no-match

    metrics = [
        ("instanovoplus_prediction_log_probability",
         "InstaNovo+ Log Probability"),
        ("delta_mass_ppm",
         "Delta Mass (ppm)"),
    ]

    # ── 6. plot ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 7))
    fig.suptitle(
        "InstaNovo+ Predictions: MATCH vs NO_MATCH\n"
        "(peptide key = stripped instanovoplus_predictions)",
        fontsize=13, fontweight="bold", y=1.01
    )

    jitter_rng = np.random.default_rng(42)

    for ax, (col, ylabel) in zip(axes, metrics):
        data_match   = match_df[col].dropna().astype(float)
        data_nomatch = nomatch_df[col].dropna().astype(float)
        data_arrays  = [data_match, data_nomatch]
        labels       = [f"MATCH\n(n={len(data_match)})",
                        f"NO_MATCH\n(n={len(data_nomatch)})"]

        bp = ax.boxplot(
            data_arrays,
            patch_artist=True,
            widths=0.45,
            medianprops=dict(color="black", linewidth=2.5),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            showfliers=False,
        )

        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        # jittered dots (subsample if large)
        for i, (vals, color) in enumerate(zip(data_arrays, colors), start=1):
            if len(vals) == 0:
                continue
            plot_vals = vals.sample(min(len(vals), 2000), random_state=42)
            x_jitter  = jitter_rng.uniform(-0.18, 0.18, size=len(plot_vals)) + i
            ax.scatter(x_jitter, plot_vals,
                       alpha=0.25, s=12, color=color, zorder=2, linewidths=0)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(labels, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(ylabel, fontsize=12, fontweight="bold")
        ax.yaxis.grid(True, linestyle="--", alpha=0.5)
        ax.set_axisbelow(True)

    # shared legend
    handles = [
        mpatches.Patch(facecolor=colors[0], alpha=0.7,
                       label=f"MATCH  (n_proteins > 0)"),
        mpatches.Patch(facecolor=colors[1], alpha=0.7,
                       label=f"NO_MATCH  (n_proteins = 0)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2,
               fontsize=11, bbox_to_anchor=(0.5, -0.04))

    plt.tight_layout()
    out_path = Path(args.out)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved → {out_path.resolve()}")
    plt.show()


if __name__ == "__main__":
    main()
