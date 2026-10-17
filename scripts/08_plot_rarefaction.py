#!/usr/bin/env python3
"""Step 8 — Alpha rarefaction curves (Shannon) by group, from the data inside
the QIIME 2 alpha-rarefaction visualization exported in step 5.

Usage:
  python scripts/08_plot_rarefaction.py results/05_diversity/alpha_rarefaction_export/shannon.csv \
      --group cage_id --out results/07_phyloseq/figures/rarefaction_shannon_cage.png
"""
import argparse
import os
import re
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", type=Path)
    ap.add_argument("--group", default="cage_id")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    iters = [c for c in df.columns if re.match(r"depth-\d+_iter-\d+", c)]
    long = df.melt(id_vars=["sample-id", args.group], value_vars=iters, var_name="key", value_name="shannon")
    long["depth"] = long.key.str.extract(r"depth-(\d+)_").astype(int)
    per_sample = long.groupby(["sample-id", args.group, "depth"], as_index=False).shannon.mean()
    curves = per_sample.groupby([args.group, "depth"]).shannon.agg(["mean", "sem"]).reset_index()

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for g, d in curves.groupby(args.group):
        ax.errorbar(d.depth, d["mean"], yerr=d["sem"], marker="o", ms=3, capsize=2, label=str(g))
    ax.set(xlabel="Sequencing depth (reads)", ylabel="Shannon index (mean ± SEM)",
           title=f"Alpha rarefaction by {args.group}")
    ax.legend(title=args.group, frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)

    final = curves[curves.depth == curves.depth.max()].sort_values("mean", ascending=False)
    print(f"Shannon at the maximum depth ({curves.depth.max()} reads):")
    print(final[[args.group, "mean"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
