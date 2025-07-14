#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
summarize_MD_comparisons.py  ▸  Aggregate & report MD comparison JSON files
===========================================================================
Reads any number of JSON summaries created by analysisMD.py and produces:

1. A tidy CSV table (one row per system) with key quantitative metrics.
2. An overall statistics block (mean / median / std of Complete–Published deltas).
3. Optional box‑plots illustrating cross‑system trends.

Typical use
-----------
    python summarize_MD_comparisons.py results/*_summary.json \
        --csv  all_systems_summary.csv \
        --plots summary_plots.pdf

Dependencies
------------
    python >= 3.8   •   numpy   •   pandas   •   matplotlib
"""

import argparse, json, sys, pathlib, statistics
from typing import Dict, Any, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ──────────────────────────────────────────────────────────────────────────────
def safe_get(d: Dict[str, Any], path: List[str], default=np.nan):
    """Traverse nested dict `d` following keys in `path` list; return default if any miss."""
    for key in path:
        if isinstance(d, dict) and key in d:
            d = d[key]
        else:
            return default
    return d


def extract_metrics(js: Dict[str, Any]) -> Dict[str, Any]:
    """Pull out the metrics we want to compare across systems.  
       Missing fields are filled with NaN so that pandas handles them gracefully.
    """
    out = {
        # Identifiers
        "system": pathlib.Path(js.get("out_prefix", "unknown")).stem,
        "label1": js.get("label1", "Structure1"),
        "label2": js.get("label2", "Structure2"),

        # Basic counts
        "n_residues_1": safe_get(js, ["n_residues", js.get("label1", "Complete Structure")]),
        "n_residues_2": safe_get(js, ["n_residues", js.get("label2", "Published Structure")]),

        # RMSF (full)
        "rmsf_full_mean_1": safe_get(js, ["rmsf", "full", "mean1"]),
        "rmsf_full_mean_2": safe_get(js, ["rmsf", "full", "mean2"]),
        "rmsf_full_ks_p":   safe_get(js, ["rmsf", "full", "ks_p"]),

        # RMSF (aligned)
        "rmsf_aln_mean_1":  safe_get(js, ["rmsf", "aligned", "mean1"]),
        "rmsf_aln_mean_2":  safe_get(js, ["rmsf", "aligned", "mean2"]),
        "rmsf_aln_ks_p":    safe_get(js, ["rmsf", "aligned", "ks_p"]),

        # Radius of gyration
        "rg_mean_1": safe_get(js, ["rg", "mean1"]),
        "rg_mean_2": safe_get(js, ["rg", "mean2"]),
        "rg_ks_p":   safe_get(js, ["rg", "ks_p"]),

        # Backbone RMSD
        "bb_rmsd_mean_1": safe_get(js, ["rmsd_backbone", js.get("label1", "Complete Structure"), "mean"]),
        "bb_rmsd_mean_2": safe_get(js, ["rmsd_backbone", js.get("label2", "Published Structure"), "mean"]),

        # Persistent H‑bonds
        "hbonds_shared":   safe_get(js, ["hbonds", "shared"]),
        "hbonds_excl_1":   safe_get(js, ["hbonds", "exclusive1"]),
        "hbonds_excl_2":   safe_get(js, ["hbonds", "exclusive2"]),

        # Secondary‑structure disagreement
        "dssp_diff_pct":   safe_get(js, ["dssp_diff_pct"]),

        # DCCM bulk metric
        "dccm_mean_delta": safe_get(js, ["dccm_overall", "mean_delta"]),

        # Ligand (if present)
        "lig_present":     "ligand" in js,
        "lig_rmsd_mean_1": safe_get(js, ["ligand", "rmsd", js.get("label1", "Complete Structure"), "mean"]),
        "lig_rmsd_mean_2": safe_get(js, ["ligand", "rmsd", js.get("label2", "Published Structure"), "mean"]),
        "lig_rmsf_mean_1": safe_get(js, ["ligand", "rmsf", js.get("label1", "Complete Structure"), "mean"]),
        "lig_rmsf_mean_2": safe_get(js, ["ligand", "rmsf", js.get("label2", "Published Structure"), "mean"]),
    }

    # Convenience deltas (Complete – Published)
    for m in ["rmsf_full_mean", "rmsf_aln_mean", "rg_mean",
              "bb_rmsd_mean", "lig_rmsd_mean", "lig_rmsf_mean"]:
        if np.isfinite(out.get(f"{m}_1", np.nan)) and np.isfinite(out.get(f"{m}_2", np.nan)):
            out[f"{m}_delta"] = out[f"{m}_1"] - out[f"{m}_2"]
        else:
            out[f"{m}_delta"] = np.nan

    return out


def overall_stats(df: pd.DataFrame, col_pattern: str) -> pd.DataFrame:
    """Return mean / median / std of all columns matching pattern '*_delta'."""
    cols = [c for c in df.columns if c.endswith(col_pattern)]
    return (df[cols]
            .agg(["mean", "median", "std", "count"])
            .round(3)
            .T
            .rename_axis("metric")
            .reset_index())


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────
def main(argv: List[str] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate *_summary.json files from analysisMD and compare Complete vs Published structures.")
    parser.add_argument("json_files", nargs="+",
                        help="One or more summary JSON files (wildcards OK)")
    parser.add_argument("--csv",  metavar="OUT.csv", help="Write per‑system table to CSV")
    parser.add_argument("--plots", metavar="OUT.pdf", help="Save matplotlib box‑plots")
    args = parser.parse_args(argv)

    # Expand wildcards for shells that don't (e.g. Windows)
    files = [str(p) for pattern in args.json_files for p in pathlib.Path(".").glob(pattern)]
    if not files:
        sys.exit("No JSON files matched!")

    # ── Build table ───────────────────────────────────────────────────────────
    records: List[Dict[str, Any]] = []
    for fp in files:
        with open(fp, "r") as fh:
            try:
                js = json.load(fh)
                records.append(extract_metrics(js))
            except Exception as exc:
                print(f"[WARN] Skipping {fp}: {exc}", file=sys.stderr)

    df = pd.DataFrame(records)
    df.sort_values("system", inplace=True)
    print("\nPer‑system metrics (first 10 rows):")
    print(df.head(10).to_markdown(index=False))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\n[OK] Per‑system table written to: {args.csv}")

    # ── Cross‑system statistics ───────────────────────────────────────────────
    delta_cols = [c for c in df.columns if c.endswith("_delta")]
    stats_df = overall_stats(df, "_delta")
    print("\n===  Cross‑system Δ(Complete – Published) statistics  ===")
    print(stats_df.to_markdown(index=False))

    # ── Optional plotting ─────────────────────────────────────────────────────
    if args.plots:
        with plt.rc_context({"figure.figsize": (7, 4)}):
            fig, ax = plt.subplots()
            bp = df.boxplot(column=["rmsf_full_mean_delta", "rg_mean_delta",
                                    "bb_rmsd_mean_delta", "dccm_mean_delta"],
                            ax=ax, return_type="both")
            ax.axhline(0, ls="--", lw=1, zorder=0)
            ax.set_title("Distribution of Complete–Published differences")
            ax.set_ylabel("Δ Metric  (Complete − Published)")
            fig.tight_layout()
            fig.savefig(args.plots)
            plt.close(fig)
            print(f"[OK] Plots written to: {args.plots}")


if __name__ == "__main__":
    main()
