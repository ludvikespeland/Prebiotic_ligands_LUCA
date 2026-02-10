#!/usr/bin/env python3
"""
17_similarity_summary_luca_vs_nonluca.py

Summarize max_sim_prebiotic distributions for LUCA vs non-LUCA binders.

Input must include:
- luca_binding_class (LUCA_binder / non_LUCA_binder)
- max_sim_prebiotic (float)

Outputs an table with:
N, mean, std, min, q1, median, q3, max,
frac_sim_ge_0.30 / 0.40 / 0.50

Optionally also writes an "enrichment" table:
ratio of LUCA fraction / non-LUCA fraction at each threshold.
"""

from __future__ import annotations

import argparse
import datetime as _dt
from pathlib import Path

import numpy as np
import pandas as pd


def ts() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def summarize_similarity(x: pd.Series, thresholds=(0.30, 0.40, 0.50)) -> dict:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        out = {
            "N": 0, "mean": np.nan, "std": np.nan, "min": np.nan, "q1": np.nan,
            "median": np.nan, "q3": np.nan, "max": np.nan
        }
        for t in thresholds:
            out[f"frac_sim_ge_{t:.2f}"] = np.nan
        return out

    out = {
        "N": int(len(x)),
        "mean": float(x.mean()),
        "std": float(x.std(ddof=1)) if len(x) > 1 else np.nan,
        "min": float(x.min()),
        "q1": float(x.quantile(0.25)),
        "median": float(x.median()),
        "q3": float(x.quantile(0.75)),
        "max": float(x.max()),
    }
    for t in thresholds:
        out[f"frac_sim_ge_{t:.2f}"] = float((x >= t).mean())
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize max similarity distributions for LUCA vs non-LUCA binders.")
    p.add_argument("--input", type=Path, required=True, help="Input similarity table CSV")
    p.add_argument("--outdir", type=Path, default=Path("data/processed/stats/similarity"), help="Output directory")
    p.add_argument("--tag", type=str, default="ge6", help="Tag used in output filenames")
    p.add_argument("--class-col", type=str, default="luca_binding_class", help="LUCA class column")
    p.add_argument("--sim-col", type=str, default="max_sim_prebiotic", help="Similarity column")
    p.add_argument("--write-enrichment", action="store_true", help="Also write LUCA/non-LUCA enrichment ratios table")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input not found: {args.input}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{ts()}] Reading: {args.input}")
    df = pd.read_csv(args.input, low_memory=False)

    if args.class_col not in df.columns:
        raise KeyError(f"Missing '{args.class_col}' in input. Columns: {list(df.columns)}")
    if args.sim_col not in df.columns:
        raise KeyError(f"Missing '{args.sim_col}' in input. Columns: {list(df.columns)}")

    # Keep only LUCA/non-LUCA
    df = df[df[args.class_col].isin(["LUCA_binder", "non_LUCA_binder"])].copy()
    df[args.sim_col] = pd.to_numeric(df[args.sim_col], errors="coerce")
    df = df.dropna(subset=[args.sim_col]).copy()

    print(f"[{ts()}] Rows after filtering LUCA/non-LUCA + non-NaN similarity: {len(df):,}")

    mapping = {
        "LUCA_binder": "LUCA_binder_CHEMBL",
        "non_LUCA_binder": "non_LUCA_binder_CHEMBL",
    }

    rows = []
    for raw, nice in mapping.items():
        sub = df.loc[df[args.class_col] == raw, args.sim_col]
        stats = summarize_similarity(sub, thresholds=(0.30, 0.40, 0.50))
        row = {"group": nice}
        row.update(stats)
        rows.append(row)

    summary_df = pd.DataFrame(rows)

    # Stable column order
    col_order = [
        "group", "N", "mean", "std", "min", "q1", "median", "q3", "max",
        "frac_sim_ge_0.30", "frac_sim_ge_0.40", "frac_sim_ge_0.50",
    ]
    # older names compatibility
    rename_map = {
        "frac_sim_ge_0.30": "frac_sim_ge_0.30",
        "frac_sim_ge_0.40": "frac_sim_ge_0.40",
        "frac_sim_ge_0.50": "frac_sim_ge_0.50",
    }
    summary_df = summary_df.rename(columns=rename_map)

    out_summary = args.outdir / f"SI_similarity_summary_LUCA_vs_nonLUCA_{args.tag}.csv"
    summary_df[col_order].to_csv(out_summary, index=False)
    print(f"[{ts()}] Wrote: {out_summary}")

    if args.write_enrichment:
        # Tail enrichment ratios: (LUCA frac) / (non-LUCA frac)
        luca = summary_df.loc[summary_df["group"] == "LUCA_binder_CHEMBL"].iloc[0]
        non = summary_df.loc[summary_df["group"] == "non_LUCA_binder_CHEMBL"].iloc[0]

        enr_rows = []
        for t in (0.30, 0.40, 0.50):
            key = f"frac_sim_ge_{t:.2f}"
            f_l = float(luca.get(key, np.nan))
            f_n = float(non.get(key, np.nan))
            enr = (f_l / f_n) if (np.isfinite(f_l) and np.isfinite(f_n) and f_n > 0) else np.nan
            enr_rows.append({
                "threshold": t,
                "frac_LUCA": f_l,
                "frac_non_LUCA": f_n,
                "enrichment_ratio_LUCA_over_nonLUCA": enr
            })

        enr_df = pd.DataFrame(enr_rows)
        out_enr = args.outdir / f"SI_similarity_tail_enrichment_{args.tag}.csv"
        enr_df.to_csv(out_enr, index=False)
        print(f"[{ts()}] Wrote: {out_enr}")

    print(f"[{ts()}] Done.")


if __name__ == "__main__":
    main()
