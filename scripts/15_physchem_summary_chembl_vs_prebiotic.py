#!/usr/bin/env python3
"""
15_physchem_summary_chembl_vs_prebiotic.py

Build an SI-ready physicochemical summary table comparing:
- ChEMBL LUCA binders (luca_binding_class == "LUCA_binder")
- ChEMBL non-LUCA binders (luca_binding_class == "non_LUCA_binder")
- Prebiotic set (all compounds as one group)

Inputs:
1) --chembl : compound-level table with physchem columns AND 'luca_binding_class'
2) --prebiotic : compound-level table with physchem columns (no LUCA split needed)

Outputs:
- <outdir>/physchem_summary_<tag>.csv

Columns:
group, property, N, mean, std, min, q1, median, q3, max
"""

from __future__ import annotations

import argparse
import datetime as _dt
from pathlib import Path

import numpy as np
import pandas as pd


PHYSCHEM_COLS = ["MW", "cLogP", "HBD", "HBA", "TPSA", "RotB", "AromaticRings", "HeavyAtoms"]


def ts() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def summarize_property(x) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "N": 0,
            "mean": np.nan,
            "std": np.nan,
            "min": np.nan,
            "q1": np.nan,
            "median": np.nan,
            "q3": np.nan,
            "max": np.nan,
        }
    return {
        "N": int(x.size),
        "mean": float(np.mean(x)),
        "std": float(np.std(x, ddof=1)) if x.size > 1 else np.nan,
        "min": float(np.min(x)),
        "q1": float(np.percentile(x, 25)),
        "median": float(np.median(x)),
        "q3": float(np.percentile(x, 75)),
        "max": float(np.max(x)),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Physchem SI summary: ChEMBL LUCA/non-LUCA vs Prebiotic (all)"
    )
    p.add_argument("--chembl", type=Path, required=True, help="ChEMBL compound-level physchem CSV")
    p.add_argument("--prebiotic", type=Path, required=True, help="Prebiotic compound-level physchem CSV")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")
    p.add_argument("--tag", type=str, default="ge6", help="Tag for output filename")
    p.add_argument(
        "--prebiotic-group-name",
        type=str,
        default="prebiotic_all",
        help="Group label for prebiotic cohort",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.chembl.exists():
        raise FileNotFoundError(f"ChEMBL input not found: {args.chembl}")
    if not args.prebiotic.exists():
        raise FileNotFoundError(f"Prebiotic input not found: {args.prebiotic}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{ts()}] Reading ChEMBL:   {args.chembl}")
    chembl_df = pd.read_csv(args.chembl, low_memory=False)

    print(f"[{ts()}] Reading Prebiotic:{args.prebiotic}")
    pre_df = pd.read_csv(args.prebiotic, low_memory=False)

    # Normalize RotB name if needed (defensive)
    if "RotB" not in pre_df.columns and "rotB" in pre_df.columns:
        pre_df = pre_df.rename(columns={"rotB": "RotB"})
    if "RotB" not in chembl_df.columns and "rotB" in chembl_df.columns:
        chembl_df = chembl_df.rename(columns={"rotB": "RotB"})

    # Check ChEMBL needs luca_binding_class
    if "luca_binding_class" not in chembl_df.columns:
        raise KeyError(
            "ChEMBL file must contain column 'luca_binding_class' "
            "(expected values include 'LUCA_binder' and 'non_LUCA_binder')."
        )

    # Slice groups
    luca_df = chembl_df[chembl_df["luca_binding_class"] == "LUCA_binder"].copy()
    non_luca_df = chembl_df[chembl_df["luca_binding_class"] == "non_LUCA_binder"].copy()
    pre_all_df = pre_df.copy()

    print(f"[{ts()}] Group sizes:")
    print(f"  LUCA_binder_CHEMBL:     {len(luca_df):,}")
    print(f"  non_LUCA_binder_CHEMBL: {len(non_luca_df):,}")
    print(f"  {args.prebiotic_group_name}:           {len(pre_all_df):,}")

    groups = [
        ("LUCA_binder_CHEMBL", luca_df),
        ("non_LUCA_binder_CHEMBL", non_luca_df),
        (args.prebiotic_group_name, pre_all_df),
    ]

    # Build summary rows
    rows = []
    for prop in PHYSCHEM_COLS:
        for group_name, df_group in groups:
            if prop not in df_group.columns:
                # Skip missing properties (but keep going)
                continue
            stats = summarize_property(df_group[prop])
            rows.append({"group": group_name, "property": prop, **stats})

    summary_df = pd.DataFrame(rows)
    summary_df = summary_df[
        ["group", "property", "N", "mean", "std", "min", "q1", "median", "q3", "max"]
    ]

    out_path = args.outdir / f"physchem_summary_{args.tag}.csv"
    summary_df.to_csv(out_path, index=False)

    print(f"[{ts()}] Wrote: {out_path}")
    print(f"[{ts()}] Preview (first 12 rows):")
    print(summary_df.head(12).to_string(index=False))
    print(f"[{ts()}] Done.")


if __name__ == "__main__":
    main()
