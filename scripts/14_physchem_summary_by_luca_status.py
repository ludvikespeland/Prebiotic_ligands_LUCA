#!/usr/bin/env python3
"""
14_physchem_summary_by_luca_status.py

Generate SI-ready physicochemical summary statistics comparing
LUCA-binding vs non-LUCA-binding compounds.

Input:
- Compound-level table with:
    - luca_binding_class
    - physicochemical descriptors

Output:
- CSV table with mean, SD, median, IQR, and range per property and group
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PHYS_CHEM_COLS = [
    "MW", "cLogP", "HBD", "HBA",
    "TPSA", "RotB", "AromaticRings", "HeavyAtoms"
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Summarize physicochemical properties by LUCA-binding status"
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Compound-level physchem + LUCA annotation CSV (from script 13)",
    )
    p.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output CSV for SI summary table",
    )
    return p.parse_args()


def median_iqr(x):
    return (
        np.nanmedian(x),
        np.nanpercentile(x, 25),
        np.nanpercentile(x, 75),
    )


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    print(f"Loading physchem table: {args.input}")
    df = pd.read_csv(args.input)

    # Restrict to clearly classified compounds
    luca_df = df[df["luca_binding_class"] == "LUCA_binder"].copy()
    non_luca_df = df[df["luca_binding_class"] == "non_LUCA_binder"].copy()

    print(f"LUCA binders:     {len(luca_df)}")
    print(f"Non-LUCA binders: {len(non_luca_df)}")

    rows = []

    for prop in PHYS_CHEM_COLS:
        for group_name, sub_df in [
            ("LUCA_binder", luca_df),
            ("non_LUCA_binder", non_luca_df),
        ]:
            vals = sub_df[prop]
            med, q1, q3 = median_iqr(vals)

            rows.append({
                "group": group_name,
                "property": prop,
                "N_non_missing": vals.notna().sum(),
                "mean": np.nanmean(vals),
                "std": np.nanstd(vals, ddof=1),
                "min": np.nanmin(vals),
                "q1_25pct": q1,
                "median": med,
                "q3_75pct": q3,
                "max": np.nanmax(vals),
            })

    summary_df = pd.DataFrame(rows)[
        [
            "group", "property", "N_non_missing",
            "mean", "std",
            "min", "q1_25pct", "median", "q3_75pct", "max",
        ]
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.out, index=False)

    print(f"\n SI physchem summary written to:\n{args.out}")


if __name__ == "__main__":
    main()
