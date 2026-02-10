#!/usr/bin/env python3
"""
09_split_physchem_matched_sets.py

Split an annotated compound–target table (already run through 04–06)
into prebiotic vs matched-control subsets, using the compound membership
table produced by 08_match_prebiotic_to_chembl_physchem.py.

Inputs
- --annotated: compound–target CSV that already contains LUCA/root_taxon annotations
               (e.g., chembl_prebiotic_physchem_matched_compound_target_with_root_taxon.csv)
- --matched-compounds: output from script 08 containing compound_id + is_prebiotic
                       (e.g., chembl_prebiotic_physchem_matched_compounds.csv)

Outputs
- <outdir>/<tag>_prebiotic.csv
- <outdir>/<tag>_matched_controls.csv
- <outdir>/<tag>_split_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Split an annotated compound–target table into prebiotic vs matched controls using script-08 membership."
    )
    p.add_argument("--annotated", type=Path, required=True, help="Annotated compound–target CSV (post 04–06)")
    p.add_argument(
        "--matched-compounds",
        type=Path,
        required=True,
        help="Matched compounds CSV from script 08 (must contain compound_id and is_prebiotic)",
    )
    p.add_argument("--compound-id-col", type=str, default="compound_id", help="Compound ID column name in both files")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")
    p.add_argument(
        "--tag",
        type=str,
        default="physchem_matched",
        help="Filename tag/prefix for outputs (default: physchem_matched)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.annotated.exists():
        raise FileNotFoundError(f"Annotated CSV not found: {args.annotated}")
    if not args.matched_compounds.exists():
        raise FileNotFoundError(f"Matched compounds CSV not found: {args.matched_compounds}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    annotated = pd.read_csv(args.annotated, low_memory=False)
    members = pd.read_csv(args.matched_compounds, low_memory=False)

    cid = args.compound_id_col
    if cid not in annotated.columns:
        raise KeyError(f"Annotated table missing '{cid}'. Columns: {list(annotated.columns)}")
    if cid not in members.columns:
        raise KeyError(f"Matched-compounds table missing '{cid}'. Columns: {list(members.columns)}")

    if "is_prebiotic" not in members.columns:
        raise KeyError(
            "Matched-compounds table must contain column 'is_prebiotic' (True/False). "
            "This comes from script 08 output."
        )

    # Normalize IDs to string for safe joins
    annotated[cid] = annotated[cid].astype(str)
    members[cid] = members[cid].astype(str)

    # Build membership map
    member_map = (
        members[[cid, "is_prebiotic"]]
        .drop_duplicates(subset=[cid])
        .copy()
    )

    # Merge membership onto annotated rows
    merged = annotated.merge(member_map, on=cid, how="left")

    # Sanity: how many annotated rows have membership?
    n_total_rows = len(merged)
    n_with_flag = merged["is_prebiotic"].notna().sum()
    n_missing_flag = n_total_rows - n_with_flag

    # Keep only rows that are part of the matched set
    in_set = merged[merged["is_prebiotic"].notna()].copy()

    prebiotic_rows = in_set[in_set["is_prebiotic"] == True].copy()
    control_rows = in_set[in_set["is_prebiotic"] == False].copy()

    # Summary counts (row-level and compound-level)
    def n_unique(df: pd.DataFrame) -> int:
        return df[cid].nunique(dropna=True)

    summary = pd.DataFrame(
        [
            {
                "annotated_file": str(args.annotated),
                "matched_compounds_file": str(args.matched_compounds),
                "rows_total_annotated": n_total_rows,
                "rows_with_membership_flag": int(n_with_flag),
                "rows_missing_membership_flag": int(n_missing_flag),
                "rows_in_matched_set": len(in_set),
                "rows_prebiotic": len(prebiotic_rows),
                "rows_controls": len(control_rows),
                "unique_compounds_in_set": n_unique(in_set),
                "unique_compounds_prebiotic": n_unique(prebiotic_rows),
                "unique_compounds_controls": n_unique(control_rows),
            }
        ]
    )

    out_pre = args.outdir / f"{args.tag}_prebiotic.csv"
    out_ctl = args.outdir / f"{args.tag}_matched_controls.csv"
    out_sum = args.outdir / f"{args.tag}_split_summary.csv"

    prebiotic_rows.to_csv(out_pre, index=False)
    control_rows.to_csv(out_ctl, index=False)
    summary.to_csv(out_sum, index=False)

    print("Wrote:")
    print("  Prebiotic:", out_pre)
    print("  Controls :", out_ctl)
    print("  Summary  :", out_sum)

    print("\nKey counts:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
