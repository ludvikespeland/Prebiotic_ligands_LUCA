#!/usr/bin/env python3
"""
11_luca_enrichment_shared_targets.py

Shared-target LUCA enrichment analysis (compound level).

Purpose
-------
Compare LUCA-target engagement between a case set (e.g., prebiotic compounds)
and a control set (e.g., clinical compounds), restricting BOTH datasets to the
intersection of targets hit by both sets ("shared targets only").


Inputs
------
--case:    CSV with at least [case-id-col, target_chembl_id, root_taxon_final]
--control: CSV with at least [control-id-col, target_chembl_id, root_taxon_final]

Outputs
-------
Writes a one-line summary CSV to --outdir and prints a readable console summary.

Example
-------
python scripts\\11_luca_enrichment_shared_targets.py --case data\\processed\\final\\pchembl_ge6\\near_exact_matches_with_root_taxon.csv --control data\\processed\\final\\clinical_pchembl_ge6\\chembl_drugs_compound_target_clinical_pchembl_ge6p0_with_root_taxon.csv --case-id-col compound_id --control-id-col molecule_chembl_id --taxon-col root_taxon_final --target-col target_chembl_id --outdir data\\processed\\stats --tag shared_targets_ge6
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


# ----------------------------
# Taxon parsing
# ----------------------------
def normalize_taxon_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower()


def add_luca_flag_from_root_taxon(df: pd.DataFrame, taxon_col: str, label: str) -> pd.DataFrame:
    if taxon_col not in df.columns:
        raise KeyError(f"[{label}] Missing required column: '{taxon_col}'")

    tax = normalize_taxon_series(df[taxon_col])

    luca_labels = {"luca"}
    non_luca_labels = {
        "bacteria",
        "bacteriota",
        "archaea",
        "archaeota",
        "eukaryota",
        "eukaryote",
        "eucaryote",
        "eucaryota",
    }
    known = luca_labels | non_luca_labels

    keep = tax.isin(known)
    n_total = len(df)
    n_drop = int((~keep).sum())

    df2 = df.loc[keep].copy()
    tax2 = normalize_taxon_series(df2[taxon_col])
    df2["is_luca_row"] = tax2.isin(luca_labels)

    print(f"[{label}] Total rows: {n_total:,}")
    print(f"[{label}] Dropped unknown taxon rows: {n_drop:,}")
    print(f"[{label}] Rows after taxon filtering: {len(df2):,}")
    vc = df2["is_luca_row"].value_counts(dropna=False)
    print(f"[{label}] is_luca_row counts:\n{vc.to_string()}\n")

    return df2


# ----------------------------
# Compound-level aggregation
# ----------------------------
def compound_level_any_luca(
    df: pd.DataFrame, compound_col: str, label: str
) -> pd.DataFrame:
    if compound_col not in df.columns:
        raise KeyError(f"[{label}] Missing required compound id column: '{compound_col}'")
    if "is_luca_row" not in df.columns:
        raise KeyError(f"[{label}] Missing 'is_luca_row' column (internal error).")

    out = (
        df.groupby(compound_col, as_index=False)["is_luca_row"]
        .any()
        .rename(columns={compound_col: "compound_id", "is_luca_row": "binds_luca"})
    )
    return out


def summarize_compounds(df: pd.DataFrame, label: str) -> Tuple[int, int, int, float]:
    n_total = int(len(df))
    n_luca = int(df["binds_luca"].sum())
    n_non = n_total - n_luca
    frac = (n_luca / n_total) if n_total > 0 else float("nan")

    print(f"=== {label} compounds (shared targets only) ===")
    print(f"Total compounds:        {n_total:,}")
    print(f"LUCA-binding compounds: {n_luca:,}")
    print(f"Non-LUCA compounds:     {n_non:,}")
    print(f"Fraction LUCA-binding:  {frac:.4f}\n")

    return n_total, n_luca, n_non, frac


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Shared-target LUCA enrichment (compound-level Fisher's exact test)."
    )
    p.add_argument("--case", type=Path, required=True, help="Case CSV (e.g., prebiotic)")
    p.add_argument("--control", type=Path, required=True, help="Control CSV (e.g., clinical)")
    p.add_argument("--case-id-col", type=str, default="compound_id", help="Case compound id column")
    p.add_argument("--control-id-col", type=str, default="molecule_chembl_id", help="Control compound id column")
    p.add_argument("--taxon-col", type=str, default="root_taxon_final", help="Taxon label column")
    p.add_argument("--target-col", type=str, default="target_chembl_id", help="Target id column")
    p.add_argument("--outdir", type=Path, default=Path("data/processed/stats"), help="Output directory")
    p.add_argument("--tag", type=str, default="", help="Tag appended to output filenames")
    p.add_argument("--write-compound-tables", action="store_true",
                   help="Also write compound-level case/control tables (after shared-target restriction)")
    return p.parse_args()


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    args = parse_args()

    if not args.case.exists():
        raise FileNotFoundError(f"Case file not found: {args.case}")
    if not args.control.exists():
        raise FileNotFoundError(f"Control file not found: {args.control}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    tag = (args.tag.strip().replace(" ", "_"))
    tag_suffix = f"_{tag}" if tag else ""

    print("Loading input files...")
    case_df = pd.read_csv(args.case, low_memory=False)
    ctrl_df = pd.read_csv(args.control, low_memory=False)

    # Add LUCA flags and filter to known taxon categories
    case_df = add_luca_flag_from_root_taxon(case_df, args.taxon_col, "CASE")
    ctrl_df = add_luca_flag_from_root_taxon(ctrl_df, args.taxon_col, "CONTROL")

    # Validate target column
    for df, lbl in [(case_df, "CASE"), (ctrl_df, "CONTROL")]:
        if args.target_col not in df.columns:
            raise KeyError(f"[{lbl}] Missing required column: '{args.target_col}'")

    # Identify shared targets
    case_targets = set(case_df[args.target_col].dropna().astype(str).unique())
    ctrl_targets = set(ctrl_df[args.target_col].dropna().astype(str).unique())
    shared_targets = case_targets & ctrl_targets

    print("Shared target statistics:")
    print(f"Targets hit by CASE:    {len(case_targets):,}")
    print(f"Targets hit by CONTROL: {len(ctrl_targets):,}")
    print(f"Shared targets:         {len(shared_targets):,}\n")

    # Restrict to shared targets only
    case_shared = case_df[case_df[args.target_col].astype(str).isin(shared_targets)].copy()
    ctrl_shared = ctrl_df[ctrl_df[args.target_col].astype(str).isin(shared_targets)].copy()

    print("Rows after restricting to shared targets:")
    print(f"CASE rows:    {len(case_shared):,}")
    print(f"CONTROL rows: {len(ctrl_shared):,}\n")

    # Aggregate to compound level (ANY LUCA)
    case_comp = compound_level_any_luca(case_shared, args.case_id_col, "CASE")
    ctrl_comp = compound_level_any_luca(ctrl_shared, args.control_id_col, "CONTROL")

    # Summaries
    case_n, case_luca, case_non, case_frac = summarize_compounds(case_comp, "CASE")
    ctrl_n, ctrl_luca, ctrl_non, ctrl_frac = summarize_compounds(ctrl_comp, "CONTROL")

    # Fisher's exact test (two-sided)
    table = np.array([[case_luca, case_non], [ctrl_luca, ctrl_non]], dtype=int)
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

    print("2x2 contingency table (shared targets only):")
    print(table)
    print("\nFisher's exact test (two-sided):")
    print(f"Odds ratio: {odds_ratio:.6g}")
    print(f"P-value:    {p_value:.6e}\n")

    # Write summary CSV
    summary = pd.DataFrame(
        [{
            "case_file": str(args.case),
            "control_file": str(args.control),
            "case_id_col": args.case_id_col,
            "control_id_col": args.control_id_col,
            "taxon_col": args.taxon_col,
            "target_col": args.target_col,
            "shared_targets_n": len(shared_targets),
            "case_rows_shared": len(case_shared),
            "control_rows_shared": len(ctrl_shared),
            "case_total_compounds": case_n,
            "case_luca_compounds": case_luca,
            "case_nonluca_compounds": case_non,
            "control_total_compounds": ctrl_n,
            "control_luca_compounds": ctrl_luca,
            "control_nonluca_compounds": ctrl_non,
            "odds_ratio": float(odds_ratio),
            "fisher_p_two_sided": float(p_value),
        }]
    )

    out_summary = args.outdir / f"luca_shared_targets_compound_level{tag_suffix}.csv"
    summary.to_csv(out_summary, index=False)
    print(f"Wrote summary CSV: {out_summary}")

    # Optional: write compound-level tables
    if args.write_compound_tables:
        out_case = args.outdir / f"case_compound_flags_shared_targets{tag_suffix}.csv"
        out_ctrl = args.outdir / f"control_compound_flags_shared_targets{tag_suffix}.csv"
        case_comp.to_csv(out_case, index=False)
        ctrl_comp.to_csv(out_ctrl, index=False)
        print(f"Wrote CASE compound flags: {out_case}")
        print(f"Wrote CONTROL compound flags: {out_ctrl}")


if __name__ == "__main__":
    main()
