#!/usr/bin/env python3
"""
12_stratified_fisher_by_target_count.py

Stratified (by target count) compound-level LUCA enrichment analysis using Fisher's exact test.

What it does
------------
1) Loads a case dataset (e.g., prebiotic) and a control dataset (e.g., clinical)
2) Filters rows to those with root_taxon_final in {LUCA, Bacteria, Archaea, Eukaryota} (case-insensitive)
3) Collapses to compound-level:
   - binds_luca: True if the compound hits ANY LUCA-classified target
   - n_targets : number of UNIQUE targets hit by the compound (nunique of target_chembl_id)
4) Bins compounds into strata by n_targets (default bins: 1, 2–3, 4+)
5) Runs Fisher's exact test within each stratum
6) Writes a tidy results CSV (and optional compound table)

Inputs must contain
-------------------
- root_taxon_final (or --taxon-col)
- target_chembl_id (or --target-col)
- a compound id column for case and control:
    case:    --case-id-col (default: compound_id)
    control: --control-id-col (default: molecule_chembl_id)

Outputs
-------
- <outdir>/<tag>_stratified_fisher_results.csv
- (optional) <outdir>/<tag>_compound_table.csv

Example
-------
python scripts\\12_stratified_fisher_by_target_count.py --case data\\processed\\final\\pchembl_ge4\\near_exact_matches_with_root_taxon.csv --control data\\processed\\final\\clinical_pchembl_ge4\\chembl_drugs_compound_target_clinical_pchembl_ge4p0_with_root_taxon.csv --case-id-col compound_id --control-id-col molecule_chembl_id --taxon-col root_taxon_final --target-col target_chembl_id --outdir data\\processed\\stats --tag stratified_ge4 --write-compound-table
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


# -----------------------------
# Taxon filtering / LUCA flag
# -----------------------------
def add_luca_flag_from_root_taxon(df: pd.DataFrame, taxon_col: str, name: str) -> pd.DataFrame:
    if taxon_col not in df.columns:
        raise ValueError(f"Column '{taxon_col}' not found in {name} dataset.")

    tax = df[taxon_col].astype(str).str.strip().str.lower()

    luca_labels = {"luca"}
    non_luca_labels = {
        "bacteria", "bacteriota",
        "archaea", "archaeota",
        "eukaryota", "eukaryote", "eucaryote", "eucaryota",
    }

    known = luca_labels | non_luca_labels
    known_mask = tax.isin(known)

    n_total = len(df)
    n_excl = int((~known_mask).sum())
    print(f"[{name}] total rows: {n_total:,}")
    print(f"[{name}] rows excluded (unknown taxon): {n_excl:,}")

    out = df.loc[known_mask].copy()
    tax2 = out[taxon_col].astype(str).str.strip().str.lower()
    out["is_luca_row"] = tax2.isin(luca_labels)

    vc = out["is_luca_row"].value_counts(dropna=False)
    print(f"[{name}] LUCA vs non-LUCA rows after filtering:")
    print(vc.to_string())
    print()
    return out


# -----------------------------
# Compound-level collapse
# -----------------------------
def collapse_to_compounds(
    df: pd.DataFrame,
    *,
    compound_col: str,
    target_col: str,
    name: str,
) -> pd.DataFrame:
    for col in (compound_col, target_col, "is_luca_row"):
        if col not in df.columns:
            raise ValueError(f"[{name}] Missing required column '{col}'.")

    sub = df.dropna(subset=[compound_col, target_col]).copy()

    comp = (
        sub.groupby(compound_col, dropna=False)
        .agg(
            binds_luca=("is_luca_row", "any"),
            n_targets=(target_col, lambda x: int(pd.Series(x).nunique())),
        )
        .reset_index()
        .rename(columns={compound_col: "compound_id"})
    )

    # Remove weird rows (shouldn't happen, but safe)
    comp = comp[comp["n_targets"] > 0].copy()
    return comp


# -----------------------------
# Stratified Fisher
# -----------------------------
def fisher_2x2(pre_luca: int, pre_non: int, ctrl_luca: int, ctrl_non: int):
    table = np.array([[pre_luca, pre_non], [ctrl_luca, ctrl_non]], dtype=int)
    # If any row is empty, fisher_exact is still defined, but you may prefer to skip
    if table.sum() == 0 or table[0].sum() == 0 or table[1].sum() == 0:
        return np.nan, np.nan, table
    orr, p = fisher_exact(table, alternative="two-sided")
    return float(orr), float(p), table


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stratified Fisher's exact test by compound target count.")
    p.add_argument("--case", type=Path, required=True, help="Case CSV (e.g., prebiotic)")
    p.add_argument("--control", type=Path, required=True, help="Control CSV (e.g., clinical)")
    p.add_argument("--case-id-col", type=str, default="compound_id", help="Compound ID column in case CSV")
    p.add_argument("--control-id-col", type=str, default="molecule_chembl_id", help="Compound ID column in control CSV")
    p.add_argument("--taxon-col", type=str, default="root_taxon_final", help="Taxon label column")
    p.add_argument("--target-col", type=str, default="target_chembl_id", help="Target ID column")
    p.add_argument("--outdir", type=Path, default=Path("data/processed/stats"), help="Output directory")
    p.add_argument("--tag", type=str, default="stratified", help="Tag used in output filenames")
    p.add_argument(
        "--bins",
        type=str,
        default="0,1,3,inf",
        help="Comma-separated bin edges for n_targets (default: 0,1,3,inf -> 1 / 2-3 / 4+)",
    )
    p.add_argument(
        "--labels",
        type=str,
        default="1,2-3,4+",
        help="Comma-separated bin labels (must be len(bins)-1). Default: 1,2-3,4+",
    )
    p.add_argument("--write-compound-table", action="store_true", help="Also write the combined compound table CSV")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.case.exists():
        raise FileNotFoundError(f"Case file not found: {args.case}")
    if not args.control.exists():
        raise FileNotFoundError(f"Control file not found: {args.control}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Parse bins/labels
    raw_bins = [x.strip().lower() for x in args.bins.split(",")]
    bins: list[float] = []
    for x in raw_bins:
        if x in {"inf", "+inf", "infty", "infinity"}:
            bins.append(float("inf"))
        else:
            bins.append(float(x))

    labels = [x.strip() for x in args.labels.split(",")]
    if len(labels) != (len(bins) - 1):
        raise ValueError(f"--labels must have len(bins)-1 labels. Got bins={bins}, labels={labels}")

    # Load
    print("Loading input files...")
    case_df = pd.read_csv(args.case, low_memory=False)
    ctrl_df = pd.read_csv(args.control, low_memory=False)
    print(f"Case shape:    {case_df.shape}")
    print(f"Control shape: {ctrl_df.shape}\n")

    # Taxon filter + LUCA flag
    case_df = add_luca_flag_from_root_taxon(case_df, args.taxon_col, "case")
    ctrl_df = add_luca_flag_from_root_taxon(ctrl_df, args.taxon_col, "control")

    # Collapse to compound-level
    case_comp = collapse_to_compounds(
        case_df, compound_col=args.case_id_col, target_col=args.target_col, name="case"
    )
    case_comp["group"] = "case"

    ctrl_comp = collapse_to_compounds(
        ctrl_df, compound_col=args.control_id_col, target_col=args.target_col, name="control"
    )
    ctrl_comp["group"] = "control"

    combined = pd.concat([case_comp, ctrl_comp], ignore_index=True)
    combined["n_targets_bin"] = pd.cut(combined["n_targets"], bins=bins, labels=labels)

    print("Compound-table summary:")
    print(combined.groupby("group")["compound_id"].nunique().rename("unique_compounds").to_string())
    print()

    print("Stratum counts (unique compounds):")
    stratum_counts = (
        combined.groupby(["n_targets_bin", "group"])["compound_id"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(labels)
    )
    print(stratum_counts.to_string())
    print()

    # Run stratified Fishers
    results = []
    for lab in labels:
        sub = combined[combined["n_targets_bin"] == lab].copy()
        if sub.empty:
            results.append(
                {
                    "bin": lab,
                    "case_total": 0,
                    "case_luca": 0,
                    "case_nonluca": 0,
                    "control_total": 0,
                    "control_luca": 0,
                    "control_nonluca": 0,
                    "odds_ratio": np.nan,
                    "fisher_p_two_sided": np.nan,
                    "table": "",
                }
            )
            continue

        case_luca = int(sub[(sub["group"] == "case") & (sub["binds_luca"])].shape[0])
        case_non = int(sub[(sub["group"] == "case") & (~sub["binds_luca"])].shape[0])
        ctrl_luca = int(sub[(sub["group"] == "control") & (sub["binds_luca"])].shape[0])
        ctrl_non = int(sub[(sub["group"] == "control") & (~sub["binds_luca"])].shape[0])

        orr, p, table = fisher_2x2(case_luca, case_non, ctrl_luca, ctrl_non)

        print(f"=== Stratum n_targets_bin = {lab} ===")
        print("2x2 table [[case_LUCA, case_nonLUCA], [control_LUCA, control_nonLUCA]]:")
        print(table)
        if np.isfinite(orr) and np.isfinite(p):
            print(f"Odds ratio: {orr:.6g}")
            print(f"P-value:    {p:.6g}")
        else:
            print("Not enough data to run Fisher (one group missing).")
        print()

        results.append(
            {
                "bin": lab,
                "case_total": int(case_luca + case_non),
                "case_luca": int(case_luca),
                "case_nonluca": int(case_non),
                "control_total": int(ctrl_luca + ctrl_non),
                "control_luca": int(ctrl_luca),
                "control_nonluca": int(ctrl_non),
                "odds_ratio": orr,
                "fisher_p_two_sided": p,
                "table": f"[[{case_luca},{case_non}],[{ctrl_luca},{ctrl_non}]]",
            }
        )

    # Write outputs
    out_res = args.outdir / f"{args.tag}_stratified_fisher_results.csv"
    meta = {
        "case_file": str(args.case),
        "control_file": str(args.control),
        "case_id_col": args.case_id_col,
        "control_id_col": args.control_id_col,
        "taxon_col": args.taxon_col,
        "target_col": args.target_col,
        "bins": args.bins,
        "labels": args.labels,
    }

    res_df = pd.DataFrame(results)
    # Attach metadata columns to every row (easy grepping later)
    for k, v in meta.items():
        res_df.insert(0, k, v)

    res_df.to_csv(out_res, index=False)
    print(f"Wrote stratified results -> {out_res}")

    if args.write_compound_table:
        out_comp = args.outdir / f"{args.tag}_compound_table.csv"
        combined.to_csv(out_comp, index=False)
        print(f"Wrote compound table -> {out_comp}")


if __name__ == "__main__":
    main()
