#!/usr/bin/env python3
"""
10_luca_enrichment_fisher_compound_level.py

Compound-level LUCA enrichment using a two-sided Fisher's exact test.

Given two input tables (case vs control) that contain:
  - a compound identifier column
  - a taxon assignment column (default: root_taxon_final)

We define:
  - a compound "binds LUCA" if ANY of its rows is labeled as LUCA
  - non-LUCA rows are (Bacteria/Archaea/Eukaryota variants)
  - rows with other/unknown taxon labels are excluded

Outputs:
  - summary CSV with counts, fractions, odds ratio, two-sided Fisher p-value
  - (optional) per-compound table for both groups

Example:
  python scripts\\09_luca_enrichment_fisher_compound_level.py ^
    --case data\\processed\\final\\pchembl_ge4\\prebiotic_with_root_taxon.csv ^
    --control data\\processed\\final\\clinical_pchembl_ge4\\clinical_with_root_taxon.csv ^
    --case-id-col compound_id ^
    --control-id-col molecule_chembl_id ^
    --outdir results\\pchembl_ge4
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compound-level LUCA enrichment (two-sided Fisher).")
    p.add_argument("--case", type=Path, required=True, help="Case CSV (e.g., prebiotic table)")
    p.add_argument("--control", type=Path, required=True, help="Control CSV (e.g., clinical or matched control)")
    p.add_argument("--case-id-col", type=str, default="compound_id", help="Compound ID column in case CSV")
    p.add_argument("--control-id-col", type=str, default="molecule_chembl_id", help="Compound ID column in control CSV")
    p.add_argument("--taxon-col", type=str, default="root_taxon_final", help="Taxon label column")
    p.add_argument("--outdir", type=Path, default=Path("results/luca_enrichment"), help="Output directory")
    p.add_argument("--tag", type=str, default="", help="Optional tag added to output filenames")
    p.add_argument("--write-compound-table", action="store_true", help="Also write per-compound LUCA flag table")
    return p.parse_args()


def normalize_taxon_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower()


def allowed_taxon_sets() -> Tuple[Set[str], Set[str]]:
    luca_labels = {"luca"}
    non_luca_labels = {
        "bacteria", "bacteriota",
        "archaea", "archaeota",
        "eukaryota", "eukaryote", "eucaryote", "eucaryota",
    }
    return luca_labels, non_luca_labels


def filter_known_taxa(df: pd.DataFrame, taxon_col: str, name: str) -> pd.DataFrame:
    if taxon_col not in df.columns:
        raise ValueError(f"Column '{taxon_col}' not found in {name} dataset.")

    luca_labels, non_luca_labels = allowed_taxon_sets()
    known = luca_labels | non_luca_labels

    tax = normalize_taxon_series(df[taxon_col])
    known_mask = tax.isin(known)

    n_total = len(df)
    n_unknown = int((~known_mask).sum())
    print(f"[{name}] Total rows: {n_total:,}")
    print(f"[{name}] Rows with unknown '{taxon_col}' (excluded): {n_unknown:,}")

    df2 = df.loc[known_mask].copy()
    tax2 = normalize_taxon_series(df2[taxon_col])
    df2["is_luca_row"] = tax2.isin(luca_labels)
    return df2


def compound_level_any_luca(df: pd.DataFrame, id_col: str, group_name: str) -> pd.DataFrame:
    if id_col not in df.columns:
        raise ValueError(f"Expected compound ID column '{id_col}' in {group_name} dataset.")

    out = (
        df.groupby(id_col, as_index=False)["is_luca_row"]
        .any()
        .rename(columns={id_col: "compound_id", "is_luca_row": "binds_luca"})
    )
    out["group"] = group_name
    return out


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """
    Returns: (odds_ratio, p_value_two_sided)
    Contingency:
        LUCA   non-LUCA
    case     a      b
    ctrl     c      d
    """
    # Odds ratio (with standard definition); handle zeros gracefully
    if b == 0 or c == 0:
        # Infinite OR if a>0 and d>0 and (b==0 or c==0). We'll still compute p via fisher.
        or_est = np.inf
    else:
        or_est = (a * d) / (b * c)

    try:
        from scipy.stats import fisher_exact  # type: ignore

        _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        return float(or_est), float(p)
    except Exception:
        # Fallback: compute two-sided Fisher via hypergeometric tail summation
        # This is slower but works without SciPy for moderate counts.
        from math import comb

        def hypergeom_p(x: int, row1: int, col1: int, n: int) -> float:
            return (comb(col1, x) * comb(n - col1, row1 - x)) / comb(n, row1)

        row1 = a + b
        row2 = c + d
        col1 = a + c
        col2 = b + d
        n = row1 + row2

        # Range of possible x values in top-left cell
        xmin = max(0, row1 - col2)
        xmax = min(row1, col1)

        p_obs = hypergeom_p(a, row1, col1, n)

        # Two-sided: sum probabilities <= observed (classic Fisher two-sided definition)
        p_two = 0.0
        for x in range(xmin, xmax + 1):
            px = hypergeom_p(x, row1, col1, n)
            if px <= p_obs + 1e-12:
                p_two += px

        return float(or_est), float(min(1.0, p_two))


def summarize_compounds(df: pd.DataFrame, label: str) -> Dict[str, float]:
    n_total = int(df.shape[0])
    n_luca = int(df["binds_luca"].sum())
    n_non = int(n_total - n_luca)
    frac = float(n_luca / n_total) if n_total else float("nan")
    print(f"\n[{label}] compounds:")
    print(f"  total: {n_total:,}")
    print(f"  LUCA-binding: {n_luca:,}")
    print(f"  non-LUCA: {n_non:,}")
    print(f"  fraction LUCA-binding: {frac:.4f}")
    return {"n_total": n_total, "n_luca": n_luca, "n_non": n_non, "frac_luca": frac}


def main() -> None:
    args = parse_args()

    if not args.case.exists():
        raise FileNotFoundError(f"Case CSV not found: {args.case}")
    if not args.control.exists():
        raise FileNotFoundError(f"Control CSV not found: {args.control}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = f"_{args.tag}" if args.tag else ""

    df_case_raw = pd.read_csv(args.case, dtype=str, low_memory=False)
    df_ctrl_raw = pd.read_csv(args.control, dtype=str, low_memory=False)

    df_case = filter_known_taxa(df_case_raw, args.taxon_col, "case")
    df_ctrl = filter_known_taxa(df_ctrl_raw, args.taxon_col, "control")

    case_comp = compound_level_any_luca(df_case, args.case_id_col, "case")
    ctrl_comp = compound_level_any_luca(df_ctrl, args.control_id_col, "control")

    s_case = summarize_compounds(case_comp, "case")
    s_ctrl = summarize_compounds(ctrl_comp, "control")

    # 2x2 table
    a = s_case["n_luca"]
    b = s_case["n_non"]
    c = s_ctrl["n_luca"]
    d = s_ctrl["n_non"]

    or_est, p_two = fisher_exact_2x2(int(a), int(b), int(c), int(d))

    print("\n2x2 table (compound-level):")
    print("            LUCA   non-LUCA")
    print(f"case      {int(a):6d} {int(b):9d}")
    print(f"control   {int(c):6d} {int(d):9d}")
    print(f"\nOdds ratio: {or_est}")
    print(f"Fisher's exact (two-sided) p-value: {p_two:.6g}")

    # Write summary
    out_summary = args.outdir / f"luca_enrichment_compound_level{tag}.csv"
    summary = pd.DataFrame(
        [
            {
                "case_file": str(args.case),
                "control_file": str(args.control),
                "case_id_col": args.case_id_col,
                "control_id_col": args.control_id_col,
                "taxon_col": args.taxon_col,
                "case_total_compounds": int(s_case["n_total"]),
                "case_luca_compounds": int(s_case["n_luca"]),
                "case_nonluca_compounds": int(s_case["n_non"]),
                "control_total_compounds": int(s_ctrl["n_total"]),
                "control_luca_compounds": int(s_ctrl["n_luca"]),
                "control_nonluca_compounds": int(s_ctrl["n_non"]),
                "odds_ratio": or_est,
                "fisher_p_two_sided": p_two,
            }
        ]
    )
    summary.to_csv(out_summary, index=False)
    print(f"\nWrote summary -> {out_summary}")

    # Optional: write per-compound table
    if args.write_compound_table:
        out_comp = args.outdir / f"compound_level_luca_flags{tag}.csv"
        pd.concat([case_comp, ctrl_comp], ignore_index=True).to_csv(out_comp, index=False)
        print(f"Wrote compound table -> {out_comp}")


if __name__ == "__main__":
    main()
