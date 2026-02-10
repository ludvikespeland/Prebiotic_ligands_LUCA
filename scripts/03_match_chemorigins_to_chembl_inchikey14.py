#!/usr/bin/env python3
"""
03_match_chemorigins_to_chembl_inchikey14.py

Match standardized ChemOrigins molecules to ChEMBL filtered actives using
InChIKey connectivity block (first 14 characters; IK14).

This script is intended to be run separately for different ChEMBL thresholds
(e.g., pChEMBL >= 6 vs pChEMBL >= 4) by pointing --chembl-csv to the relevant
processed file.

Inputs:
- ChemOrigins bins (CSV) containing an InChIKey-like column (e.g. inchi_key, inchikey)
- ChEMBL aggregated compound-target table (CSV) containing parent_inchi_key (preferred)

Output:
- CSV of near-exact matches (IK14) at the compound–target level

Environment:
- stats (pandas)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

import pandas as pd


INCHIKEY_CANDIDATES = ["inchikey", "InChIKey", "INCHIKEY", "inchi_key", "standard_inchi_key", "parent_inchi_key"]


def find_inchikey_col(df: pd.DataFrame, preferred: Optional[str] = None) -> str:
    if preferred and preferred in df.columns:
        return preferred
    for c in INCHIKEY_CANDIDATES:
        if c in df.columns:
            return c
    raise ValueError(f"No InChIKey-like column found. Tried: {INCHIKEY_CANDIDATES}")


def prep_ik14(series: pd.Series) -> pd.DataFrame:
    ik = series.astype(str).str.strip().str.upper()
    ik = ik[ik.str.len() >= 14]
    return pd.DataFrame(
        {
            "IK": ik,
            "IK14": ik.str.slice(0, 14),
        }
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Match ChemOrigins bins to ChEMBL actives by InChIKey IK14.")
    p.add_argument(
        "--prim-bins",
        type=Path,
        nargs="+",
        required=True,
        help="Paths to ChemOrigins bin CSVs (e.g. bin_100_300_molecules.csv ...)",
    )
    p.add_argument(
        "--chembl-csv",
        type=Path,
        required=True,
        help="Path to ChEMBL aggregated compound-target CSV (e.g. data/processed/chembl_36/pchembl_ge6/chembl_actives_compound_target.csv)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("data/processed/matches/near_exact_matches_inchikey14.csv"),
        help="Output CSV path",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Load and concatenate ChemOrigins bins
    prim_parts: List[pd.DataFrame] = []
    for fp in args.prim_bins:
        if not fp.exists():
            raise FileNotFoundError(f"ChemOrigins bin not found: {fp}")
        prim_parts.append(pd.read_csv(fp))
    prim = pd.concat(prim_parts, ignore_index=True)

    pcol = find_inchikey_col(prim)
    prim_ik = prep_ik14(prim[pcol])
    prim = prim.loc[prim_ik.index].copy()
    prim["PRIM_IK"] = prim_ik["IK"].values
    prim["PRIM_IK14"] = prim_ik["IK14"].values

    # Deduplicate by connectivity block
    prim = prim.drop_duplicates(subset=["PRIM_IK14"]).copy()

    # Load ChEMBL aggregated table
    if not args.chembl_csv.exists():
        raise FileNotFoundError(f"ChEMBL CSV not found: {args.chembl_csv}")
    chem = pd.read_csv(args.chembl_csv)

    ccol = find_inchikey_col(chem, preferred="parent_inchi_key")
    chem_ik = prep_ik14(chem[ccol])
    chem = chem.loc[chem_ik.index].copy()
    chem["CHEMBL_IK"] = chem_ik["IK"].values
    chem["CHEMBL_IK14"] = chem_ik["IK14"].values

    # Keep one row per (IK14, target) to avoid redundant duplicates
    if "target_chembl_id" not in chem.columns:
        raise ValueError("Expected column 'target_chembl_id' in ChEMBL CSV.")
    chem = chem.drop_duplicates(subset=["CHEMBL_IK14", "target_chembl_id"]).copy()

    # Join on IK14
    near = prim.merge(
        chem,
        left_on="PRIM_IK14",
        right_on="CHEMBL_IK14",
        how="inner",
        suffixes=("_prim", "_chembl"),
    )

    # Tidy output columns (keep what exists)
    keep = [
        "PRIM_IK",
        "PRIM_IK14",
        "std_smiles",
        "MW",
        "TPSA",
        "cLogP",
        "HBA",
        "HBD",
        "rotB",
        "compound_id",
        "parent_inchi_key",
        "canonical_smiles",
        "target_chembl_id",
        "target_pref_name",
        "pchembl_median",
        "pchembl_best",
        "n_measurements",
    ]
    keep = [c for c in keep if c in near.columns]

    near = near.drop_duplicates(subset=["PRIM_IK14", "target_chembl_id"]).copy()

    # Ensure output directory exists
    args.out.parent.mkdir(parents=True, exist_ok=True)
    near[keep].to_csv(args.out, index=False)

    # Summary
    print(f"Primordial unique IK14: {prim['PRIM_IK14'].nunique():,}")
    print(f"ChEMBL unique IK14:     {chem['CHEMBL_IK14'].nunique():,}")
    print(f"Near-exact matches (IK14): {len(near):,}")
    print("Wrote:", args.out)


if __name__ == "__main__":
    main()
