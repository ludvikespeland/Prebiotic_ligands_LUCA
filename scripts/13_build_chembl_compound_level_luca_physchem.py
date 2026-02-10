#!/usr/bin/env python3
"""
13_build_compound_level_luca_table.py

Build a compound-level table from a LUCA-annotated compound–target table
(prebiotic or ChEMBL), and compute RDKit physicochemical descriptors per compound.

REQUIRED input columns (row-level / compound–target):
- compound_id
- target_chembl_id
- root_taxon_final
- canonical_smiles  (or pass --smiles-col std_smiles)

Compound-level LUCA class (luca_binding_class):
- LUCA_binder      : compound has >=1 row with root_taxon_final == "LUCA" (case-insensitive)
- non_LUCA_binder  : compound has >=1 row with root_taxon_final in {Bacteria, Archaea, Eukaryota} (case-insensitive)
- unknown          : otherwise (e.g., only None/NaN/unclassified)

Outputs (filenames are NOT hardcoded; derived from input and/or --name + --tag):
- <base>_<tag>_compound_level.csv
- <base>_<tag>_compound_level_with_physchem.csv
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


PHYSCHEM_COLS = ["MW", "cLogP", "HBD", "HBA", "TPSA", "RotB", "AromaticRings", "HeavyAtoms"]
NON_LUCA_ROOTS = {"bacteria", "archaea", "eukaryota"}


def ts() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def sanitize_token(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s.strip("_")


def unique_join(values: Iterable, sep: str = ";") -> str:
    vals = []
    for v in values:
        if pd.isna(v):
            continue
        sv = str(v).strip()
        if sv == "" or sv.lower() in {"nan", "none"}:
            continue
        vals.append(sv)
    if not vals:
        return ""
    return sep.join(sorted(set(vals)))


def classify_compound_taxa(taxa_series: pd.Series) -> str:
    taxa = set(str(x).strip().lower() for x in taxa_series.dropna().astype(str))
    if "luca" in taxa:
        return "LUCA_binder"
    if len(taxa.intersection(NON_LUCA_ROOTS)) > 0:
        return "non_LUCA_binder"
    return "unknown"


def calc_physchem(smiles: str) -> dict:
    if not isinstance(smiles, str) or smiles.strip() == "":
        return {k: np.nan for k in PHYSCHEM_COLS}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {k: np.nan for k in PHYSCHEM_COLS}
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Descriptors.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "AromaticRings": Lipinski.NumAromaticRings(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build compound-level LUCA classes + RDKit physchem descriptors.")
    p.add_argument("--input", type=Path, required=True, help="Input LUCA-annotated compound–target CSV")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")
    p.add_argument(
        "--name",
        type=str,
        default="",
        help="Base name for outputs (default: derived from input filename stem)",
    )
    p.add_argument(
        "--tag",
        type=str,
        default="",
        help="Optional tag appended to filenames (e.g., ge6, prebiotic, etc.)",
    )
    p.add_argument("--compound-id-col", type=str, default="compound_id")
    p.add_argument("--target-col", type=str, default="target_chembl_id")
    p.add_argument("--taxon-col", type=str, default="root_taxon_final")
    p.add_argument("--smiles-col", type=str, default="canonical_smiles")
    p.add_argument("--inchi-col", type=str, default="parent_inchi_key")
    p.add_argument("--pchembl-best-col", type=str, default="pchembl_best")
    p.add_argument("--n-measurements-col", type=str, default="n_measurements")
    p.add_argument("--progress-every", type=int, default=20000)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input not found: {args.input}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    base = sanitize_token(args.name) if args.name.strip() else sanitize_token(args.input.stem)
    tag = sanitize_token(args.tag)
    suffix = f"_{tag}" if tag else ""

    out_compound = args.outdir / f"{base}{suffix}_compound_level.csv"
    out_physchem = args.outdir / f"{base}{suffix}_compound_level_with_physchem.csv"

    print(f"[{ts()}] Reading: {args.input}")
    df = pd.read_csv(args.input, low_memory=False)

    # Required columns
    required = [args.compound_id_col, args.target_col, args.taxon_col, args.smiles_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}\nFound columns: {list(df.columns)}")

    # Optional columns present?
    inchi_col = args.inchi_col if args.inchi_col in df.columns else None
    pchembl_best_col = args.pchembl_best_col if args.pchembl_best_col in df.columns else None
    n_meas_col = args.n_measurements_col if args.n_measurements_col in df.columns else None

    print(f"[{ts()}] Input rows: {len(df):,}")
    print(f"[{ts()}] Unique compounds: {df[args.compound_id_col].nunique(dropna=True):,}")

    # Normalize taxon (case-insensitive robust handling)
    tax_norm = df[args.taxon_col].astype(str).str.strip().str.lower()
    df["_taxon_norm"] = tax_norm
    df["_is_luca_row"] = df["_taxon_norm"].eq("luca")

    # Compound-level class
    print(f"[{ts()}] Classifying compounds (LUCA_binder / non_LUCA_binder / unknown)...")
    compound_class = (
        df.groupby(args.compound_id_col, dropna=False)[args.taxon_col]
        .apply(classify_compound_taxa)
        .rename("luca_binding_class")
    )

    # Compound-level aggregation
    print(f"[{ts()}] Aggregating to compound level...")
    gb = df.groupby(args.compound_id_col, dropna=False)

    compounds = pd.DataFrame({
        "canonical_smiles": gb[args.smiles_col].first(),
        "parent_inchi_key": gb[inchi_col].first() if inchi_col else np.nan,
        "pchembl_best_max": gb[pchembl_best_col].max() if pchembl_best_col else np.nan,
        "n_rows": gb.size(),
        "n_measurements_total": gb[n_meas_col].sum() if n_meas_col else 0,
        "n_targets_total": gb[args.target_col].nunique(dropna=True),
        "n_LUCA_rows": gb["_is_luca_row"].sum(),
        "n_LUCA_targets": gb.apply(lambda g: g.loc[g["_is_luca_row"], args.target_col].nunique(dropna=True)),
        "target_chembl_ids": gb[args.target_col].apply(unique_join),
        "root_taxa_seen": gb[args.taxon_col].apply(unique_join),
    })

    compounds = compounds.join(compound_class, how="left")
    compounds = compounds.reset_index().rename(columns={args.compound_id_col: "compound_id"})
    compounds["is_LUCA_binder"] = compounds["luca_binding_class"].eq("LUCA_binder")

    print(f"[{ts()}] Compound-level rows: {len(compounds):,}")
    print(f"[{ts()}] LUCA class counts:\n{compounds['luca_binding_class'].value_counts(dropna=False).to_string()}")

    compounds.to_csv(out_compound, index=False)
    print(f"[{ts()}] Wrote: {out_compound}")

    # Physchem
    print(f"[{ts()}] Computing RDKit physchem for {len(compounds):,} compounds...")
    phys_rows = []
    total = len(compounds)

    for i, smi in enumerate(compounds["canonical_smiles"].tolist(), 1):
        phys_rows.append(calc_physchem(smi))
        if args.progress_every > 0 and (i % args.progress_every == 0):
            print(f"[{ts()}]  ...processed {i:,}/{total:,}")

    phys_df = pd.DataFrame(phys_rows)
    compounds_phys = pd.concat([compounds, phys_df], axis=1)

    compounds_phys.to_csv(out_physchem, index=False)
    print(f"[{ts()}] Wrote: {out_physchem}")

    missing_smiles = compounds["canonical_smiles"].isna().sum()
    all_nan_phys = compounds_phys[PHYSCHEM_COLS].isna().all(axis=1).sum()
    print(f"[{ts()}] Missing SMILES: {missing_smiles:,} | All-NaN physchem rows: {all_nan_phys:,}")
    print(f"[{ts()}] Done.")


if __name__ == "__main__":
    main()
