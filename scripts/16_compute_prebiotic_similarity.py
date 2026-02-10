#!/usr/bin/env python3
"""
15_compute_prebiotic_similarity.py

Compute, for each ChEMBL compound, the MAX Tanimoto similarity to a set of
prebiotic compounds using Morgan fingerprints (ECFP-like).

- Prebiotics: fingerprints generated once (from canonical_smiles, fallback std_smiles)
- ChEMBL: for each compound canonical_smiles -> fingerprint -> max similarity
- Output: input ChEMBL table + max_sim_prebiotic column (or a compact subset if requested)

Default fingerprint: Morgan radius=2, fpSize=2048 (ECFP4-like).

Example:
python scripts\\15_compute_prebiotic_similarity.py ^
  --chembl data\\processed\\final\\physchem\\chembl_compound_level_all_chembl_ge6_with_physchem.csv ^
  --prebiotic data\\processed\\final\\physchem\\near_exact_matches_ge6_prebiotic_compound_level_with_physchem.csv ^
  --outdir data\\processed\\stats\\similarity ^
  --tag ge6
"""

from __future__ import annotations

import argparse
import datetime as _dt
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


def ts() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def first_valid_smiles(row: pd.Series, cols: List[str]) -> Optional[str]:
    for c in cols:
        if c in row and isinstance(row[c], str) and row[c].strip():
            return row[c].strip()
    return None


def smiles_to_fp(smiles: str, gen) -> Optional[object]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return gen.GetFingerprint(mol)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute max Tanimoto similarity to prebiotic set (Morgan/ECFP).")
    p.add_argument("--chembl", type=Path, required=True, help="ChEMBL compound-level CSV with canonical_smiles")
    p.add_argument("--prebiotic", type=Path, required=True, help="Prebiotic compound-level CSV with smiles columns")
    p.add_argument("--outdir", type=Path, default=Path("data/processed/stats/similarity"), help="Output directory")
    p.add_argument("--tag", type=str, default="ge6", help="Tag used in output filename")
    p.add_argument("--chembl-smiles-col", type=str, default="canonical_smiles", help="SMILES column in ChEMBL table")
    p.add_argument(
        "--prebiotic-smiles-cols",
        type=str,
        default="canonical_smiles,std_smiles",
        help="Comma-separated SMILES columns to try for prebiotics (first valid used)",
    )
    p.add_argument("--radius", type=int, default=2, help="Morgan radius (default 2 ~ ECFP4)")
    p.add_argument("--nbits", type=int, default=2048, help="Fingerprint size (default 2048)")
    p.add_argument("--progress-every", type=int, default=50000, help="Print progress every N ChEMBL rows")
    p.add_argument(
        "--compact",
        action="store_true",
        help="Write a compact output table (selected columns) instead of copying all columns",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.chembl.exists():
        raise FileNotFoundError(f"ChEMBL CSV not found: {args.chembl}")
    if not args.prebiotic.exists():
        raise FileNotFoundError(f"Prebiotic CSV not found: {args.prebiotic}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{ts()}] Reading ChEMBL:    {args.chembl}")
    chembl_df = pd.read_csv(args.chembl, low_memory=False)
    print(f"[{ts()}] Reading prebiotic: {args.prebiotic}")
    pre_df = pd.read_csv(args.prebiotic, low_memory=False)

    if args.chembl_smiles_col not in chembl_df.columns:
        raise KeyError(f"ChEMBL missing smiles column '{args.chembl_smiles_col}'. Columns: {list(chembl_df.columns)}")

    pre_smiles_cols = [c.strip() for c in args.prebiotic_smiles_cols.split(",") if c.strip()]
    if not any(c in pre_df.columns for c in pre_smiles_cols):
        raise KeyError(
            f"Prebiotic file has none of the SMILES columns {pre_smiles_cols}. Columns: {list(pre_df.columns)}"
        )

    # Morgan generator
    gen = GetMorganGenerator(radius=args.radius, fpSize=args.nbits)

    # Build prebiotic fingerprints
    print(f"[{ts()}] Building prebiotic fingerprints (radius={args.radius}, nbits={args.nbits})...")
    prebiotic_fps = []
    n_bad = 0
    for _, row in pre_df.iterrows():
        smi = first_valid_smiles(row, pre_smiles_cols)
        if smi is None:
            n_bad += 1
            continue
        fp = smiles_to_fp(smi, gen)
        if fp is None:
            n_bad += 1
            continue
        prebiotic_fps.append(fp)

    prebiotic_fps = list(prebiotic_fps)
    if len(prebiotic_fps) == 0:
        raise RuntimeError("No valid prebiotic fingerprints were generated (check SMILES columns).")

    print(f"[{ts()}] Prebiotic fingerprints: {len(prebiotic_fps):,} (skipped {n_bad:,} rows with missing/invalid SMILES)")

    # Compute max similarity for each ChEMBL compound
    print(f"[{ts()}] Computing max similarity for {len(chembl_df):,} ChEMBL rows...")
    max_sims: List[float] = []
    n_rows = len(chembl_df)

    for i, smi in enumerate(chembl_df[args.chembl_smiles_col].tolist(), start=1):
        fp = smiles_to_fp(smi, gen)
        if fp is None:
            max_sims.append(np.nan)
        else:
            sims = DataStructs.BulkTanimotoSimilarity(fp, prebiotic_fps)
            max_sims.append(float(max(sims)) if sims else np.nan)

        if args.progress_every > 0 and (i % args.progress_every == 0):
            print(f"[{ts()}]  ...processed {i:,}/{n_rows:,}")

    chembl_df["max_sim_prebiotic"] = max_sims

    # Output
    out_file = args.outdir / f"chembl_compound_level_maxsim_prebiotic_{args.tag}.csv"

    if args.compact:
        cols = [
            "compound_id",
            args.chembl_smiles_col,
            "luca_binding_class",
            "is_LUCA_binder",
            "pchembl_best_max",
            "MW", "cLogP", "HBD", "HBA", "TPSA", "RotB", "AromaticRings", "HeavyAtoms",
            "max_sim_prebiotic",
        ]
        cols = [c for c in cols if c in chembl_df.columns]
        chembl_df[cols].to_csv(out_file, index=False)
        print(f"[{ts()}] Wrote (compact): {out_file}")
    else:
        chembl_df.to_csv(out_file, index=False)
        print(f"[{ts()}] Wrote: {out_file}")

    n_nan = int(pd.isna(chembl_df["max_sim_prebiotic"]).sum())
    print(f"[{ts()}] Done. Missing similarities (NaN): {n_nan:,}/{n_rows:,}")


if __name__ == "__main__":
    main()
