#!/usr/bin/env python3
"""
08_match_prebiotic_to_chembl_physchem.py

Physicochemical matching of prebiotic compounds to ChEMBL controls.

- Computes RDKit physicochemical descriptors
- Performs 1:3 nearest-neighbor matching (z-scored, control-based)
- Outputs:
    (1) matched compound table (cases + matched controls)
    (2) matched compound–target activity table

"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


# ---------------------------------------------------------------------
# Descriptor calculation
# ---------------------------------------------------------------------
def calc_descriptors(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcTPSA(mol),
        Descriptors.NumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        mol.GetNumHeavyAtoms(),
    ]


DESCRIPTOR_COLUMNS = [
    "MW", "cLogP", "HBD", "HBA",
    "TPSA", "RotB", "AromaticRings", "HeavyAtoms"
]


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Physicochemical matching of prebiotic compounds to ChEMBL controls"
    )
    p.add_argument(
        "--chembl",
        type=Path,
        required=True,
        help="ChEMBL compound–target CSV containing compound_id and canonical_smiles",
    )
    p.add_argument(
        "--prebiotic",
        type=Path,
        required=True,
        help="CSV listing prebiotic compound_ids (must contain column 'compound_id')",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("data/processed/matched_controls"),
        help="Output directory",
    )
    p.add_argument(
        "--controls-per-case",
        type=int,
        default=3,
        help="Number of matched controls per prebiotic compound (default: 3)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    if not args.chembl.exists():
        raise FileNotFoundError(f"ChEMBL CSV not found: {args.chembl}")
    if not args.prebiotic.exists():
        raise FileNotFoundError(f"Prebiotic CSV not found: {args.prebiotic}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    chembl = pd.read_csv(args.chembl)
    prebiotic = pd.read_csv(args.prebiotic)

    if "compound_id" not in prebiotic.columns:
        raise KeyError("Prebiotic CSV must contain column 'compound_id'")
    if not {"compound_id", "canonical_smiles"}.issubset(chembl.columns):
        raise KeyError("ChEMBL CSV must contain 'compound_id' and 'canonical_smiles'")

    prebiotic_ids = set(prebiotic["compound_id"])

    compounds = (
        chembl[["compound_id", "canonical_smiles"]]
        .drop_duplicates()
        .dropna(subset=["canonical_smiles"])
        .assign(is_prebiotic=lambda x: x["compound_id"].isin(prebiotic_ids))
        .reset_index(drop=True)
    )

    print("Compound counts:")
    print(compounds["is_prebiotic"].value_counts(), "\n")

    print("Calculating RDKit descriptors...")
    desc = compounds["canonical_smiles"].apply(calc_descriptors)
    desc_df = pd.DataFrame(desc.tolist(), columns=DESCRIPTOR_COLUMNS)

    compounds = pd.concat([compounds, desc_df], axis=1).dropna()

    cases = compounds[compounds["is_prebiotic"]].reset_index(drop=True)
    controls = compounds[~compounds["is_prebiotic"]].reset_index(drop=True)

    print(f"Prebiotic compounds: {len(cases)}")
    print(f"Control pool size:   {len(controls)}\n")

    # -----------------------------------------------------------------
    # Z-score features (control-based)
    # -----------------------------------------------------------------
    mu = controls[DESCRIPTOR_COLUMNS].mean()
    sd = controls[DESCRIPTOR_COLUMNS].std()

    cases_X = (cases[DESCRIPTOR_COLUMNS] - mu) / sd
    controls_X = (controls[DESCRIPTOR_COLUMNS] - mu) / sd

    # -----------------------------------------------------------------
    # Nearest-neighbor matching
    # -----------------------------------------------------------------
    print(f"Matching {args.controls_per_case} controls per prebiotic compound...")
    matched_control_ids = set()

    for i in range(len(cases_X)):
        dists = np.linalg.norm(
            controls_X.values - cases_X.iloc[i].values,
            axis=1
        )
        nearest = controls.iloc[
            np.argsort(dists)[: args.controls_per_case]
        ]["compound_id"]
        matched_control_ids.update(nearest)

    matched_controls = controls[
        controls["compound_id"].isin(matched_control_ids)
    ]

    print(f"Matched controls: {len(matched_controls)}\n")

    # -----------------------------------------------------------------
    # Output: matched compounds
    # -----------------------------------------------------------------
    matched_compounds = pd.concat([cases, matched_controls], ignore_index=True)
    matched_compounds["matched_group"] = np.where(
        matched_compounds["is_prebiotic"],
        "prebiotic",
        "matched_control",
    )

    out_compounds = args.outdir / "chembl_prebiotic_physchem_matched_compounds.csv"
    matched_compounds.to_csv(out_compounds, index=False)

    # -----------------------------------------------------------------
    # Output: matched activities
    # -----------------------------------------------------------------
    matched_ids = set(matched_compounds["compound_id"])
    matched_activities = chembl[
        chembl["compound_id"].isin(matched_ids)
    ]

    out_activities = args.outdir / "chembl_prebiotic_physchem_matched_compound_target.csv"
    matched_activities.to_csv(out_activities, index=False)

    # -----------------------------------------------------------------
    # Sanity check
    # -----------------------------------------------------------------
    print("Descriptor means (sanity check):")
    print(
        matched_compounds
        .groupby("matched_group")[DESCRIPTOR_COLUMNS]
        .mean()
        .round(3)
    )

    print("\nDONE")


if __name__ == "__main__":
    main()
