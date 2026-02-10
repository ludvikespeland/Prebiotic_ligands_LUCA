#!/usr/bin/env python3
"""
01_chemorigins_standardize_bins_fragments.py

Standardize ChemOrigins molecules, filter to organic/non-metal species,
compute basic properties, bin by molecular weight, and (optionally)
compute Murcko scaffolds and BRICS fragments per bin.

Note:
- Murcko/BRICS fragmentation outputs are included for transparency / exploratory
  characterization and are not required to reproduce the main manuscript analyses.

Environment:
- chem (RDKit; see environment/chem.yml)

Inputs:
- data/raw/chemorigins/molecule-annotations.json  (or user-specified path)

Outputs (default):
- data/processed/chemorigins_bins/
    bin_lt_100_molecules.csv
    bin_100_300_molecules.csv
    bin_300_500_molecules.csv
    bin_500_800_molecules.csv
    bin_gt_800_molecules.csv
    bin_lt_100_unique_fragments.csv
    bin_100_300_unique_fragments.csv
    bin_300_500_unique_fragments.csv
    bin_500_800_unique_fragments.csv
    bin_gt_800_unique_fragments.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize as molst
from rdkit.Chem.Scaffolds import MurckoScaffold

# Silence RDKit logs
RDLogger.DisableLog("rdApp.*")


# ----------------------------
# Defaults (repo-relative)
# ----------------------------
DEFAULT_INPUT = Path("data/raw/chemorigins/molecule-annotations.json")
DEFAULT_OUTDIR = Path("data/processed/chemorigins_bins")


# ----------------------------
# I/O helpers
# ----------------------------
def load_molecule_json(fp: Path) -> pd.DataFrame:
    with fp.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if isinstance(data, list):
        df = pd.DataFrame(data)
    elif isinstance(data, dict) and "molecules" in data:
        df = pd.DataFrame(data["molecules"])
    else:
        raise ValueError(
            "Unexpected JSON structure. Expected a list or a dict with key 'molecules'."
        )

    if "smiles" not in df.columns:
        raise ValueError("No 'smiles' column found in the JSON data.")

    return df


def save_table(df: pd.DataFrame, outdir: Path, name: str) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    fp = outdir / f"{name}.csv"
    df.to_csv(fp, index=False)
    return fp


# ----------------------------
# RDKit standardization
# ----------------------------
_normalizer = molst.Normalizer()
_reionizer = molst.Reionizer()
_uncharger = molst.Uncharger()
_metal_disc = molst.MetalDisconnector()
_lfc = molst.LargestFragmentChooser()

try:
    _taut_enum = molst.TautomerEnumerator()
except Exception:
    _taut_enum = None


def standardize_parent(mol: Chem.Mol, use_tautomer: bool = True) -> Chem.Mol:
    """
    Apply a conservative RDKit standardization pipeline.
    """
    m = _normalizer.normalize(mol)
    m = _reionizer.reionize(m)
    m = _metal_disc.Disconnect(m)
    m = _lfc.choose(m)
    try:
        # newer RDKit builds
        m = molst.ChargeParent(m)
    except AttributeError:
        # fallback
        m = _uncharger.uncharge(m)

    if use_tautomer and _taut_enum is not None:
        m = _taut_enum.Canonicalize(m)

    return m


# ----------------------------
# Filters & properties
# ----------------------------
METAL_ATOMIC_NUMS = set(
    [3, 11, 19, 37, 55, 87, 4, 12, 20, 38, 56, 88]  # alkali/alkaline earth
    + list(range(21, 31))  # transition row 1
    + list(range(39, 49))  # transition row 2
    + list(range(72, 81))  # transition row 3
    + list(range(57, 72))  # lanthanides
    + list(range(89, 104))  # actinides (approx)
    + [13, 31, 49, 50, 81, 82, 83]  # common post-transition metals
)


def is_organic_no_metal(mol: Chem.Mol) -> bool:
    atoms = list(mol.GetAtoms())
    has_carbon = any(a.GetAtomicNum() == 6 for a in atoms)
    has_metal = any(a.GetAtomicNum() in METAL_ATOMIC_NUMS for a in atoms)
    return has_carbon and not has_metal


def compute_props(mol: Chem.Mol) -> Dict[str, Any]:
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Crippen.MolLogP(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "rotB": Lipinski.NumRotatableBonds(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "std_smiles": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
    }


# ----------------------------
# Fragmentation helpers
# ----------------------------
def murcko_core_from_smiles(smi: str) -> Optional[str]:
    try:
        core = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(
            smi, includeChirality=False
        )
        return core if core else None
    except Exception:
        return None


def brics_frags_from_smiles(smi: str) -> List[str]:
    m = Chem.MolFromSmiles(smi)
    if not m:
        return []
    try:
        return sorted(BRICS.BRICSDecompose(m))
    except Exception:
        return []


def list_to_semicolon(lst: List[str]) -> str:
    if not lst:
        return ""
    return ";".join(sorted({x for x in lst if x}))


def annotate_bin(
    df_bin: pd.DataFrame, do_fragments: bool = True
) -> pd.DataFrame:
    """
    Optionally annotate a bin with Murcko scaffold and BRICS fragments.
    """
    df_bin = df_bin.copy()
    if not do_fragments:
        # keep columns present but empty for consistency
        df_bin["murcko"] = None
        df_bin["brics_list"] = [[] for _ in range(len(df_bin))]
        df_bin["brics_unique"] = ""
        return df_bin

    if df_bin.empty:
        df_bin["murcko"] = None
        df_bin["brics_list"] = []
        df_bin["brics_unique"] = ""
        return df_bin

    df_bin["murcko"] = [murcko_core_from_smiles(s) for s in df_bin["std_smiles"]]
    df_bin["brics_list"] = [brics_frags_from_smiles(s) for s in df_bin["std_smiles"]]
    df_bin["brics_unique"] = [list_to_semicolon(frags) for frags in df_bin["brics_list"]]
    return df_bin


def unique_frags_table(df_bin: pd.DataFrame, label: str) -> pd.DataFrame:
    murckos = sorted({m for m in df_bin.get("murcko", []) if m})
    brics = sorted({f for lst in df_bin.get("brics_list", []) for f in (lst or []) if f})
    return pd.DataFrame(
        {
            "bin": [label] * (len(murckos) + len(brics)),
            "type": (["murcko"] * len(murckos)) + (["brics"] * len(brics)),
            "fragment_smiles": murckos + brics,
        }
    )


# ----------------------------
# Binning
# ----------------------------
def split_bins(clean: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    return {
        "lt_100": clean[clean["MW"] < 100].copy(),
        "100_300": clean[(clean["MW"] >= 100) & (clean["MW"] < 300)].copy(),
        "300_500": clean[(clean["MW"] >= 300) & (clean["MW"] < 500)].copy(),
        "500_800": clean[(clean["MW"] >= 500) & (clean["MW"] <= 800)].copy(),
        "gt_800": clean[clean["MW"] > 800].copy(),
    }


# ----------------------------
# Main
# ----------------------------
def main(
    input_path: Path,
    outdir: Optional[Path],
    use_tautomer: bool,
    single_component_only: bool,
    do_fragments: bool,
) -> None:
    df_in = load_molecule_json(input_path)
    print(f"Loaded input records: {len(df_in):,}")
    print("Columns:", ", ".join(df_in.columns[:12]), "..." if len(df_in.columns) > 12 else "")

    rows: List[Dict[str, Any]] = []
    skipped_parse = 0
    skipped_std = 0

    for _, r in df_in.iterrows():
        smi = r.get("smiles")
        if not isinstance(smi, str) or not smi.strip():
            skipped_parse += 1
            continue

        m = Chem.MolFromSmiles(smi)
        if m is None:
            skipped_parse += 1
            continue

        try:
            m_std = standardize_parent(m, use_tautomer=use_tautomer)
        except Exception:
            skipped_std += 1
            continue

        # If requested, drop multi-component post-standardization SMILES
        if single_component_only:
            std_tmp = Chem.MolToSmiles(m_std, canonical=True, isomericSmiles=True)
            if "." in std_tmp:
                continue

        if not is_organic_no_metal(m_std):
            continue

        props = compute_props(m_std)
        rows.append({**r.to_dict(), **props})

    clean = pd.DataFrame(rows)
    if clean.empty:
        raise RuntimeError("No molecules survived standardization + organic/no-metal filters.")

    # Deduplicate (prefer InChIKey if present under common names)
    inchikey_col = None
    for c in ["inchikey", "InChIKey", "standard_inchi_key", "inchi_key"]:
        if c in clean.columns:
            inchikey_col = c
            break

    if inchikey_col:
        clean = clean.drop_duplicates(subset=[inchikey_col])
    else:
        clean = clean.drop_duplicates(subset=["std_smiles"])

    print(
        f"After standardization + organic/no-metal filtering: {len(clean):,} "
        f"(skipped parse: {skipped_parse:,}, skipped std: {skipped_std:,})"
    )

    bins = split_bins(clean)
    print(
        "MW bins:",
        f"<100: {len(bins['lt_100']):,} | "
        f"100–300: {len(bins['100_300']):,} | "
        f"300–500: {len(bins['300_500']):,} | "
        f"500–800: {len(bins['500_800']):,} | "
        f">800: {len(bins['gt_800']):,}"
    )

    # Annotate bins with fragments (optional)
    for k in list(bins.keys()):
        bins[k] = annotate_bin(bins[k], do_fragments=do_fragments)

    # Unique fragments per bin
    uniq = {
        "lt_100": unique_frags_table(bins["lt_100"], "<100"),
        "100_300": unique_frags_table(bins["100_300"], "100-300"),
        "300_500": unique_frags_table(bins["300_500"], "300-500"),
        "500_800": unique_frags_table(bins["500_800"], "500-800"),
        "gt_800": unique_frags_table(bins["gt_800"], ">800"),
    }

    # Preview
    preview = bins["100_300"]
    show_cols = [c for c in ["title", "name", "std_smiles", "MW", "TPSA", "murcko", "brics_unique"] if c in preview.columns]
    if not preview.empty and show_cols:
        print("\nPreview — 100–300 bin (first 5 rows):")
        print(preview[show_cols].head(5).to_string(index=False))

    # Save
    if outdir is not None:
        outdir = Path(outdir)
        save_table(bins["lt_100"], outdir, "bin_lt_100_molecules")
        save_table(bins["100_300"], outdir, "bin_100_300_molecules")
        save_table(bins["300_500"], outdir, "bin_300_500_molecules")
        save_table(bins["500_800"], outdir, "bin_500_800_molecules")
        save_table(bins["gt_800"], outdir, "bin_gt_800_molecules")

        save_table(uniq["lt_100"], outdir, "bin_lt_100_unique_fragments")
        save_table(uniq["100_300"], outdir, "bin_100_300_unique_fragments")
        save_table(uniq["300_500"], outdir, "bin_300_500_unique_fragments")
        save_table(uniq["500_800"], outdir, "bin_500_800_unique_fragments")
        save_table(uniq["gt_800"], outdir, "bin_gt_800_unique_fragments")

        print(f"\nSaved CSVs to: {outdir.resolve()}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Standardize ChemOrigins molecules, bin by MW, and optionally compute Murcko/BRICS fragments."
    )
    p.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Path to ChemOrigins JSON (default: {DEFAULT_INPUT})",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Output directory for CSVs (default: {DEFAULT_OUTDIR}). Use --no-save to skip saving.",
    )
    p.add_argument(
        "--no-save",
        action="store_true",
        help="If set, do not write any CSV files.",
    )
    p.add_argument(
        "--no-tautomer",
        action="store_true",
        help="Disable canonical tautomer enumeration (if supported by your RDKit build).",
    )
    p.add_argument(
        "--single-component-only",
        action="store_true",
        help="Drop molecules that remain multi-component after standardization (SMILES contains '.').",
    )
    p.add_argument(
        "--no-fragments",
        action="store_true",
        help="Skip Murcko/BRICS fragmentation (faster).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        input_path=args.input,
        outdir=None if args.no_save else args.outdir,
        use_tautomer=not args.no_tautomer,
        single_component_only=args.single_component_only,
        do_fragments=not args.no_fragments,
    )
