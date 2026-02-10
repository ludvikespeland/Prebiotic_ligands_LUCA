#!/usr/bin/env python3
"""
07_query_chembl_clinical.py

Extract clinical and/or approved small-molecule binding activities
from ChEMBL (SQLite), restricted to SINGLE PROTEIN targets.

Supports pChEMBL thresholds corresponding to 1 µM (>=6)
and 100 µM (>=4).

Outputs:
- Raw activity table
- Aggregated compound × target table
- SMILES file for compounds
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Query ChEMBL for clinical / approved binding compounds.")
    p.add_argument("--db", type=Path, required=True, help="Path to chembl_36.db")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")
    p.add_argument("--pchembl-min", type=float, required=True, help="Minimum pChEMBL value (e.g. 6 or 4)")
    p.add_argument("--approved-only", action="store_true", help="Restrict to approved drugs only")
    p.add_argument("--confidence-min", type=int, default=8)
    p.add_argument("--chunksize", type=int, default=100_000)
    return p.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    if not args.db.exists():
        raise FileNotFoundError(f"ChEMBL DB not found: {args.db}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    phase_clause = "md.max_phase = 4" if args.approved_only else "md.max_phase >= 1"
    tag = "approved" if args.approved_only else "clinical"
    tag += f"_pchembl_ge{str(args.pchembl_min).replace('.', 'p')}"

    out_acts = args.outdir / f"chembl_drugs_activities_{tag}.csv"
    out_roll = args.outdir / f"chembl_drugs_compound_target_{tag}.csv"
    out_smi = args.outdir / f"chembl_drugs_compounds_{tag}.smi"

    print("Using database:", args.db)
    print("Output directory:", args.outdir)
    print("Subset:", tag)

    con = sqlite3.connect(str(args.db))

    sql = f"""
    WITH drugs AS (
      SELECT molregno, chembl_id AS molecule_chembl_id
      FROM molecule_dictionary md
      WHERE {phase_clause}
    ),
    filt_acts AS (
      SELECT a.activity_id, a.assay_id, a.molregno, a.pchembl_value
      FROM activities a
      JOIN assays s ON a.assay_id = s.assay_id
      JOIN target_dictionary t ON s.tid = t.tid
      WHERE s.assay_type = 'B'
        AND s.confidence_score >= {args.confidence_min}
        AND t.target_type = 'SINGLE PROTEIN'
        AND a.pchembl_value >= {args.pchembl_min}
    )
    SELECT
      fa.activity_id,
      d.molecule_chembl_id,
      fa.pchembl_value,
      t.chembl_id AS target_chembl_id,
      t.pref_name AS target_pref_name,
      cs.accession AS uniprot_accession
    FROM filt_acts fa
    JOIN drugs d ON fa.molregno = d.molregno
    JOIN assays s ON fa.assay_id = s.assay_id
    JOIN target_dictionary t ON s.tid = t.tid
    LEFT JOIN target_components tc ON t.tid = tc.tid
    LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id;
    """

    if out_acts.exists():
        out_acts.unlink()

    print("Streaming activities...")
    total = 0
    first = True
    for chunk in pd.read_sql_query(sql, con, chunksize=args.chunksize):
        chunk.to_csv(out_acts, mode="a", header=first, index=False)
        total += len(chunk)
        first = False
        print(f"  wrote {total:,} rows")

    if total == 0:
        con.close()
        raise RuntimeError("Query returned zero rows.")

    print("Aggregating compound × target...")
    acts = pd.read_csv(out_acts)
    agg = (
        acts
        .groupby(["molecule_chembl_id", "target_chembl_id", "target_pref_name", "uniprot_accession"], dropna=False)
        .agg(
            pchembl_median=("pchembl_value", "median"),
            pchembl_best=("pchembl_value", "max"),
            n_measurements=("pchembl_value", "size")
        )
        .reset_index()
    )
    agg.to_csv(out_roll, index=False)
    print("Wrote:", out_roll)

    print("Exporting SMILES...")
    mol_ids = tuple(agg["molecule_chembl_id"].dropna().unique())
    q = f"""
    SELECT md.chembl_id, cs.canonical_smiles
    FROM molecule_dictionary md
    JOIN compound_structures cs ON md.molregno = cs.molregno
    WHERE md.chembl_id IN {mol_ids};
    """
    smi = pd.read_sql_query(q, con)
    smi[["canonical_smiles", "chembl_id"]].dropna().to_csv(
        out_smi, sep="\t", index=False, header=False
    )

    con.close()
    print("Done.")


if __name__ == "__main__":
    main()
