#!/usr/bin/env python3
"""
02_query_chembl_binding.py

Offline ChEMBL 36 (SQLite) extraction of high-confidence binding activities
for SINGLE PROTEIN targets, with parent-compound aggregation.

Filters:
- Binding assays (assay_type = 'B')
- Assay confidence >= 8
- Target type = 'SINGLE PROTEIN'
- pChEMBL >= user-specified threshold

Environment:
- stats / chembl (pandas, sqlite3; no RDKit required)

Typical usage:
- pChEMBL >= 6  (≈ 1 µM)
- pChEMBL >= 4  (≈ 100 µM)
"""

from __future__ import annotations

import argparse
import sqlite3
import time
from pathlib import Path
from typing import List

import pandas as pd


# ----------------------------
# Argument parsing
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Query ChEMBL SQLite for filtered binding data.")
    p.add_argument("--db", type=Path, required=True, help="Path to chembl_36.db")
    p.add_argument("--outdir", type=Path, default=Path("data/processed/chembl_36"),
                   help="Base output directory")
    p.add_argument("--pchembl-min", type=float, required=True,
                   help="Minimum pChEMBL value (e.g. 6 or 4)")
    p.add_argument("--confidence-min", type=int, default=8,
                   help="Minimum assay confidence score")
    p.add_argument("--chunksize-sql", type=int, default=250_000)
    p.add_argument("--batch-mol-ids", type=int, default=10_000)
    p.add_argument("--skip-rowcount", action="store_true")
    return p.parse_args()


# ----------------------------
# SQLite helpers
# ----------------------------
def connect_sqlite(db_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(db_path))
    con.execute("PRAGMA journal_mode=OFF;")
    con.execute("PRAGMA synchronous=OFF;")
    con.execute("PRAGMA temp_store=MEMORY;")
    return con


def fetch_mols_batch(con: sqlite3.Connection, ids: List[str]) -> pd.DataFrame:
    qmarks = ",".join(["?"] * len(ids))
    sql = f"""
    SELECT
      m.chembl_id AS molecule_chembl_id,
      cs.canonical_smiles,
      cs.standard_inchi_key AS inchi_key,
      COALESCE(mp.chembl_id, m.chembl_id) AS parent_chembl_id,
      COALESCE(csp.standard_inchi_key, cs.standard_inchi_key) AS parent_inchi_key
    FROM molecule_dictionary m
    LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
    LEFT JOIN molecule_hierarchy  mh ON m.molregno = mh.molregno
    LEFT JOIN molecule_dictionary mp ON mh.parent_molregno = mp.molregno
    LEFT JOIN compound_structures csp ON mh.parent_molregno = csp.molregno
    WHERE m.chembl_id IN ({qmarks})
    """
    return pd.read_sql_query(sql, con, params=ids)


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    args = parse_args()

    if not args.db.exists():
        raise FileNotFoundError(f"ChEMBL database not found: {args.db}")

    tag = f"pchembl_ge{str(args.pchembl_min).replace('.', 'p')}"
    outdir = args.outdir / tag
    outdir.mkdir(parents=True, exist_ok=True)

    raw_csv = outdir / "chembl_actives_activities.csv"
    agg_csv = outdir / "chembl_actives_compound_target.csv"
    smi_path = outdir / "chembl_actives_compounds.smi"

    print("Using database:", args.db)
    print("Output directory:", outdir)
    print(f"Filters: pChEMBL >= {args.pchembl_min}, confidence >= {args.confidence_min}")

    con = connect_sqlite(args.db)

    sql_select = f"""
    SELECT
      a.activity_id,
      m.chembl_id AS molecule_chembl_id,
      a.pchembl_value,
      ass.assay_type,
      ass.confidence_score AS assay_confidence_score,
      t.chembl_id AS target_chembl_id,
      t.pref_name AS target_pref_name,
      t.target_type
    FROM activities a
    JOIN assays ass ON a.assay_id = ass.assay_id
    JOIN target_dictionary t ON ass.tid = t.tid
    JOIN molecule_dictionary m ON a.molregno = m.molregno
    WHERE ass.assay_type = 'B'
      AND ass.confidence_score >= {args.confidence_min}
      AND a.pchembl_value >= {args.pchembl_min}
      AND t.target_type = 'SINGLE PROTEIN'
    """

    if not args.skip_rowcount:
        t0 = time.time()
        sql_count = f"SELECT COUNT(1) FROM ({sql_select})"
        est = pd.read_sql_query(sql_count, con).iloc[0, 0]
        print(f"Estimated matching rows: {est:,} ({time.time()-t0:.1f}s)")

    if raw_csv.exists():
        raw_csv.unlink()

    print("Streaming filtered activities...")
    total = 0
    first = True
    for chunk in pd.read_sql_query(sql_select, con, chunksize=args.chunksize_sql):
        chunk.to_csv(raw_csv, mode="a", header=first, index=False)
        total += len(chunk)
        first = False
        print(f"  wrote {total:,} rows")

    if total == 0:
        con.close()
        raise RuntimeError("Query returned zero rows.")

    print("Collecting unique molecules...")
    mol_ids = set()
    for chunk in pd.read_csv(raw_csv, usecols=["molecule_chembl_id"], chunksize=200_000):
        mol_ids.update(chunk["molecule_chembl_id"].dropna().astype(str))

    mol_ids = sorted(mol_ids)
    print(f"Unique molecules: {len(mol_ids):,}")

    print("Fetching molecule metadata...")
    parts = []
    for i in range(0, len(mol_ids), args.batch_mol_ids):
        parts.append(fetch_mols_batch(con, mol_ids[i:i + args.batch_mol_ids]))
    df_mol = pd.concat(parts, ignore_index=True)

    print("Aggregating per compound–target...")
    df = pd.read_csv(raw_csv)
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df[df["pchembl_value"].notna()]

    df = df.merge(df_mol, on="molecule_chembl_id", how="left")
    df["compound_id"] = df["parent_chembl_id"].fillna(df["molecule_chembl_id"])

    agg = (
        df.groupby(["compound_id", "target_chembl_id"], as_index=False)
          .agg(
              n_measurements=("pchembl_value", "count"),
              pchembl_median=("pchembl_value", "median"),
              pchembl_best=("pchembl_value", "max"),
          )
          .merge(
              df[["compound_id", "canonical_smiles", "parent_inchi_key"]]
              .drop_duplicates("compound_id"),
              on="compound_id",
              how="left",
          )
          .merge(
              df[["target_chembl_id", "target_pref_name"]]
              .drop_duplicates("target_chembl_id"),
              on="target_chembl_id",
              how="left",
          )
    )

    agg.to_csv(agg_csv, index=False)
    print("Wrote aggregated table:", agg_csv)

    uniq = agg[["compound_id", "canonical_smiles", "parent_inchi_key"]].drop_duplicates("compound_id")
    uniq = uniq[uniq["canonical_smiles"].notna()]

    with smi_path.open("w", encoding="utf-8") as f:
        for _, r in uniq.iterrows():
            f.write(f"{r['canonical_smiles']}\t{r['compound_id']}\n")

    print("Wrote SMILES:", smi_path)

    con.close()
    print("Finished.")


if __name__ == "__main__":
    main()
