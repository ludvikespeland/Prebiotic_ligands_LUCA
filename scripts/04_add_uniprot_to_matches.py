#!/usr/bin/env python3
"""
04_add_uniprot_to_matches.py

Add UniProt accessions to a match table that
contains ChEMBL target identifiers (target_chembl_id), using an offline ChEMBL
SQLite database.

Outputs:
- exploded: one row per (input row, uniprot accession) where multiple accessions exist
- collapsed: one row per input row with multi-accession fields joined by '|'
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
from typing import List, Optional

import pandas as pd


DEFAULT_INPUT = Path("data/processed/matches/near_exact_matches_inchikey14_pchembl_ge6.csv")
DEFAULT_DB = Path("data/external/chembl_36/chembl_36.db")  # recommended placeholder; override with --db
DEFAULT_OUTDIR = Path("data/processed/matches")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Add UniProt accessions to match tables via ChEMBL SQLite.")
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input CSV containing target_chembl_id")
    p.add_argument("--db", type=Path, required=True, help="Path to ChEMBL SQLite database (chembl_36.db)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory for augmented CSVs")
    p.add_argument("--single-protein-only", action="store_true", help="Keep only targets with target_type='SINGLE PROTEIN'")
    p.add_argument("--human-only", action="store_true", help="Keep only components with organism='Homo sapiens'")
    p.add_argument("--batch-size", type=int, default=900, help="SQLite IN() batch size (keep < 999)")
    return p.parse_args()


def collapse_unique(values: pd.Series) -> pd.NA | str:
    vals = [str(x) for x in pd.Series(values).dropna().unique().tolist()]
    return "|".join(vals) if vals else pd.NA


def fetch_batch(
    conn: sqlite3.Connection,
    ids_batch: List[str],
    single_protein_only: bool,
    human_only: bool,
) -> pd.DataFrame:
    placeholders = ",".join(["?"] * len(ids_batch))
    single_clause = "AND td.target_type = 'SINGLE PROTEIN'" if single_protein_only else ""
    human_clause = "AND cs.organism = 'Homo sapiens'" if human_only else ""

    sql = f"""
        SELECT
          td.chembl_id       AS target_chembl_id,
          td.pref_name       AS target_pref_name,
          td.target_type     AS target_type,
          cs.accession       AS uniprot_accession,
          cs.organism        AS organism,
          cs.description     AS description
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE td.chembl_id IN ({placeholders})
          {single_clause}
          {human_clause}
    """
    return pd.read_sql_query(sql, conn, params=ids_batch)


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input CSV not found: {args.input}")
    if not args.db.exists():
        raise FileNotFoundError(f"ChEMBL SQLite DB not found: {args.db}")

    df_in = pd.read_csv(args.input)
    if "target_chembl_id" not in df_in.columns:
        raise ValueError("Input CSV must contain a column named 'target_chembl_id'.")

    target_ids = df_in["target_chembl_id"].dropna().astype(str).unique().tolist()
    print(f"Loaded {len(df_in):,} rows; mapping {len(target_ids):,} unique target_chembl_id values.")

    conn = sqlite3.connect(str(args.db))

    frames: List[pd.DataFrame] = []
    for i in range(0, len(target_ids), args.batch_size):
        batch = target_ids[i : i + args.batch_size]
        df_map = fetch_batch(conn, batch, args.single_protein_only, args.human_only)
        frames.append(df_map)
        print(f"  batch {i // args.batch_size + 1}: {len(batch)} IDs -> {len(df_map):,} rows")

    conn.close()

    if frames:
        mapping = pd.concat(frames, ignore_index=True)
    else:
        mapping = pd.DataFrame(
            columns=[
                "target_chembl_id",
                "target_pref_name",
                "target_type",
                "uniprot_accession",
                "organism",
                "description",
            ]
        )

    n_input = len(set(target_ids))
    n_with_uniprot = mapping["target_chembl_id"].nunique()
    print("Summary:")
    print(f"  Input unique targets     : {n_input:,}")
    print(f"  Targets with any UniProt : {n_with_uniprot:,}")
    print(f"  Unmapped targets         : {n_input - n_with_uniprot:,}")

    # Exploded merge (one row per (input row, accession))
    df_exploded = df_in.merge(mapping, on="target_chembl_id", how="left")

    # Collapsed mapping (combine multiple values per target with '|')
    collapsed = (
        mapping.groupby("target_chembl_id", as_index=False)
        .agg(
            {
                "target_pref_name": collapse_unique,
                "target_type": collapse_unique,
                "uniprot_accession": collapse_unique,
                "organism": collapse_unique,
                "description": collapse_unique,
            }
        )
    )
    df_collapsed = df_in.merge(collapsed, on="target_chembl_id", how="left")

    # Write outputs
    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = args.input.stem
    out_exploded = args.outdir / f"{stem}_with_uniprot_exploded.csv"
    out_collapsed = args.outdir / f"{stem}_with_uniprot_collapsed.csv"

    df_exploded.to_csv(out_exploded, index=False)
    df_collapsed.to_csv(out_collapsed, index=False)

    print("Wrote:")
    print(f"  Exploded : {len(df_exploded):,} rows -> {out_exploded}")
    print(f"  Collapsed: {len(df_collapsed):,} rows -> {out_collapsed}")


if __name__ == "__main__":
    main()
