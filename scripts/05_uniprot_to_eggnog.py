#!/usr/bin/env python3
"""
05_uniprot_to_eggnog.py

Extract EggNOG cross-references for UniProt accessions using
the UniProt idmapping.dat(.gz) file.

Note:
- The idmapping file also contains OMA and OrthoDB identifiers;
  these are parsed but were not used in downstream analyses.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import re
import time
from pathlib import Path
from typing import Dict

import pandas as pd


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def base_uniprot(acc: str) -> str:
    """Normalize UniProt accession (drop isoform suffixes)."""
    if acc is None:
        return ""
    s = str(acc).strip().split()[0]
    return re.sub(r"-\d+$", "", s)


def open_maybe_gzip(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


DB_KEYS: Dict[str, set[str]] = {
    "EggNOG": {"eggnog", "eggnog_id"},
    "OMA": {"oma"},
    "OrthoDB": {"orthodb", "orthodb_id"},
}


def normalize_db_label(db_raw: str) -> str | None:
    key = (db_raw or "").strip().lower()
    for canon, variants in DB_KEYS.items():
        if key == canon.lower() or key in {v.lower() for v in variants}:
            return canon
    return None


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Map UniProt IDs to EggNOG / OMA / OrthoDB using UniProt idmapping.")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--uniprot-col", type=str, default="uniprot_accession")
    p.add_argument("--idmapping", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=Path("data/processed/xrefs"))
    p.add_argument("--merge", action="store_true")
    p.add_argument("--log-every", type=int, default=2_000_000)
    return p.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input CSV not found: {args.input}")
    if not args.idmapping.exists():
        raise FileNotFoundError(f"idmapping file not found: {args.idmapping}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    df_in = pd.read_csv(args.input, dtype=str)
    if args.uniprot_col not in df_in.columns:
        raise KeyError(f"Column '{args.uniprot_col}' not found in input CSV.")

    df_in[args.uniprot_col] = df_in[args.uniprot_col].astype(str).str.strip().str.split().str[0]
    df_in["_uniprot_base"] = df_in[args.uniprot_col].apply(base_uniprot)

    targets = set(df_in["_uniprot_base"].dropna().unique().tolist())
    print(f"[{ts()}] Loaded {len(targets):,} unique UniProt accessions", flush=True)

    collected = {u: {k: set() for k in DB_KEYS} for u in targets}
    long_rows: list[tuple[str, str, str]] = []

    line_count = 0
    kept_count = 0
    start = time.time()

    with open_maybe_gzip(args.idmapping) as fh:
        # 🔴 CRITICAL PROOF-OF-LIFE LINE
        print(f"[{ts()}] idmapping file opened, starting scan", flush=True)

        for line in fh:
            line_count += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            from_id, to_db, to_id = parts[0], parts[1], parts[2]
            if from_id not in targets:
                continue

            canon = normalize_db_label(to_db)
            if not canon:
                continue

            collected[from_id][canon].add(to_id)
            long_rows.append((from_id, canon, to_id))
            kept_count += 1

            if line_count % args.log_every == 0:
                print(f"[{ts()}] Read {line_count:,} lines | kept {kept_count:,}", flush=True)

    elapsed = (time.time() - start) / 60
    print(f"[{ts()}] Finished scan: {line_count:,} lines in {elapsed:.1f} min", flush=True)

    # -----------------------------------------------------------------
    # Write long TSV
    # -----------------------------------------------------------------
    out_long = args.outdir / "uniprot_idmapping_eggnog_oma_orthodb.tsv"
    with out_long.open("w", encoding="utf-8") as fh:
        fh.write("uniprot_base\tdb\txref_id\n")
        for u, db, xid in long_rows:
            fh.write(f"{u}\t{db}\t{xid}\n")

    print(f"[{ts()}] Wrote long mapping TSV -> {out_long}", flush=True)

    # -----------------------------------------------------------------
    # Build wide table
    # -----------------------------------------------------------------
    rows = []
    for u in sorted(targets):
        rows.append({
            "uniprot_base": u,
            "eggnog_ids": ";".join(sorted(collected[u]["EggNOG"])),
            "oma_ids": ";".join(sorted(collected[u]["OMA"])),
            "orthodb_ids": ";".join(sorted(collected[u]["OrthoDB"])),
            "has_eggnog": bool(collected[u]["EggNOG"]),
            "has_oma": bool(collected[u]["OMA"]),
            "has_orthodb": bool(collected[u]["OrthoDB"]),
        })

    wide = pd.DataFrame(rows)
    out_wide = args.outdir / "uniprot_to_eggnog_oma_orthodb.csv"
    wide.to_csv(out_wide, index=False)
    print(f"[{ts()}] Wrote wide mapping CSV -> {out_wide}", flush=True)

    # -----------------------------------------------------------------
    # Optional merge
    # -----------------------------------------------------------------
    if args.merge:
        merged = df_in.merge(wide, left_on="_uniprot_base", right_on="uniprot_base", how="left")
        merged.drop(columns=["_uniprot_base", "uniprot_base"], inplace=True, errors="ignore")
        out_merged = args.outdir / f"{args.input.stem}_with_xrefs.csv"
        merged.to_csv(out_merged, index=False)
        print(f"[{ts()}] Wrote merged CSV -> {out_merged}", flush=True)

    # -----------------------------------------------------------------
    # Coverage summary
    # -----------------------------------------------------------------
    tot = len(targets)
    print(f"[{ts()}] Coverage summary (n={tot} UniProt IDs):", flush=True)
    for k in ["EggNOG", "OMA", "OrthoDB"]:
        c = wide[f"has_{k.lower()}"].sum()
        print(f"  {k:8s}: {c}/{tot} ({100*c/tot:.1f}%)", flush=True)

    any_hit = (wide["has_eggnog"] | wide["has_oma"] | wide["has_orthodb"]).sum()
    print(f"  Any of the three: {any_hit}/{tot} ({100*any_hit/tot:.1f}%)", flush=True)


if __name__ == "__main__":
    main()
