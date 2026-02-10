#!/usr/bin/env python3
"""
06_assign_root_taxon_luca.py

Assign a root taxon (LUCA / Bacteria / Archaea / Eukaryota) to EggNOG orthologous
groups using STRING protein orthology and NCBI taxonomy.

This script implements a conservative, multi-pass strategy:
- Requires cross-domain evidence for LUCA assignment
- Counts DISTINCT member proteins per domain
- Uses EggNOG/STRING orthology + nodes.dmp taxonomy

Outputs the input table augmented with:
- root_taxon_final
- taxid_anchor_final
- rule_applied
- domains_present
- proteins_per_domain
- levels_seen_count
- og_prefix / og_suffix
- is_direct_luca_og
"""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Assign LUCA/root taxon to EggNOG orthologous groups.")
    p.add_argument("--input", type=Path, required=True, help="Input CSV containing EggNOG IDs")
    p.add_argument("--og-col", type=str, default="eggnog_ids", help="Column containing EggNOG OG IDs")
    p.add_argument("--protein-orthology", type=Path, required=True, help="STRING protein.orthology.v12.0.txt")
    p.add_argument("--nodes-dmp", type=Path, required=True, help="NCBI nodes.dmp taxonomy file")
    p.add_argument("--out", type=Path, required=True, help="Output CSV (parent directory will be created)")
    p.add_argument("--chunk-size", type=int, default=2_000_000, help="Chunk size for streaming large files")
    return p.parse_args()


# ---------------------------------------------------------------------
# Taxonomy helpers
# ---------------------------------------------------------------------
CELLULAR = {2, 2157, 2759}
ANCHOR_NAME = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota"}


def load_nodes(nodes_path: Path) -> Dict[int, int]:
    parent_of: Dict[int, int] = {}
    with nodes_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.replace("\t|", "|").split("|")]
            if len(parts) < 2:
                continue
            try:
                taxid = int(parts[0])
                parent = int(parts[1])
            except ValueError:
                continue
            parent_of[taxid] = parent
    return parent_of


def to_anchor(taxid: str, parent_of: Dict[int, int]) -> int | None:
    """Map any taxid to nearest of {2, 2157, 2759}, or None."""
    try:
        x = int(taxid)
    except Exception:
        return None

    seen = set()
    while True:
        if x in CELLULAR:
            return x
        seen.add(x)
        p = parent_of.get(x)
        if p is None or p == x or p in seen:
            return None
        x = p


# ---------------------------------------------------------------------
# OG utilities
# ---------------------------------------------------------------------
PREFIX_RE = re.compile(r"^([A-Za-z]+)")
SUFFIX_RE = re.compile(r"^[A-Za-z]+(.*)$")
DIGITS_RE = re.compile(r"^\d+$")
CANON_PREFIXES = {"kog", "cog", "arcog", "nog"}


def og_prefix(og: str) -> str:
    m = PREFIX_RE.match(og or "")
    return m.group(1) if m else ""


def og_suffix(og: str) -> str:
    m = SUFFIX_RE.match(og or "")
    return m.group(1) if m else ""


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    for p in (args.input, args.protein_orthology, args.nodes_dmp):
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    df = pd.read_csv(args.input, dtype=str, low_memory=False)
    if args.og_col not in df.columns:
        raise KeyError(f"Column '{args.og_col}' not found in input CSV.")

    wanted_ogs: Set[str] = set(df[args.og_col].dropna().astype(str))
    print(f"Input rows: {len(df):,}")
    print(f"Unique OGs requested: {len(wanted_ogs):,}")

    parent_of = load_nodes(args.nodes_dmp)

    # -----------------------------------------------------------------
    # PASS 0: global scan (LUCA OGs + suffix/prefix structure)
    # -----------------------------------------------------------------
    LUCA_OGS = set()
    suffix_to_prefixes = defaultdict(set)
    suffixes_with_direct_luca = set()

    reader = pd.read_csv(
        args.protein_orthology,
        sep="\t",
        dtype=str,
        chunksize=args.chunk_size,
        usecols=["taxonomy_level", "orthologous_group_or_ortholog"],
        low_memory=False,
        on_bad_lines="skip",
    )

    for ch in reader:
        mask = ch["taxonomy_level"].isin(["1", "131567"])
        if mask.any():
            luca_ogs_chunk = set(ch.loc[mask, "orthologous_group_or_ortholog"])
            LUCA_OGS |= luca_ogs_chunk
            for og in luca_ogs_chunk:
                s = og_suffix(og)
                if s:
                    suffixes_with_direct_luca.add(s)

        for og in ch["orthologous_group_or_ortholog"].astype(str):
            p = og_prefix(og)
            s = og_suffix(og)
            if p and s:
                suffix_to_prefixes[s].add(p)

    print(f"PASS0: LUCA OGs = {len(LUCA_OGS):,}")

    # -----------------------------------------------------------------
    # PASS 1: OG → member proteins
    # -----------------------------------------------------------------
    og_to_proteins = defaultdict(set)

    reader = pd.read_csv(
        args.protein_orthology,
        sep="\t",
        dtype=str,
        chunksize=args.chunk_size,
        usecols=["orthologous_group_or_ortholog", "#protein"],
        low_memory=False,
        on_bad_lines="skip",
    )

    for ch in reader:
        sub = ch[ch["orthologous_group_or_ortholog"].isin(wanted_ogs)]
        for og, prot in zip(
            sub["orthologous_group_or_ortholog"].astype(str),
            sub["#protein"].astype(str),
        ):
            og_to_proteins[og].add(prot)

    all_proteins = set().union(*og_to_proteins.values())
    print(f"Distinct member proteins: {len(all_proteins):,}")

    # -----------------------------------------------------------------
    # PASS 2: protein → domain anchors
    # -----------------------------------------------------------------
    prot_to_domains = defaultdict(set)

    reader = pd.read_csv(
        args.protein_orthology,
        sep="\t",
        dtype=str,
        chunksize=args.chunk_size,
        usecols=["#protein", "taxonomy_level"],
        low_memory=False,
        on_bad_lines="skip",
    )

    for ch in reader:
        sub = ch[ch["#protein"].isin(all_proteins)]
        for prot, lv in zip(sub["#protein"].astype(str), sub["taxonomy_level"].astype(str)):
            a = to_anchor(lv, parent_of)
            if a is not None:
                prot_to_domains[prot].add(a)

    # -----------------------------------------------------------------
    # Decision per OG
    # -----------------------------------------------------------------
    rows = []

    for og in wanted_ogs:
        member_prots = og_to_proteins.get(og, set())
        per_domain = {2: set(), 2157: set(), 2759: set()}
        for p in member_prots:
            for a in prot_to_domains.get(p, set()):
                per_domain[a].add(p)

        proteins_per_domain = {a: len(s) for a, s in per_domain.items()}
        domains_present = {a for a, c in proteins_per_domain.items() if c > 0}

        has_cross_domain = len(domains_present) >= 2

        pfx = og_prefix(og)
        sfx = og_suffix(og)
        numeric_suffix = bool(DIGITS_RE.match(sfx or ""))
        interesting_prefixes = {
            p for p in suffix_to_prefixes.get(sfx, set()) if p.lower() in CANON_PREFIXES
        }

        if og in LUCA_OGS and has_cross_domain:
            root = "LUCA"
            taxid = 1
            rule = "direct_LUCA_OG_crossdomain"
            is_direct = True
        elif numeric_suffix and sfx in suffixes_with_direct_luca and len(interesting_prefixes) >= 2 and has_cross_domain:
            root = "LUCA"
            taxid = 1
            rule = "suffix_cross_prefix_crossdomain"
            is_direct = False
        elif has_cross_domain:
            root = "LUCA"
            taxid = 1
            rule = ">=2_superkingdoms_members"
            is_direct = False
        elif 2 in domains_present:
            root, taxid, rule, is_direct = "Bacteria", 2, "single_domain", False
        elif 2157 in domains_present:
            root, taxid, rule, is_direct = "Archaea", 2157, "single_domain", False
        elif 2759 in domains_present:
            root, taxid, rule, is_direct = "Eukaryota", 2759, "single_domain", False
        else:
            root, taxid, rule, is_direct = None, None, "no_evidence", False

        rows.append(
            (
                og,
                root,
                taxid,
                rule,
                ";".join(ANCHOR_NAME[a] for a in sorted(domains_present)),
                proteins_per_domain,
                sum(proteins_per_domain.values()),
                pfx,
                sfx,
                is_direct,
            )
        )

    anc_df = pd.DataFrame(
        rows,
        columns=[
            "og_id",
            "root_taxon_final",
            "taxid_anchor_final",
            "rule_applied",
            "domains_present",
            "proteins_per_domain",
            "levels_seen_count",
            "og_prefix",
            "og_suffix",
            "is_direct_luca_og",
        ],
    )

    out = df.merge(anc_df, left_on=args.og_col, right_on="og_id", how="left").drop(columns=["og_id"])

    # ---- FIX: ensure output directory exists ----
    args.out.parent.mkdir(parents=True, exist_ok=True)

    out.to_csv(args.out, index=False)
    print("Wrote:", args.out)
    print("\nFinal label counts:")
    print(out["root_taxon_final"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
