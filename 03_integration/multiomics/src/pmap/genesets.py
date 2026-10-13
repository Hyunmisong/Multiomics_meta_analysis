"""Download and cache pathway gene sets (KEGG, Reactome, Hallmark).

Gene sets come from the Enrichr GMT endpoint, which serves current,
HGNC-symbol gene sets for all three collections without a login. The raw GMT
is cached verbatim under data/genesets/ so a run is reproducible offline.
"""
from __future__ import annotations

import urllib.request
from typing import Dict

import pandas as pd

from .paths import GENESETS, load_config

ENRICHR_GMT = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}"


def fetch_gmt(library: str, force: bool = False) -> str:
    """Return the GMT text for an Enrichr library, caching it on disk."""
    cache = GENESETS / f"{library}.gmt"
    if cache.exists() and not force:
        return cache.read_text(encoding="utf-8")
    with urllib.request.urlopen(ENRICHR_GMT.format(lib=library), timeout=120) as resp:
        text = resp.read().decode("utf-8")
    cache.write_text(text, encoding="utf-8")
    return text


def parse_gmt(text: str) -> Dict[str, set]:
    sets = {}
    for line in text.strip().split("\n"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        name, genes = parts[0], {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
        if genes:
            sets[name] = genes
    return sets


def load_collections(force: bool = False) -> pd.DataFrame:
    """Return a long table: collection, pathway, gene (one row per membership)."""
    cfg = load_config("analysis.yaml")["genesets"]
    rows = []
    for collection, library in cfg["sources"].items():
        sets = parse_gmt(fetch_gmt(library, force=force))
        for pathway, genes in sets.items():
            if not (cfg["min_size"] <= len(genes) <= cfg["max_size"]):
                continue
            for g in genes:
                rows.append((collection, pathway, g))
    gs = pd.DataFrame(rows, columns=["collection", "pathway", "gene"])
    return gs.drop_duplicates()


def geneset_dict(gs: pd.DataFrame) -> Dict[tuple, set]:
    """{(collection, pathway): {genes}} from the long table."""
    return {k: set(v) for k, v in gs.groupby(["collection", "pathway"])["gene"]}
