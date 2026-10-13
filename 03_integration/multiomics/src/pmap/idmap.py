"""UniProt accession -> HGNC symbol mapping, cached on disk.

Only PG2 (proteome) ships accessions with no gene column; everything else
already carries symbols and only needs whitespace/case normalisation.
"""
from __future__ import annotations

import urllib.request
from typing import Dict

import pandas as pd

from .paths import GENESETS

UNIPROT_STREAM = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29"
    "&fields=accession%2Cgene_primary%2Cgene_synonym&format=tsv"
)
CACHE = GENESETS / "uniprot_human_idmap.tsv"


def _download() -> None:
    with urllib.request.urlopen(UNIPROT_STREAM, timeout=300) as resp:
        CACHE.write_bytes(resp.read())


def load_uniprot_map(force: bool = False) -> Dict[str, str]:
    if force or not CACHE.exists():
        _download()
    df = pd.read_csv(CACHE, sep="\t")
    df = df.rename(columns={"Entry": "accession", "Gene Names (primary)": "symbol"})
    df = df.dropna(subset=["symbol"])
    return dict(zip(df["accession"].astype(str), df["symbol"].astype(str).str.upper()))


def clean_symbol(s) -> str:
    """Normalise a gene symbol field; returns '' when unusable."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    s = str(s).strip()
    if not s or s.lower() in {"nan", "none", "-"}:
        return ""
    # multi-gene fields such as "GENEA;GENEB" keep the first entry
    for sep in (";", ",", "|", " "):
        if sep in s:
            s = s.split(sep)[0]
            break
    return s.upper()


def strip_isoform(acc) -> str:
    if acc is None or (isinstance(acc, float) and pd.isna(acc)):
        return ""
    return str(acc).strip().split("-")[0].split(";")[0]
