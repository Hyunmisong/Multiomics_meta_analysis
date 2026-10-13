"""Project paths and config loading."""
from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"
DATA = ROOT / "data"
GENESETS = DATA / "genesets"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

for _d in (GENESETS, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)


def load_config(name: str) -> dict:
    with open(CONFIG / name, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def resolve_data_path(raw: str) -> Path:
    """Resolve a path from datasets.yaml, relative to the project root."""
    p = Path(raw)
    if p.is_absolute():
        return p
    return (ROOT / p).resolve()
