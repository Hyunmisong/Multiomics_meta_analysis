"""Cross-dataset evidence counts and pathway-level concentration statistics.

Outputs (results/):
  gene_evidence_<scheme>.csv    per gene: how many arms of each layer, which way
  pathway_evidence.csv          per pathway: in how many arms enriched
  concentration_<scheme>.csv    permutation test on cross-cohort concentration
"""
import _bootstrap  # noqa: F401

import json

import pandas as pd

from pmap import aggregate, genesets
from pmap.paths import RESULTS, load_config


def main():
    cfg = load_config("analysis.yaml")
    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    ora = pd.read_parquet(RESULTS / "ora.parquet")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)
    with open(RESULTS / "universes.json") as fh:
        uni = {k: set(v) for k, v in json.load(fh).items()}

    pe = aggregate.pathway_evidence(ora)
    pe.to_csv(RESULTS / "pathway_evidence.csv", index=False)
    print(f"pathway_evidence: {len(pe):,} rows")

    for scheme in cfg["threshold_schemes"]:
        ev = aggregate.gene_evidence(gcalls, scheme)
        ev.to_csv(RESULTS / f"gene_evidence_{scheme}.csv", index=False)
        print(f"gene_evidence_{scheme}: {len(ev):,} genes")

    # Concentration: proteome cohorts, and each phospho normalisation on its own
    # so the two can be compared on the same footing.
    conc = []
    for scheme in cfg["threshold_schemes"]:
        for direction in ["any", "up", "down"]:
            conc.append(aggregate.concentration_test(
                ora, uni, sets, scheme, direction, layer="proteome",
                alpha=cfg["concentration"]["per_cohort_alpha"],
                use_padj=cfg["concentration"]["per_cohort_use_padj"]))
            for nrm in ["raw", "prot_adj"]:
                conc.append(aggregate.concentration_test(
                    ora, uni, sets, scheme, direction, layer="phospho", norm=nrm,
                    alpha=cfg["concentration"]["per_cohort_alpha"],
                    use_padj=cfg["concentration"]["per_cohort_use_padj"]))
    # genome runs under its own schemes and has only the "down" direction
    for scheme in (cfg.get("genome_schemes") or {}):
        for direction in ["any", "down"]:
            conc.append(aggregate.concentration_test(
                ora, uni, sets, scheme, direction, layer="genome",
                alpha=cfg["concentration"]["per_cohort_alpha"],
                use_padj=cfg["concentration"]["per_cohort_use_padj"]))
    conc = pd.concat([c for c in conc if len(c)], ignore_index=True)
    conc.to_csv(RESULTS / "concentration.csv", index=False)
    hits = conc[conc["conc_padj"] <= 0.05]
    print(f"concentration: {len(conc):,} rows, {len(hits):,} with FDR<0.05")
    print(hits.groupby(["layer", "norm", "scheme", "direction"]).size().to_string())


if __name__ == "__main__":
    main()
