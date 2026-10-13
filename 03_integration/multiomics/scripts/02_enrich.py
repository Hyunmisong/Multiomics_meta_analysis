"""Download gene sets and run over-representation analysis for every
arm x threshold scheme x direction.

Outputs (results/):
  genesets_long.csv.gz   collection / pathway / gene membership actually used
  ora.parquet            full ORA table
  ora_sig.csv            FDR < 0.05 rows only, for eyeballing
"""
import _bootstrap  # noqa: F401

import pandas as pd

from pmap import enrich, genesets
from pmap.paths import RESULTS, load_config


def main():
    cfg = load_config("analysis.yaml")
    print("fetching gene sets ...")
    gs = genesets.load_collections()
    print(gs.groupby("collection")["pathway"].nunique().to_string())
    gs.to_csv(RESULTS / "genesets_long.csv.gz", index=False)
    sets = genesets.geneset_dict(gs)

    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    import json
    with open(RESULTS / "universes.json") as fh:
        uni = {k: set(v) for k, v in json.load(fh).items()}

    print("running ORA ...")
    ora = enrich.run_ora(
        gcalls, uni, sets,
        direction_modes=cfg["direction_modes"],
        min_overlap=cfg["enrichment"]["min_overlap"],
        fdr_method=cfg["enrichment"]["fdr_method"],
    )
    ora.to_parquet(RESULTS / "ora.parquet", index=False)
    sig = ora[ora["padj"] <= 0.05]
    sig.sort_values(["collection", "scheme", "arm", "padj"]).to_csv(
        RESULTS / "ora_sig.csv", index=False)
    print(f"  {len(ora):,} tests, {len(sig):,} enriched at FDR<0.05")
    print(sig.groupby(["scheme", "direction"]).size().to_string())


if __name__ == "__main__":
    main()
