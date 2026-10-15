"""Combine proteome, scRNA and genome on each gene, then judge pathways.

Outputs (results/):
  gene_multiomic.csv       per gene per direction: each layer's cohort count,
                           its surprisal, the combined z, and how many layers
                           support the call
  pathway_multiomic.csv    6 analyses: 2 universes x 3 pathway tests
  pathway_multiomic_summary.csv  how much the 6 agree
"""
import _bootstrap  # noqa: F401

import json

import pandas as pd

from pmap import genesets, multiomic
from pmap.paths import RESULTS, load_config

LAYER_SCHEMES = {"proteome": "fdr05", "scrna": "fdr05", "genome": "g_top200_n2"}


def main():
    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)
    with open(RESULTS / "universes.json") as fh:
        universes = {k: set(v) for k, v in json.load(fh).items()}

    cfg = load_config("analysis.yaml")
    LAYER_SCHEMES["genome"] = (cfg.get("cross_layer_scheme") or {}).get(
        "genome", LAYER_SCHEMES["genome"])

    print("layer schemes:", LAYER_SCHEMES)
    ev = multiomic.gene_evidence(gcalls, universes, LAYER_SCHEMES)
    ev.to_csv(RESULTS / "gene_multiomic.csv", index=False)
    print(f"\ngene evidence: {len(ev):,} gene x direction rows")
    print(ev.groupby(["direction", "n_layers_supporting"]).size()
          .unstack(fill_value=0).to_string())

    u3 = multiomic.universes_by_coverage(universes, 3)
    u2 = multiomic.universes_by_coverage(universes, 2)
    print(f"\nuniverses: all three layers {len(u3):,} | at least two {len(u2):,}")

    out = []
    for uname, uni in [("three_layers", u3), ("two_or_more", u2)]:
        for d in ["up", "down"]:
            r = multiomic.pathway_tests(ev, sets, uni, uname, d)
            if len(r):
                out.append(r)
    res = pd.concat(out, ignore_index=True)
    res.to_csv(RESULTS / "pathway_multiomic.csv", index=False)

    print(f"\n{len(res):,} pathway tests")
    print("\npathways at FDR <= 0.05, by analysis:")
    sig = res[res["padj"] <= 0.05]
    print(sig.pivot_table(index=["universe", "method"], columns="direction",
                          values="pathway", aggfunc="count",
                          fill_value=0).to_string())

    # how much do the six analyses agree?
    rows = []
    for (u, m, d), s in sig.groupby(["universe", "method", "direction"]):
        rows.append({"universe": u, "method": m, "direction": d,
                     "n_sig": len(s),
                     "top": "; ".join(s.nsmallest(5, "padj")["pathway"])})
    summ = pd.DataFrame(rows)
    summ.to_csv(RESULTS / "pathway_multiomic_summary.csv", index=False)

    core = {}
    for (u, m), s in sig[sig.direction == "down"].groupby(["universe", "method"]):
        core[(u, m)] = set(s["collection"] + " | " + s["pathway"])
    if len(core) > 1:
        shared = set.intersection(*core.values())
        print(f"\ndown-regulated pathways found by all {len(core)} analyses: {len(shared)}")
        for k in sorted(shared)[:12]:
            print("   ", k)


if __name__ == "__main__":
    main()
