"""Load every registered dataset, apply all threshold schemes, write tidy tables.

Outputs (results/):
  features_long.parquet   one row per feature per arm per scheme
  gene_calls.parquet      significant features collapsed to genes
  arm_summary.csv         how many genes each arm calls, per scheme
  universes.json          background gene list per arm
"""
import _bootstrap  # noqa: F401

import json

import pandas as pd

from pmap import harmonize
from pmap.paths import RESULTS, load_config


def main():
    cfg = load_config("analysis.yaml")
    print("loading datasets ...")
    long = harmonize.load_all()
    uni = harmonize.universes(long)

    print(f"  {len(long):,} feature rows across {long['arm'].nunique()} arms")
    for arm, n in long.groupby("arm").size().items():
        print(f"    {arm:28s} {n:7,d} features  {long[long.arm==arm].gene.nunique():6,d} genes"
              f"  universe {len(uni[arm]):6,d}")

    scored, gcalls = [], []
    for name, scheme in cfg["threshold_schemes"].items():
        s = harmonize.apply_scheme(long, name, scheme)
        scored.append(s)
        gcalls.append(harmonize.gene_level(s))
        print(f"  scheme {name:12s} -> {int(s['is_sig'].sum()):,} significant features")
    # the genome layer is scored by its own schemes, run alongside
    for name, scheme in (cfg.get("genome_schemes") or {}).items():
        s = harmonize.apply_genome_scheme(long, name, scheme)
        if s.empty:
            continue
        scored.append(s)
        gcalls.append(harmonize.gene_level(s))
        print(f"  genome {name:12s} -> {int(s['is_sig'].sum()):,} genes called "
              f"across {s.loc[s['is_sig'], 'dataset'].nunique()} cohorts")

    scored = pd.concat(scored, ignore_index=True)
    gcalls = pd.concat(gcalls, ignore_index=True)

    scored.to_parquet(RESULTS / "features_long.parquet", index=False)
    gcalls.to_parquet(RESULTS / "gene_calls.parquet", index=False)

    from pmap.aggregate import arm_totals
    arm_totals(gcalls).to_csv(RESULTS / "arm_summary.csv", index=False)
    with open(RESULTS / "universes.json", "w") as fh:
        json.dump({k: sorted(v) for k, v in uni.items()}, fh)

    print(f"written -> {RESULTS}")


if __name__ == "__main__":
    main()
