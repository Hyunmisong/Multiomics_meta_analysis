"""Combine the layers on one gene, then ask whether a pathway's members carry
more combined evidence than the background.

This is the other direction from `aggregate.concentration_test`, which keeps
the layers apart and counts cohorts per pathway. Here the evidence is combined
on each gene first, and the pathway is judged from its members.

Why counts are never added across layers: per-cohort call rates differ by two
orders of magnitude (genome 0.1-0.7%, scRNA 2-38%, proteome 14-88%), so a gene
called in 1 of 4 genome cohorts is a far rarer event than one called in 4 of 8
proteome cohorts. Each layer's count is therefore turned into a tail
probability against that layer's own rates before anything is combined.

The genome layer can only ever support "down": a truncating variant is loss of
function. Up-regulation therefore has at most two supporting layers.
"""
from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm
from statsmodels.stats.multitest import multipletests

from .aggregate import _poisson_binomial_sf

LAYERS = ("proteome", "scrna", "genome")
DIRECTIONS = ("up", "down")


def layer_universes(universes: Dict[str, set], layer: str) -> Dict[str, set]:
    """{cohort: measured genes} for one layer."""
    return {k.split("|")[1]: v for k, v in universes.items()
            if k.split("|")[0] == layer}


def gene_evidence(gene_calls: pd.DataFrame, universes: Dict[str, set],
                  layer_schemes: Dict[str, str]) -> pd.DataFrame:
    """Per gene: how many cohorts of each layer call it, and how surprising that is.

    k   cohorts calling the gene in that direction
    m   cohorts of that layer that measured it
    z   Phi^-1(1 - P(K >= k)), the Poisson-binomial tail against each cohort's
        own observed call rate for that direction
    """
    rows = {}
    for layer in LAYERS:
        scheme = layer_schemes.get(layer)
        if scheme is None:
            continue
        uni = layer_universes(universes, layer)
        if not uni:
            continue
        calls = gene_calls[(gene_calls["layer"] == layer)
                           & (gene_calls["scheme"] == scheme)]
        for d in DIRECTIONS:
            if layer == "genome" and d == "up":
                continue            # loss of function has one direction
            called = {c: set(calls.loc[(calls["dataset"] == c)
                                       & (calls["direction"] == d), "gene"])
                      for c in uni}
            rate = {c: len(called[c]) / max(len(uni[c]), 1) for c in uni}
            for g in set().union(*uni.values()):
                here = [c for c in uni if g in uni[c]]
                if not here:
                    continue
                k = sum(1 for c in here if g in called[c])
                r = rows.setdefault((g, d), {})
                r[f"{layer}_m"] = len(here)
                r[f"{layer}_k"] = k
                r[f"{layer}_z"] = (
                    float(norm.isf(_poisson_binomial_sf(k, [rate[c] for c in here])))
                    if k else 0.0)
                if k:
                    r.setdefault("support", []).append(layer)

    out = []
    for (g, d), r in rows.items():
        zs = [r[f"{L}_z"] for L in LAYERS if f"{L}_z" in r]
        rec = {"gene": g, "direction": d, "n_layers_measured": len(zs),
               "n_layers_supporting": len(r.get("support", [])),
               "layers_supporting": ";".join(r.get("support", [])),
               "z_combined": float(np.sum(zs) / np.sqrt(len(zs))) if zs else 0.0}
        for L in LAYERS:
            for s in ("k", "m", "z"):
                rec[f"{L}_{s}"] = r.get(f"{L}_{s}", np.nan)
        out.append(rec)
    ev = pd.DataFrame(out)
    ev["p_combined"] = norm.sf(ev["z_combined"])
    return ev.sort_values(["direction", "z_combined"], ascending=[True, False])


def universes_by_coverage(universes: Dict[str, set], min_layers: int) -> set:
    """Genes measured in at least `min_layers` of the three layers."""
    per = [set().union(*layer_universes(universes, L).values())
           for L in LAYERS if layer_universes(universes, L)]
    from collections import Counter
    c = Counter()
    for s in per:
        c.update(s)
    return {g for g, n in c.items() if n >= min_layers}


# ------------------------------------------------------------ pathway tests
def _bh(df: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
    df = df.copy()
    df["padj"] = np.nan
    for _, idx in df.groupby(keys).groups.items():
        df.loc[idx, "padj"] = multipletests(df.loc[idx, "pval"], method="fdr_bh")[1]
    return df


def pathway_tests(ev: pd.DataFrame, sets: Dict[tuple, set], universe: set,
                  universe_name: str, direction: str,
                  z_cut: float = 1.645, min_members: int = 5,
                  n_perm: int = 2000, seed: int = 0) -> pd.DataFrame:
    """Three competitive tests on the same combined evidence.

    wilcoxon     members' z against every other gene's z. No cutoff, so weak
                 but consistent evidence still counts.
    ora          genes above z_cut are hits; one-sided Fisher. Comparable with
                 the rest of the pipeline, which is cutoff-based.
    permutation  the members' median z against random gene sets of the same
                 size. Corrects for how the universe is composed, but not for
                 genes within a pathway being correlated.
    """
    e = ev[(ev["direction"] == direction) & (ev["gene"].isin(universe))]
    z = dict(zip(e["gene"], e["z_combined"]))
    allz = np.array([z.get(g, 0.0) for g in sorted(universe)])
    hits_all = {g for g in universe if z.get(g, 0.0) >= z_cut}
    rng = np.random.default_rng(seed)

    rows = []
    for (coll, pw), genes in sets.items():
        mem = sorted(genes & universe)
        if len(mem) < min_members:
            continue
        mz = np.array([z.get(g, 0.0) for g in mem])
        others = np.array([z.get(g, 0.0) for g in sorted(universe - set(mem))])
        base = {"collection": coll, "pathway": pw, "direction": direction,
                "universe": universe_name, "n_members": len(mem),
                "n_with_evidence": int((mz > 0).sum()),
                "median_z": float(np.median(mz)),
                "max_z": float(mz.max()) if len(mz) else np.nan}

        u = mannwhitneyu(mz, others, alternative="greater")
        rows.append({**base, "method": "wilcoxon", "stat": float(u.statistic),
                     "pval": float(u.pvalue)})

        a = len(set(mem) & hits_all)
        from scipy.stats import fisher_exact
        orr, p = fisher_exact(
            [[a, len(mem) - a],
             [len(hits_all) - a, len(universe) - len(mem) - len(hits_all) + a]],
            alternative="greater")
        rows.append({**base, "method": "ora", "stat": float(orr), "pval": float(p),
                     "n_hits": a})

        obs = np.median(mz)
        draws = np.median(rng.choice(allz, size=(n_perm, len(mem)), replace=True), axis=1)
        rows.append({**base, "method": "permutation", "stat": obs,
                     "pval": float((np.sum(draws >= obs) + 1) / (n_perm + 1))})
    df = pd.DataFrame(rows)
    return _bh(df, ["collection", "method", "universe", "direction"]) if len(df) else df
