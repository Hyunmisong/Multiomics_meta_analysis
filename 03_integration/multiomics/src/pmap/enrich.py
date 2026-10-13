"""Over-representation analysis, per dataset arm.

Each arm is tested against its own measured background, which matters here:
the phospho arms cover ~3,600 genes while the proteome arms cover ~10,000, and
a shared background would manufacture enrichment for whatever the phospho assay
happens to detect.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def _ora_one(hits: set, universe: set, sets: Dict[tuple, set], min_overlap: int) -> pd.DataFrame:
    rows = []
    n_uni = len(universe)
    n_hit = len(hits)
    if n_hit == 0 or n_uni == 0:
        return pd.DataFrame(columns=[
            "collection", "pathway", "n_pathway_in_universe", "n_hit",
            "n_overlap", "expected", "odds_ratio", "pval", "overlap_genes"])
    for (collection, pathway), genes in sets.items():
        pw = genes & universe
        k = len(pw)
        if k < min_overlap:
            continue
        ov = hits & pw
        a = len(ov)
        if a < min_overlap:
            continue
        table = [[a, n_hit - a], [k - a, n_uni - n_hit - k + a]]
        orr, p = fisher_exact(table, alternative="greater")
        rows.append((collection, pathway, k, n_hit, a, n_hit * k / n_uni, orr, p,
                     ";".join(sorted(ov))))
    return pd.DataFrame(rows, columns=[
        "collection", "pathway", "n_pathway_in_universe", "n_hit",
        "n_overlap", "expected", "odds_ratio", "pval", "overlap_genes"])


def run_ora(gene_calls: pd.DataFrame,
            universes: Dict[str, set],
            sets: Dict[tuple, set],
            direction_modes: Dict[str, list],
            min_overlap: int = 3,
            fdr_method: str = "fdr_bh") -> pd.DataFrame:
    """ORA for every (arm x scheme x direction) cell.

    `direction` is 'any' (DEP vs not), 'up' or 'down'. Genes called 'mixed'
    (phospho sites disagreeing within a gene) join only the 'any' set.
    """
    out = []
    keys = ["arm", "layer", "dataset", "norm", "scheme"]
    for key, grp in gene_calls.groupby(keys):
        arm = key[0]
        uni = universes.get(arm, set())
        for mode, dirs in direction_modes.items():
            for d in dirs:
                hits = set(grp["gene"]) if d == "any" else set(grp.loc[grp["direction"] == d, "gene"])
                hits &= uni
                res = _ora_one(hits, uni, sets, min_overlap)
                if res.empty:
                    continue
                for col, val in zip(keys, key):
                    res[col] = val
                res["direction_mode"] = mode
                res["direction"] = d
                res["n_universe"] = len(uni)
                out.append(res)
    if not out:
        return pd.DataFrame()
    res = pd.concat(out, ignore_index=True)

    # FDR within each independent family of tests
    res["padj"] = np.nan
    fam = ["collection", "dataset", "norm", "scheme", "direction"]
    for _, idx in res.groupby(fam).groups.items():
        res.loc[idx, "padj"] = multipletests(res.loc[idx, "pval"], method=fdr_method)[1]
    res["neglog10_padj"] = -np.log10(res["padj"].clip(lower=1e-300))
    return res
