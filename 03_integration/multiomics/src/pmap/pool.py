"""Pool cohorts into one result, instead of counting how many agree.

The rest of the pipeline asks "in how many cohorts does this reproduce". That
question needs cohorts kept apart. A different question -- "what does this
tumour type look like when all its cohorts are taken together" -- needs them
combined, and that is what this module does.

Two combination rules, because the cohorts do not all carry the same evidence:

  Stouffer   combines per-gene p-values into one signed Z. Uses effect size and
             significance, and needs a p-value, so a cohort that ships only
             direction calls (PE3) cannot take part.
  Vote       counts how many cohorts call the gene in each direction, against a
             null where each cohort flags genes at its own observed rate. Uses
             only the direction, so every cohort takes part.

Running both is deliberate: if they disagree, the disagreement is attributable
to the cohorts that only one of them can use.

Pooling is not a union of gene lists. Cohort depth here spans a factor of nine
(PG1 calls 763 genes, PG4 calls 6,961), so a union would simply be the deepest
cohort with decoration. Combining the statistics instead lets a gene that is
significant in one deep cohort and null in three others be pulled back down.
"""
from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from .aggregate import _poisson_binomial_sf

PCLIP = (1e-300, 1 - 1e-16)


def gene_stats(long: pd.DataFrame, layer: str, cohorts: Iterable[str]) -> pd.DataFrame:
    """One row per gene per cohort: the most significant feature of that gene.

    Phospho-style tables carry several features per gene; the proteome tables
    occasionally do too (isoforms). The strongest evidence is kept, which is
    the same rule `gene_level` uses elsewhere.
    """
    df = long[(long["layer"] == layer) & (long["dataset"].isin(list(cohorts)))].copy()
    df = df.drop_duplicates(["dataset", "feature_id", "gene"])
    df = df.sort_values("pval", na_position="last")
    return df.groupby(["dataset", "gene"], as_index=False).first()


def pooled_universe(universes: Dict[str, set], arms: Iterable[str],
                    min_cohorts: int = 2) -> set:
    """Genes measured in at least `min_cohorts` of the given arms."""
    counts: Dict[str, int] = {}
    for arm in arms:
        for g in universes.get(arm, set()):
            counts[g] = counts.get(g, 0) + 1
    return {g for g, n in counts.items() if n >= min_cohorts}


def stouffer(per_cohort: pd.DataFrame, one_sided: bool = False,
             min_cohorts: int = 2) -> pd.DataFrame:
    """Equal-weight Stouffer combination of per-gene p-values.

    Equal weights, one vote per cohort, deliberately: patient counts span
    25-508 in the genome layer, and weighting by sample size would let one
    large cohort speak for the whole tumour type. A sample-size weighted run is
    the obvious sensitivity check.

    Two-sided input is converted to a signed z so that cohorts disagreeing on
    direction cancel rather than reinforce.
    """
    df = per_cohort.dropna(subset=["pval"]).copy()
    p = df["pval"].clip(*PCLIP)
    if one_sided:
        df["z"] = norm.isf(p)                       # large z = significant
    else:
        sign = np.sign(df["logFC"].fillna(0.0))
        sign[sign == 0] = 1.0
        df["z"] = norm.isf(p / 2) * sign

    g = df.groupby("gene").agg(n_cohorts_used=("z", "size"), z_sum=("z", "sum"))
    g = g[g["n_cohorts_used"] >= min_cohorts]
    g["Z"] = g["z_sum"] / np.sqrt(g["n_cohorts_used"])
    if one_sided:
        g["p_pooled"] = norm.sf(g["Z"])
    else:
        g["p_pooled"] = 2 * norm.sf(g["Z"].abs())
    g["padj_pooled"] = multipletests(g["p_pooled"], method="fdr_bh")[1]
    g["direction"] = np.where(g["Z"] >= 0, "up", "down")
    if one_sided:
        g["direction"] = "down"                     # LOF is loss, by definition
    return g.drop(columns="z_sum").reset_index()


def vote(per_cohort: pd.DataFrame, calls: pd.DataFrame,
         universes: Dict[str, set], layer: str,
         min_cohorts: int = 2) -> pd.DataFrame:
    """Direction-vote test: is a gene called the same way more often than the
    cohorts' own flagging rates predict?

    Null: cohort d calls any given gene "up" with probability
    (genes called up in d) / (genes measured in d), independently. The tail is
    Poisson-binomial over the cohorts that measured the gene, which keeps a
    cohort calling 6,961 genes from dominating one calling 763.
    """
    cohorts = sorted(per_cohort["dataset"].unique())
    rates = {}
    for d in cohorts:
        arm = f"{layer}|{d}|na"
        n_uni = max(len(universes.get(arm, set())), 1)
        sub = calls[calls["dataset"] == d]
        rates[d] = {
            "up": sub.loc[sub["direction"] == "up", "gene"].nunique() / n_uni,
            "down": sub.loc[sub["direction"] == "down", "gene"].nunique() / n_uni,
        }

    measured = {d: universes.get(f"{layer}|{d}|na", set()) for d in cohorts}
    called = {d: dict(zip(calls.loc[calls["dataset"] == d, "gene"],
                          calls.loc[calls["dataset"] == d, "direction"]))
              for d in cohorts}

    genes = sorted(set().union(*measured.values())) if measured else []
    rows = []
    for gene in genes:
        here = [d for d in cohorts if gene in measured[d]]
        if len(here) < min_cohorts:
            continue
        n_up = sum(1 for d in here if called[d].get(gene) == "up")
        n_dn = sum(1 for d in here if called[d].get(gene) == "down")
        if n_up == 0 and n_dn == 0:
            continue
        p_up = _poisson_binomial_sf(n_up, [rates[d]["up"] for d in here]) if n_up else 1.0
        p_dn = _poisson_binomial_sf(n_dn, [rates[d]["down"] for d in here]) if n_dn else 1.0
        rows.append({"gene": gene, "n_measured": len(here), "n_up": n_up,
                     "n_down": n_dn, "p_up": p_up, "p_down": p_dn,
                     "p_vote": min(p_up, p_dn),
                     "direction": "up" if p_up <= p_dn else "down"})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["padj_vote"] = multipletests(out["p_vote"], method="fdr_bh")[1]
    return out


def pooled_gene_set(pooled: pd.DataFrame, universe: set, alpha: float = 0.05,
                    pcol: str = "padj_pooled", direction: str | None = None,
                    top_n: int | None = None) -> set:
    """Significant genes of a pooled result, restricted to the shared universe.

    `top_n` takes the strongest genes by |Z| instead of everything under
    `alpha`, and is the right choice when comparing pooled results to each
    other. Meta-analysing four well-powered proteome cohorts puts about 70% of
    the measured universe under FDR 0.05 -- every small but consistent change
    clears it -- and a gene "set" that is 70% of the background makes an
    over-representation test close to meaningless. Taking a fixed number of
    genes also removes set size as a confounder when two groups are compared,
    which matters because the groups have different depth.
    """
    s = pooled[pooled["gene"].isin(universe)]
    if direction:
        s = s[s["direction"] == direction]
    if top_n is not None:
        s = s.reindex(s["Z"].abs().sort_values(ascending=False).index).head(top_n)
    else:
        s = s[s[pcol] <= alpha]
    return set(s["gene"])
