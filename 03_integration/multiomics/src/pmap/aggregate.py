"""Cross-dataset evidence counting and pathway-level concentration stats.

Nothing here is causal. The question is only: does evidence pile up on the
same pathway across independent cohorts more than you would expect if each
cohort's hits were scattered at random over the genes it actually measured?
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# -------------------------------------------------------------- gene level
def gene_evidence(gene_calls: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """Per gene: how many arms of each layer call it, and which way.

    One row per gene per scheme. Phospho is split by normalisation so a gene
    can be, e.g., phospho-up in 2/3 raw arms but 0/3 protein-adjusted arms --
    which is the whole point of the prototype.
    """
    g = gene_calls[gene_calls["scheme"] == scheme]
    rows = []
    for gene, sub in g.groupby("gene"):
        rec = {"gene": gene, "scheme": scheme}
        for layer, lsub in sub.groupby("layer"):
            if layer == "phospho":
                groups = [(f"phospho_{n}", nsub) for n, nsub in lsub.groupby("norm")]
                groups.append(("phospho_any", lsub))
            else:
                groups = [(layer, lsub)]
            for tag, s in groups:
                rec[f"n_{tag}"] = s["arm"].nunique()
                rec[f"n_{tag}_up"] = s.loc[s["direction"] == "up", "arm"].nunique()
                rec[f"n_{tag}_down"] = s.loc[s["direction"] == "down", "arm"].nunique()
        rows.append(rec)
    ev = pd.DataFrame(rows).fillna(0)
    num = [c for c in ev.columns if c.startswith("n_")]
    ev[num] = ev[num].astype(int)
    if "n_proteome_up" in ev:
        ev["proteome_consistent"] = np.where(
            ev["n_proteome_down"] == 0, "Up",
            np.where(ev["n_proteome_up"] == 0, "Down", "conflict"))
        ev.loc[ev["n_proteome"] == 0, "proteome_consistent"] = ""
    return ev.sort_values("gene").reset_index(drop=True)


def arm_totals(gene_calls: pd.DataFrame) -> pd.DataFrame:
    return (gene_calls.groupby(["scheme", "arm", "layer", "dataset", "norm"])
            .agg(n_genes_called=("gene", "nunique"),
                 n_up=("direction", lambda s: int((s == "up").sum())),
                 n_down=("direction", lambda s: int((s == "down").sum())),
                 n_mixed=("direction", lambda s: int((s == "mixed").sum())))
            .reset_index())


# ----------------------------------------------------------- pathway level
def pathway_evidence(ora: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Per pathway x scheme x direction: in how many arms of each layer is it
    enriched, and how strongly."""
    sig = ora[ora["padj"] <= alpha].copy()
    keys = ["collection", "pathway", "scheme", "direction"]
    rows = []
    for key, sub in sig.groupby(keys):
        rec = dict(zip(keys, key))
        rec["n_arms_total"] = sub["arm"].nunique()
        rec["median_neglog10_padj"] = sub["neglog10_padj"].median()
        rec["max_odds_ratio"] = sub["odds_ratio"].max()
        for layer, lsub in sub.groupby("layer"):
            if layer == "phospho":
                for n, nsub in lsub.groupby("norm"):
                    rec[f"n_phospho_{n}"] = nsub["arm"].nunique()
                rec["n_phospho"] = lsub["arm"].nunique()
            else:
                rec[f"n_{layer}"] = lsub["arm"].nunique()
        rec["arms"] = ";".join(sorted(sub["arm"].unique()))
        rows.append(rec)
    pe = pd.DataFrame(rows).fillna(0)
    num = [c for c in pe.columns if c.startswith("n_")]
    pe[num] = pe[num].astype(int)
    return pe.sort_values(["collection", "n_arms_total", "median_neglog10_padj"],
                          ascending=[True, False, False]).reset_index(drop=True)




# ------------------------------------------------- cross-cohort concentration
def _poisson_binomial_sf(k: int, probs) -> float:
    """P(K >= k) for independent Bernoulli trials with different p."""
    if k <= 0:
        return 1.0
    dist = np.zeros(len(probs) + 1)
    dist[0] = 1.0
    for p in probs:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= (1 - p)
    return float(dist[k:].sum())


def concentration_test(ora: pd.DataFrame,
                       universes: dict,
                       sets: dict,
                       scheme: str,
                       direction: str = "any",
                       layer: str = "proteome",
                       norm: str | None = None,
                       alpha: float = 0.05,
                       min_overlap: int = 3,
                       use_padj: bool = False) -> pd.DataFrame:
    """Does evidence for a pathway concentrate across independent cohorts?

    Statistic: in how many of the N cohorts is the pathway over-represented on
    its own (one-sided Fisher p < alpha against that cohort's own background).

    Null: each cohort d flags a pathway independently with its own observed
    flagging rate pi_d = (#pathways flagged in d) / (#pathways testable in d).
    P(K >= k) is then the Poisson-binomial tail over the cohorts that could
    actually test the pathway. This absorbs the fact that a cohort calling
    4,000 genes flags far more pathways than one calling 300, which is exactly
    the artefact a naive "k of N cohorts" count would read as biology.

    Also reported is a Stouffer combination of the per-cohort Fisher p-values,
    which keeps sub-threshold but consistent evidence instead of discarding it.

    This is an association statistic over cohorts measured on different
    patients. It says evidence piles up on the same pathway; it says nothing
    about direction of causation.
    """
    from scipy.stats import norm as _norm
    from statsmodels.stats.multitest import multipletests

    sub = ora[(ora["scheme"] == scheme) & (ora["direction"] == direction)
              & (ora["layer"] == layer)]
    if norm is not None:
        sub = sub[sub["norm"] == norm]
    if sub.empty:
        return pd.DataFrame()

    pcol = "padj" if use_padj else "pval"
    arms = sorted(sub["arm"].unique())
    arm_ds = {a: sub.loc[sub["arm"] == a, "dataset"].iloc[0] for a in arms}

    # which pathways could each arm test at all (enough members measured)?
    testable = {a: {k for k, g in sets.items()
                    if len(g & universes.get(a, set())) >= min_overlap} for a in arms}
    flagged = {a: set(zip(sub.loc[(sub["arm"] == a) & (sub[pcol] < alpha), "collection"],
                          sub.loc[(sub["arm"] == a) & (sub[pcol] < alpha), "pathway"]))
               for a in arms}
    pi = {a: (len(flagged[a]) / len(testable[a])) if testable[a] else 0.0 for a in arms}
    pmap = {a: dict(zip(zip(sub.loc[sub["arm"] == a, "collection"],
                            sub.loc[sub["arm"] == a, "pathway"]),
                        sub.loc[sub["arm"] == a, pcol])) for a in arms}

    rows = []
    for key in sets:
        tested = [a for a in arms if key in testable[a]]
        if len(tested) < 2:
            continue
        hit_ds = {arm_ds[a] for a in tested if key in flagged[a]}
        k = len(hit_ds)
        # one trial per cohort; when a cohort contributes several arms its
        # mean flagging rate is used so the cohort still counts once
        by_ds = {}
        for a in tested:
            by_ds.setdefault(arm_ds[a], []).append(pi[a])
        probs = [float(np.mean(v)) for v in by_ds.values()]
        pvals = [min((pmap[a].get(key, 1.0) for a in tested if arm_ds[a] == d), default=1.0)
                 for d in by_ds]
        pvals = np.clip(np.array(pvals, dtype=float), 1e-300, 1 - 1e-16)
        z = float(_norm.isf(pvals).sum() / np.sqrt(len(pvals)))
        rows.append({
            "collection": key[0], "pathway": key[1],
            "n_cohorts_testable": len(by_ds), "n_cohorts_enriched": k,
            "expected_cohorts": float(np.sum(probs)),
            "conc_p": _poisson_binomial_sf(k, probs),
            "stouffer_z": z, "stouffer_p": float(_norm.sf(z)),
        })

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    for col, out in [("conc_p", "conc_padj"), ("stouffer_p", "stouffer_padj")]:
        res[out] = np.nan
        for _, idx in res.groupby("collection").groups.items():
            res.loc[idx, out] = multipletests(res.loc[idx, col], method="fdr_bh")[1]
    res["scheme"], res["direction"] = scheme, direction
    res["layer"], res["norm"] = layer, norm or "all"
    return res.sort_values(["collection", "conc_p", "stouffer_p"]).reset_index(drop=True)
