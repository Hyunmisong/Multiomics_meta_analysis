"""Does the gastric / oesophageal split matter, or can the cohorts be pooled?

This is a variance-partitioning question, not a differential test: are cohorts
of the same tumour type more like each other than like cohorts of the other
type? Two ways of asking it are implemented.

  Cohort clustering   each cohort is a point in pathway-enrichment space;
                      ANOSIM and PERMANOVA ask whether the tumour-type labels
                      explain the layout. Answers "can the proteome cohorts be
                      pooled".
  Layer concordance   pool each tumour type's cohorts into one result per
                      layer, then ask whether the genome and the proteome of
                      the SAME tumour type agree more than crossed pairs do.
                      Answers "does splitting make the layers tell a more
                      coherent story".

Permutation is by exhaustive enumeration of the tumour-type labellings, not
random sampling, because the number of cohorts is small enough to enumerate --
and small enough that the attainable p-value is bounded. With 4 + 4 proteome
cohorts there are C(8,4)/2 = 35 distinct labellings, so p cannot go below
0.029; with 2 + 3 genome cohorts there are 10, so p cannot go below 0.10. The
genome layer therefore reports effect sizes without a significance claim.
"""
from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr


# ----------------------------------------------------------------- profiles
def enrichment_profile(ora: pd.DataFrame, layer: str, scheme: str,
                       collections: Sequence[str] | None = None,
                       min_cohorts: int = 2) -> pd.DataFrame:
    """Pathways x cohorts of a signed enrichment score.

    score = -log10(p_up) - (-log10(p_down)), so a pathway enriched among
    up-regulated genes is positive, among down-regulated negative, and one
    enriched in both directions sits near zero. Using the signed score means
    two cohorts that agree a pathway is involved but disagree on the direction
    are far apart, which is what "the same biology" has to mean here.
    """
    o = ora[(ora["layer"] == layer) & (ora["scheme"] == scheme)]
    if collections:
        o = o[o["collection"].isin(list(collections))]
    key = o["collection"] + " | " + o["pathway"]
    o = o.assign(key=key)

    up = o[o["direction"] == "up"].pivot_table(index="key", columns="dataset",
                                               values="pval", aggfunc="min")
    dn = o[o["direction"] == "down"].pivot_table(index="key", columns="dataset",
                                                 values="pval", aggfunc="min")
    idx = up.index.union(dn.index)
    cols = sorted(set(up.columns) | set(dn.columns))
    up = up.reindex(index=idx, columns=cols).fillna(1.0)
    dn = dn.reindex(index=idx, columns=cols).fillna(1.0)
    prof = -np.log10(up.clip(lower=1e-300)) + np.log10(dn.clip(lower=1e-300))
    # a pathway must be testable in enough cohorts to compare them on it
    tested = (up < 1.0) | (dn < 1.0)
    prof = prof[tested.sum(axis=1) >= min_cohorts]
    return prof


def profile_distance(prof: pd.DataFrame) -> pd.DataFrame:
    """1 - Spearman rho between cohort profiles."""
    cols = list(prof.columns)
    rho = np.ones((len(cols), len(cols)))
    for i, a in enumerate(cols):
        for j in range(i + 1, len(cols)):
            r = spearmanr(prof[a], prof[cols[j]]).statistic
            r = 0.0 if np.isnan(r) else r
            rho[i, j] = rho[j, i] = r
    return pd.DataFrame(1.0 - rho, index=cols, columns=cols)


# ------------------------------------------------------------------- tests
def _label_splits(n: int, k: int, symmetric: bool) -> List[tuple]:
    """Every way of choosing which k of n objects carry the first label."""
    all_ = list(combinations(range(n), k))
    if symmetric and 2 * k == n:
        # swapping the two labels gives the same partition; keep one of each
        seen, out = set(), []
        for c in all_:
            comp = tuple(sorted(set(range(n)) - set(c)))
            if comp in seen:
                continue
            seen.add(c)
            out.append(c)
        return out
    return all_


def anosim(D: pd.DataFrame, labels: Sequence[str]) -> dict:
    """ANOSIM R with an exact permutation p-value.

    R = (mean rank of between-group distances - mean rank of within-group
    distances) / (M/2), M = number of pairs. R > 0 means cohorts of the same
    type are closer; R near 0 means the labels explain nothing.
    """
    labels = np.asarray(labels)
    n = len(labels)
    d = squareform(D.values, checks=False)
    ranks = rankdata(d)
    pairs = np.array(list(combinations(range(n), 2)))

    def stat(lab):
        same = lab[pairs[:, 0]] == lab[pairs[:, 1]]
        m = len(ranks)
        return (ranks[~same].mean() - ranks[same].mean()) / (m / 2)

    obs = stat(labels)
    uniq = sorted(set(labels))
    k = int((labels == uniq[0]).sum())
    splits = _label_splits(n, k, symmetric=True)
    null = []
    for c in splits:
        lab = np.array([uniq[1]] * n, dtype=object)
        lab[list(c)] = uniq[0]
        null.append(stat(lab))
    null = np.array(null)
    return {"R": float(obs), "p": float((null >= obs).sum() / len(null)),
            "n_permutations": len(null), "min_attainable_p": 1.0 / len(null)}


def permanova(D: pd.DataFrame, labels: Sequence[str]) -> dict:
    """PERMANOVA pseudo-F and R2 with an exact permutation p-value.

    R2 is the share of total multivariate variation the labels explain -- the
    number that decides whether splitting is worth the loss of cohorts. A
    significant p with a small R2 means the split is real but not useful.
    """
    labels = np.asarray(labels)
    n = len(labels)
    d2 = D.values ** 2
    ss_total = d2[np.triu_indices(n, 1)].sum() / n

    def ss_within(lab):
        s = 0.0
        for g in set(lab):
            idx = np.where(lab == g)[0]
            if len(idx) < 2:
                continue
            sub = d2[np.ix_(idx, idx)]
            s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s

    def stat(lab):
        ssw = ss_within(lab)
        ssa = ss_total - ssw
        a = len(set(lab))
        f = (ssa / (a - 1)) / (ssw / (n - a)) if ssw > 0 else np.inf
        return f, ssa / ss_total if ss_total > 0 else np.nan

    obs_f, obs_r2 = stat(labels)
    uniq = sorted(set(labels))
    k = int((labels == uniq[0]).sum())
    null = []
    for c in _label_splits(n, k, symmetric=True):
        lab = np.array([uniq[1]] * n, dtype=object)
        lab[list(c)] = uniq[0]
        null.append(stat(lab)[0])
    null = np.array(null)
    return {"pseudo_F": float(obs_f), "R2": float(obs_r2),
            "p": float((null >= obs_f).sum() / len(null)),
            "n_permutations": len(null), "min_attainable_p": 1.0 / len(null)}


def silhouette(D: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    """Per-cohort silhouette against the tumour-type labels.

    Negative means that cohort sits closer to the other tumour type than to its
    own -- the most direct way to see which cohorts, if any, drive a split.
    """
    labels = np.asarray(labels)
    out = []
    for i, c in enumerate(D.index):
        own = labels == labels[i]
        own[i] = False
        if own.sum() == 0:
            out.append({"cohort": c, "tissue": labels[i], "silhouette": np.nan})
            continue
        a = D.values[i, own].mean()
        b = min(D.values[i, labels == g].mean() for g in set(labels) if g != labels[i])
        out.append({"cohort": c, "tissue": labels[i],
                    "silhouette": (b - a) / max(a, b)})
    return pd.DataFrame(out)


def pcoa(D: pd.DataFrame, k: int = 2) -> pd.DataFrame:
    """Classical multidimensional scaling of the distance matrix."""
    n = len(D)
    A = -0.5 * D.values ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    G = J @ A @ J
    w, v = np.linalg.eigh(G)
    order = np.argsort(w)[::-1][:k]
    coords = v[:, order] * np.sqrt(np.clip(w[order], 0, None))
    var = np.clip(w, 0, None)
    return pd.DataFrame(coords, index=D.index,
                        columns=[f"PCo{i+1}" for i in range(k)]).assign(
        **{f"PCo{i+1}_var": var[o] / var.sum() for i, o in enumerate(order)})
