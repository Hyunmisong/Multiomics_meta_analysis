"""The subproject's actual question: does dividing phosphosite ratios by the
parent protein change which pathways come out enriched?

Three levels of comparison:

1. site level  -- only possible where both arms ship the full table (PG4).
                  How many sites change call, how many flip direction, how much
                  does the effect size shrink.
2. gene level  -- overlap of the regulated gene sets per cohort.
3. pathway level -- overlap and rank agreement of the enriched pathway sets,
                  plus which pathways are gained or lost on normalisation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return np.nan
    return len(a & b) / len(a | b)


def site_level(long_scheme: pd.DataFrame) -> pd.DataFrame:
    """Paired per-site comparison for cohorts with full tables in both arms."""
    ph = long_scheme[long_scheme["layer"] == "phospho"]
    rows = []
    for (dataset, scheme), sub in ph.groupby(["dataset", "scheme"]):
        arms = set(sub["norm"])
        if not {"raw", "prot_adj"} <= arms:
            continue
        raw = sub[sub["norm"] == "raw"].set_index("feature_id")
        adj = sub[sub["norm"] == "prot_adj"].set_index("feature_id")
        common = raw.index.intersection(adj.index)
        # Fewer than ~1000 shared sites means at least one arm is a filtered
        # subset; the paired statistics would be conditioned on significance.
        paired_ok = len(common) >= 1000
        r, a = raw.loc[common], adj.loc[common]
        both = r["is_sig"] & a["is_sig"]
        flip = both & (r["direction"] != a["direction"])
        fc_r, fc_a = r["logFC"], a["logFC"]
        ok = fc_r.notna() & fc_a.notna()
        try:
            w_p = wilcoxon(fc_r[ok].abs(), fc_a[ok].abs()).pvalue if ok.sum() > 10 else np.nan
        except ValueError:
            w_p = np.nan
        rows.append({
            "dataset": dataset, "scheme": scheme,
            "paired_comparison_valid": paired_ok,
            "n_sites_raw": len(raw), "n_sites_adj": len(adj), "n_sites_common": len(common),
            "n_sig_raw": int(r["is_sig"].sum()), "n_sig_adj": int(a["is_sig"].sum()),
            "n_sig_both": int(both.sum()),
            "n_sig_raw_only": int((r["is_sig"] & ~a["is_sig"]).sum()),
            "n_sig_adj_only": int((~r["is_sig"] & a["is_sig"]).sum()),
            "n_direction_flips": int(flip.sum()),
            "median_absFC_raw": float(fc_r[ok].abs().median()),
            "median_absFC_adj": float(fc_a[ok].abs().median()),
            "absFC_shrinkage": float(fc_r[ok].abs().median() - fc_a[ok].abs().median()),
            "spearman_logFC": float(spearmanr(fc_r[ok], fc_a[ok]).statistic) if ok.sum() > 10 else np.nan,
            "wilcoxon_absFC_p": w_p,
        })
    return pd.DataFrame(rows)


def gene_level(gene_calls: pd.DataFrame) -> pd.DataFrame:
    ph = gene_calls[gene_calls["layer"] == "phospho"]
    rows = []
    for (dataset, scheme), sub in ph.groupby(["dataset", "scheme"]):
        for d in ["any", "up", "down"]:
            s = sub if d == "any" else sub[sub["direction"] == d]
            raw = set(s.loc[s["norm"] == "raw", "gene"])
            adj = set(s.loc[s["norm"] == "prot_adj", "gene"])
            rows.append({"dataset": dataset, "scheme": scheme, "direction": d,
                         "n_genes_raw": len(raw), "n_genes_adj": len(adj),
                         "n_shared": len(raw & adj), "jaccard": _jaccard(raw, adj),
                         "raw_only": len(raw - adj), "adj_only": len(adj - raw)})
    return pd.DataFrame(rows)


def pathway_level(ora: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Per cohort x collection x scheme x direction: how much of the enriched
    pathway set survives protein normalisation."""
    ph = ora[ora["layer"] == "phospho"]
    rows = []
    keys = ["dataset", "collection", "scheme", "direction"]
    for key, sub in ph.groupby(keys):
        raw = sub[sub["norm"] == "raw"]
        adj = sub[sub["norm"] == "prot_adj"]
        if raw.empty or adj.empty:
            continue
        sraw = set(raw.loc[raw["padj"] <= alpha, "pathway"])
        sadj = set(adj.loc[adj["padj"] <= alpha, "pathway"])
        m = raw.merge(adj, on="pathway", suffixes=("_raw", "_adj"))
        rho = (spearmanr(m["neglog10_padj_raw"], m["neglog10_padj_adj"]).statistic
               if len(m) > 10 else np.nan)
        rows.append(dict(zip(keys, key)) | {
            "n_enriched_raw": len(sraw), "n_enriched_adj": len(sadj),
            "n_shared": len(sraw & sadj), "jaccard": _jaccard(sraw, sadj),
            "lost_on_adjust": len(sraw - sadj), "gained_on_adjust": len(sadj - sraw),
            "spearman_neglog10padj": rho, "n_pathways_compared": len(m),
            "lost_examples": ";".join(sorted(sraw - sadj)[:8]),
            "gained_examples": ";".join(sorted(sadj - sraw)[:8]),
        })
    return pd.DataFrame(rows)


def scheme_stability(ora: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Which threshold scheme gives the most reproducible pathway calls?

    Reported per layer x collection x direction:
      - mean pairwise Jaccard of enriched pathway sets across datasets
        (higher = cohorts agree more)
      - median number of enriched pathways per arm (a scheme calling everything
        enriched scores a high Jaccard for the wrong reason, so read both)
    """
    rows = []
    for (layer, collection, scheme, direction), sub in ora.groupby(
            ["layer", "collection", "scheme", "direction"]):
        sets_by_arm = {arm: set(s.loc[s["padj"] <= alpha, "pathway"])
                       for arm, s in sub.groupby("arm")}
        sets_by_arm = {k: v for k, v in sets_by_arm.items() if v}
        arms = sorted(sets_by_arm)
        js = [_jaccard(sets_by_arm[a], sets_by_arm[b])
              for i, a in enumerate(arms) for b in arms[i + 1:]]
        rows.append({
            "layer": layer, "collection": collection, "scheme": scheme,
            "direction": direction, "n_arms_with_hits": len(arms),
            "median_n_enriched": float(np.median([len(v) for v in sets_by_arm.values()]))
            if sets_by_arm else 0.0,
            "mean_pairwise_jaccard": float(np.mean(js)) if js else np.nan,
        })
    return pd.DataFrame(rows)
