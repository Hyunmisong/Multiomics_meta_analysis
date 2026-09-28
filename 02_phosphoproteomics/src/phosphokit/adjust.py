"""Removing the protein-abundance component from phosphosite signals.

A phosphosite signal rises either because the site is phosphorylated more often or because there is
more of the protein. Two ways to separate them:

Method A (primary)
    Subtract the protein value from the phospho value in every sample, then test the differences.
    Because both values come from the same sample, the within-sample covariance is handled
    automatically: Var(X - Y) = Var(X) + Var(Y) - 2*Cov(X, Y).

Method B (sensitivity, MSstatsPTM-style)
    Test phospho and protein separately and subtract the estimates, combining the standard errors
    as sqrt(SE_p^2 + SE_prot^2). That drops the covariance term, so with positively correlated
    layers (median r was about 0.4-0.5 in these datasets) it overestimates the SE and is
    conservative. It is the right choice when the two layers come from different aliquots and
    cannot be paired, and it uses every sample where the protein was measured, which helps for
    sites seen in few samples.

Both are reported: they answer the same question with different assumptions, and agreement between
them is a useful robustness flag.
"""
import numpy as np
import pandas as pd

from .config import DEFAULTS
from .stats import bh_adjust, direction, t_pvalue, welch_satterthwaite

CATEGORY_ORDER = ["Phospho-specific", "Both, same direction", "Both, opposite direction",
                  "Explained by protein", "Protein only", "Not significant"]
CATEGORY_COLORS = dict(zip(CATEGORY_ORDER,
                           ["#C44E52", "#8172B3", "#DD8452", "#4C72B0", "#55A868", "lightgray"]))


def adjusted_matrix(phospho_log, protein_log, site_to_protein):
    """
    log2(phospho) - log2(protein of the same gene/accession) per sample.

    site_to_protein: Series mapping each phospho feature to a row of protein_log. Sites whose
    protein was not quantified are dropped here and keep only an unadjusted result.
    """
    keep = site_to_protein[site_to_protein.isin(protein_log.index)]
    protein_by_site = protein_log.loc[keep.values].set_axis(keep.index)
    common = [c for c in phospho_log.columns if c in protein_by_site.columns]
    return phospho_log.loc[keep.index, common] - protein_by_site[common]


def method_b(res_phospho, res_protein, site_to_protein):
    """Combine separate phospho and protein fits (see module docstring)."""
    b = res_phospho[["logFC", "SE", "df"]].join(site_to_protein.rename("_protein"))
    b = b.join(res_protein[["logFC", "SE", "df"]], on="_protein", rsuffix="_prot", how="inner")
    out = pd.DataFrame(index=b.index)
    out["logFC"] = b["logFC"] - b["logFC_prot"]
    out["SE"] = np.sqrt(b["SE"] ** 2 + b["SE_prot"] ** 2)
    out["df"] = welch_satterthwaite(b["SE"], b["SE_prot"], b["df"], b["df_prot"])
    out["t"] = out["logFC"] / out["SE"]
    out["P"] = t_pvalue(out["logFC"], out["SE"], out["df"])
    out["adj.P"] = bh_adjust(out["P"])
    return out


def categorize(protein_dir, phospho_dir, adjusted_dir):
    """One category per site; `adjusted_dir` should come from Method A."""
    if adjusted_dir and not protein_dir:
        return "Phospho-specific"
    if adjusted_dir and protein_dir:
        return "Both, same direction" if adjusted_dir == protein_dir else "Both, opposite direction"
    if phospho_dir:
        return "Explained by protein"
    if protein_dir:
        return "Protein only"
    return "Not significant"


def category_table(res_protein, res_phospho, res_adjusted, site_to_protein, settings=DEFAULTS):
    """Site-level categories for one comparison."""
    out = pd.DataFrame(index=res_adjusted.index)
    out["protein"] = site_to_protein.reindex(out.index).map(direction(res_protein, settings)).fillna("")
    out["phospho"] = direction(res_phospho, settings).reindex(out.index).fillna("")
    out["adjusted"] = direction(res_adjusted, settings)
    out["category"] = [categorize(*row) for row in out[["protein", "phospho", "adjusted"]].values]
    return out
