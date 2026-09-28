"""Normalization of log2 matrices (features x samples)."""
import numpy as np
import pandas as pd

from .config import DEFAULTS
from .rbridge import limma, to_r


def quantile_normalize(lg):
    """limma::normalizeQuantiles; it handles missing values by interpolation."""
    out = limma.normalizeQuantiles(to_r(lg.values.astype(float)))
    return pd.DataFrame(np.array(out), index=lg.index, columns=lg.columns)


def normalize(lg, method=DEFAULTS.normalization, common_frac=DEFAULTS.norm_common_frac):
    """
    method:
      "quantile_common" (default) quantile-normalize the features observed in >= common_frac of
                        samples, then shift the remaining features by the same per-sample offset
      "quantile"        plain quantile normalization on all features
      "median_common"   per-sample median centring on the common features
      "median"          per-sample median centring on all observed values
      "none"

    The "_common" variants exist because detection depth differs between groups; see Settings.
    """
    common = lg.notna().mean(axis=1) >= common_frac
    if method == "none":
        return lg.copy()
    if method == "median":
        return lg - lg.median()
    if method == "median_common":
        return lg - lg[common].median()
    if method == "quantile":
        return quantile_normalize(lg)
    if method == "quantile_common":
        q = quantile_normalize(lg[common])
        out = lg + (q.median() - lg[common].median())
        out.loc[common] = q
        return out
    raise ValueError(method)


def normalization_bias(lg, groups, reference, common_frac=DEFAULTS.norm_common_frac,
                       methods=("none", "median", "quantile", "median_common", "quantile_common")):
    """
    Median log2FC of the common features for each group against `reference`, per method.

    Most features should not change between groups, so a good normalization leaves this near 0.
    Use it to justify the choice of `normalization` rather than taking it on faith.
    """
    common = lg.index[lg.notna().mean(axis=1) >= common_frac]
    ref_cols = groups.index[groups == reference]
    rows = []
    for m in methods:
        nm = normalize(lg, m, common_frac)
        for g in [g for g in pd.unique(groups) if g != reference]:
            diff = (nm.loc[common, groups.index[groups == g]].mean(axis=1)
                    - nm.loc[common, ref_cols].mean(axis=1))
            rows.append({"method": m, "group": g, "median_log2FC": diff.median()})
    return pd.DataFrame(rows).pivot(index="method", columns="group", values="median_log2FC")
