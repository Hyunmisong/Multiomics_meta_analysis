"""Missing-value handling."""
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from .config import DEFAULTS


def impute(mat, method=DEFAULTS.missing_method, min_valid_frac=DEFAULTS.impute_min_valid_frac,
           settings=DEFAULTS):
    """
    Impute a log2 matrix. Only features observed in >= min_valid_frac of samples are filled;
    the rest keep their NaNs and are removed later by the per-group filter.

    "none"      no imputation (limma tests each feature on the samples where it was observed)
    "downshift" Perseus-style: draw from a narrow normal below each sample's observed distribution,
                which assumes missing = below the detection limit
    "half_min"  feature minimum minus 1 on log2 scale
    "knn"       neighbours are features with a similar profile across samples

    Returns (imputed matrix, mask of imputed cells).
    """
    empty_mask = pd.DataFrame(False, index=mat.index, columns=mat.columns)
    if method == "none":
        return mat.copy(), empty_mask

    eligible = mat.notna().mean(axis=1) >= min_valid_frac
    sub = mat[eligible].copy()
    rng = np.random.default_rng(settings.seed)

    if method == "downshift":
        for c in sub.columns:
            obs = sub[c].dropna()
            miss = sub[c].isna()
            sub.loc[miss, c] = rng.normal(obs.mean() - settings.downshift * obs.std(),
                                          settings.downshift_width * obs.std(), miss.sum())
    elif method == "half_min":
        fill = sub.min(axis=1) - 1
        sub = sub.apply(lambda col: col.fillna(fill))
    elif method == "knn":
        sub = pd.DataFrame(KNNImputer(n_neighbors=settings.knn_k).fit_transform(sub.values),
                           index=sub.index, columns=sub.columns)
    else:
        raise ValueError(method)

    out = mat.copy()
    out.loc[eligible] = sub
    return out, mat.isna() & out.notna()
