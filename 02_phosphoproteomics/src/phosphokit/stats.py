"""limma tests.

Design notes that apply to every dataset here:

* Values are tumor/normal comparisons, so the question is always whether a contrast differs from 0.
  Two shapes occur: ratio data with one value per patient (intercept-only model) and intensity data
  with separate tumor and normal samples (two-group model).
* Patients often contribute several samples. Samples of the same patient in the same group are
  averaged first: several lesions of one patient are not independent replicates of that group.
  Across groups, a patient can appear twice, which duplicateCorrelation(block = patient) handles.
* eBayes shrinks each feature's variance toward a prior fitted on all features, which matters when
  many phosphosites are quantified in only 10-20 samples.
"""
import numpy as np
import pandas as pd
from scipy import stats as sps

from .config import DEFAULTS
from .rbridge import limma, to_r, r_list_to_dict, ro

_FIT = ro.r("""
function(Y, design, block) {
    rho <- NA_real_
    rho_gene <- numeric(0)
    if (!is.null(block) && any(duplicated(block))) {
        dc <- duplicateCorrelation(Y, design, block = block)
        rho <- dc$consensus.correlation
        rho_gene <- tanh(dc$atanh.correlations)
        fit <- lmFit(Y, design, block = block, correlation = rho)
    } else {
        fit <- lmFit(Y, design)
    }
    fit <- eBayes(fit)
    k <- ncol(design)
    list(logFC = fit$coefficients[, k],
         SE = sqrt(fit$s2.post) * fit$stdev.unscaled[, k],
         t = fit$t[, k], df = fit$df.total, P = fit$p.value[, k],
         rho = rho, rho_gene = rho_gene, d0 = fit$df.prior, s0 = fit$s2.prior)
}
""")


def bh_adjust(p):
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    adj = np.empty_like(pv)
    adj[order] = np.minimum(ranked, 1)
    out[ok] = adj
    return out


def limma_test(mat, design, block=None, keep=None, extra=None):
    """
    Fit `design` (n_samples x n_coefficients) with limma and test its LAST coefficient.

    `block` gives the patient of each sample; duplicateCorrelation uses it when a patient appears
    more than once. `keep` is a boolean Series selecting the features to test - applying it before
    eBayes keeps the variance prior from being learned on features that are never reported.
    Returns (results, info).
    """
    Y = mat if keep is None else mat[keep]
    r_block = ro.StrVector(list(block)) if block is not None else ro.NULL
    out = r_list_to_dict(_FIT(to_r(Y.values), to_r(np.asarray(design, float)), r_block))

    res = pd.DataFrame({k: out[k] for k in ["logFC", "SE", "t", "df", "P"]}, index=Y.index)
    if extra is not None:
        for name, values in extra.items():
            res[name] = np.asarray(values)[keep.values] if keep is not None else np.asarray(values)
    res = res.dropna(subset=["P"])
    res["adj.P"] = bh_adjust(res["P"])
    info = {"rho": float(out["rho"][0]), "rho_per_feature": out["rho_gene"],
            "d0": float(np.ravel(out["d0"])[0]), "s0": float(np.ravel(out["s0"])[0]),
            "n_tested": len(res)}
    return res, info


def one_sample_test(mat, settings=DEFAULTS, block=None):
    """
    Intercept-only model for data that are already log2 tumor/normal ratios.
    H0: mean log ratio = 0.
    """
    n = mat.notna().sum(axis=1)
    keep = n >= settings.min_per_group
    design = np.ones((mat.shape[1], 1))
    return limma_test(mat, design, block=block, keep=keep, extra={"n": n})


def collapse_patients(mat, meta, how="mean"):
    """
    Average the samples of one patient within one group.

    meta: DataFrame indexed by sample with columns Case and Group. Returns the collapsed matrix
    and a unit table (Case, Group) with one column per patient and group.
    """
    meta = meta.loc[[s for s in mat.columns if s in meta.index]]
    key = meta["Case"].astype(str) + "|" + meta["Group"].astype(str)
    units = mat[meta.index].T.groupby(key.values).agg(how).T
    info = pd.DataFrame({"Case": units.columns.str.split("|").str[0],
                         "Group": units.columns.str.split("|").str[1]}, index=units.columns)
    return units, info


def two_group_test(mat, meta, test, ref, settings=DEFAULTS, collapse=True, how="mean"):
    """
    Test `test` groups against `ref` groups (logFC > 0 = higher in `test`).

    With collapse=True the samples of one patient in one group are averaged first. Patients present
    in both arms are blocked, so the paired part of the design is used without dropping the
    patients that appear in only one arm.
    """
    meta = meta[meta["Group"].isin(list(test) + list(ref))]
    mat = mat[[c for c in mat.columns if c in meta.index]]
    meta = meta.loc[mat.columns]
    arm = pd.Series(np.where(meta["Group"].isin(test), "test", "ref"), index=meta.index)

    if collapse:
        units, unit_info = collapse_patients(mat, meta.assign(Group=arm), how=how)
        arm = unit_info["Group"]
        block = unit_info["Case"]
    else:
        units, block = mat, meta["Case"]

    x = (arm == "test").astype(float).values
    n_ref = units.loc[:, x == 0].notna().sum(axis=1)
    n_test = units.loc[:, x == 1].notna().sum(axis=1)
    keep = (n_ref >= settings.min_per_group) & (n_test >= settings.min_per_group)

    design = np.column_stack([np.ones(len(x)), x])
    res, info = limma_test(units, design, block=block, keep=keep,
                           extra={"n_ref": n_ref, "n_test": n_test})
    info["patients_both_arms"] = int((block.groupby(block.values).size() > 1).sum())
    return res, info


def is_significant(res, settings=DEFAULTS):
    p = res["adj.P"] if settings.use_adj_p else res["P"]
    return (p < settings.p_cut) & (res["logFC"].abs() > settings.lfc_cut)


def direction(res, settings=DEFAULTS):
    """Up / Down / '' (tested but not significant)."""
    sig = is_significant(res, settings)
    return pd.Series(np.where(~sig, "", np.where(res["logFC"] > 0, "Up", "Down")), index=res.index)


def welch_satterthwaite(se1, se2, df1, df2):
    return (se1 ** 2 + se2 ** 2) ** 2 / (se1 ** 4 / df1 + se2 ** 4 / df2)


def t_pvalue(logfc, se, df):
    return 2 * sps.t.sf(np.abs(logfc / se), df)
