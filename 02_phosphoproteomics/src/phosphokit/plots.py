"""QC and result figures. Every function draws on an axis or returns a figure; saving is left to
the notebook so that figure paths stay dataset-specific."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import DEFAULTS
from .stats import direction

UP, DOWN, NS = "#C44E52", "#4C72B0", "lightgray"


def sample_boxplot(mat, ax, groups=None, colors=None, title=""):
    """Per-sample distribution of log2 values; groups colour the boxes."""
    order = mat.columns if groups is None else groups.loc[mat.columns].sort_values().index
    bp = ax.boxplot([mat[c].dropna() for c in order], showfliers=False, widths=0.7,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, col in zip(bp["boxes"], order):
        patch.set_facecolor("gray" if groups is None else colors[groups[col]])
        patch.set_edgecolor("none")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks([])
    ax.set_title(title)
    ax.set_ylabel("log2")


def phospho_protein_correlation(phospho_log, protein_by_site, ax, min_pairs=10):
    """
    Per-site correlation between phospho and protein across samples.

    A distribution clearly shifted to the right is the reason the adjustment is needed at all, and
    the reason Method B is conservative (it ignores this correlation).
    """
    both = phospho_log.notna() & protein_by_site.notna()
    ok = both.index[both.sum(axis=1) >= min_pairs]
    r = phospho_log.loc[ok].T.corrwith(protein_by_site.loc[ok].T)
    ax.hist(r.dropna(), bins=50, color="gray")
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(r.median(), color="red", ls="--", label=f"median r = {r.median():.2f}")
    ax.set_xlabel(f"Pearson r across samples (>= {min_pairs} paired values)")
    ax.set_title("Phospho vs protein correlation")
    ax.legend()
    return r


def volcano(res, ax, title="", labels=None, n_label=10, settings=DEFAULTS):
    p = res["adj.P"] if settings.use_adj_p else res["P"]
    y = -np.log10(p)
    d = direction(res, settings)
    for key, color, size in [("", NS, 4), ("Down", DOWN, 9), ("Up", UP, 9)]:
        m = d == key
        ax.scatter(res.loc[m, "logFC"], y[m], c=color, s=size, alpha=0.7, lw=0)
    ax.axhline(-np.log10(settings.p_cut), color="black", ls="--", lw=0.6)
    for v in (-settings.lfc_cut, settings.lfc_cut):
        ax.axvline(v, color="black", ls="--", lw=0.6)
    for idx, row in res[d != ""].assign(_y=y).nlargest(n_label, "_y").iterrows():
        ax.annotate(labels.get(idx, idx) if labels is not None else idx,
                    (row["logFC"], row["_y"]), fontsize=6, xytext=(3, 2), textcoords="offset points")
    ax.set_title(f"{title}\nUp {(d == 'Up').sum()} / Down {(d == 'Down').sum()} (n = {len(res)})",
                 fontsize=10)
    ax.set_xlabel("log2FC")
    ax.set_ylabel(f"-log10 {'adj. ' if settings.use_adj_p else ''}p")


def save_volcanoes(results, out_dir, labels=None, formats=("png", "pdf"), n_label=15,
                   settings=DEFAULTS, figsize=(6, 5.5)):
    """One file per analysis, large enough for a figure panel. `results`: {name: table}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, res in results.items():
        fig, ax = plt.subplots(figsize=figsize)
        lab = labels.get(name) if isinstance(labels, dict) else labels
        volcano(res, ax, title=name, labels=lab, n_label=n_label, settings=settings)
        fig.tight_layout()
        for fmt in formats:
            fig.savefig(out_dir / f"volcano_{name}.{fmt}", bbox_inches="tight")
        plt.close(fig)


def pvalue_histograms(results, figsize=None):
    """Flat with a spike near 0 is healthy; a hump near 1 points to a model or normalization problem."""
    fig, axes = plt.subplots(1, len(results), figsize=figsize or (4 * len(results), 3.2), squeeze=False)
    for ax, (name, res) in zip(axes[0], results.items()):
        ax.hist(res["P"], bins=40, color="gray")
        ax.set_title(f"{name} p-values", fontsize=10)
        ax.set_xlabel("p")
    return fig


def method_comparison(res_a, res_b, ax_logfc, ax_se, site_correlation=None):
    """Method A vs Method B: the estimates agree, the standard errors do not."""
    cmp = res_a[["logFC", "SE"]].join(res_b[["logFC", "SE"]], lsuffix="_A", rsuffix="_B", how="inner")
    ax_logfc.scatter(cmp["logFC_A"], cmp["logFC_B"], s=3, alpha=0.2)
    lim = np.nanpercentile(np.abs(cmp[["logFC_A", "logFC_B"]].values), 99.9)
    ax_logfc.plot([-lim, lim], [-lim, lim], color="red", lw=0.8)
    ax_logfc.set_xlabel("Method A")
    ax_logfc.set_ylabel("Method B")
    ax_logfc.set_title(f"Adjusted logFC (r = {cmp['logFC_A'].corr(cmp['logFC_B']):.3f})", fontsize=10)

    ratio = np.log2(cmp["SE_B"] / cmp["SE_A"])
    if site_correlation is not None:
        ax_se.scatter(site_correlation.reindex(cmp.index), ratio, s=3, alpha=0.2)
        ax_se.set_xlabel("phospho-protein r")
    else:
        ax_se.hist(ratio.dropna(), bins=50, color="gray")
        ax_se.set_xlabel("log2(SE_B / SE_A)")
    ax_se.axhline(0, color="red", lw=0.8) if site_correlation is not None else None
    ax_se.set_title("SE inflation of Method B", fontsize=10)
    return cmp


def embedding_plot(mat, groups, colors, ax_pca=None, methods=("PCA", "UMAP"), seed=DEFAULTS.seed,
                   min_valid=0.7, fill="median", perplexity=30, title_prefix=""):
    """
    PCA / UMAP / t-SNE of the samples.

    Only fairly complete features are used and the remaining gaps are filled here, for the plot
    only, so the picture does not depend on the imputation setting of the analysis.
    fill="median" (row median) suits adjusted values, which are differences rather than intensities;
    fill="downshift" suits raw intensities.
    """
    from sklearn.decomposition import PCA

    sub = mat[mat.notna().mean(axis=1) >= min_valid]
    if fill == "median":
        filled = sub.apply(lambda row: row.fillna(row.median()), axis=1)
    else:
        from .impute import impute
        filled, _ = impute(sub, method="downshift", min_valid_frac=0)
    X = filled.T.values
    X = (X - X.mean(axis=0)) / X.std(axis=0)

    pca = PCA(n_components=min(10, X.shape[0] - 1), random_state=seed)
    pcs = pca.fit_transform(X)
    embeddings = {"PCA": (pcs[:, :2], pca.explained_variance_ratio_[:2])}
    if "UMAP" in methods:
        try:
            import umap
            embeddings["UMAP"] = (umap.UMAP(n_neighbors=15, min_dist=0.3,
                                            random_state=seed).fit_transform(pcs), None)
        except ImportError:
            print("umap-learn not installed; skipping UMAP")
    if "t-SNE" in methods:
        from sklearn.manifold import TSNE
        embeddings["t-SNE"] = (TSNE(n_components=2, perplexity=perplexity, init="pca",
                                    random_state=seed).fit_transform(pcs), None)

    fig, axes = plt.subplots(1, len(embeddings), figsize=(5 * len(embeddings), 4.3), squeeze=False)
    point_colors = groups.loc[filled.columns].map(colors)
    for ax, (name, (xy, var)) in zip(axes[0], embeddings.items()):
        ax.scatter(xy[:, 0], xy[:, 1], c=point_colors, s=20)
        ax.set_title(f"{title_prefix}{name} ({sub.shape[0]} features)", fontsize=10)
        if var is not None:
            ax.set_xlabel(f"PC1 ({var[0]:.1%})")
            ax.set_ylabel(f"PC2 ({var[1]:.1%})")
        else:
            # axis values of UMAP and t-SNE carry no meaning
            ax.set_xticks([]); ax.set_yticks([])
    for g, c in colors.items():
        axes[0][0].scatter([], [], color=c, label=g)
    axes[0][0].legend(fontsize=8)
    return fig, pcs
