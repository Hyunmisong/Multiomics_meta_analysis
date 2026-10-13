"""Summary panels for the phospho-normalisation question.

  A  significant phosphosites per cohort, raw vs protein-adjusted
  B  effect-size shrinkage on the paired cohort (PG4)
  C  agreement of the enriched pathway sets (Jaccard), raw vs protein-adjusted
  D  how many pathways concentrate across cohorts, by layer and threshold scheme

Output: figures/norm_summary.png / .pdf
"""
import _bootstrap  # noqa: F401

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pmap.paths import FIGURES, RESULTS

RAW, ADJ = "#c0392b", "#2a6fb0"
SCHEME = "fdr05"


def main():
    site = pd.read_csv(RESULTS / "norm_site_level.csv")
    pw = pd.read_csv(RESULTS / "norm_pathway_level.csv")
    conc = pd.read_csv(RESULTS / "concentration.csv")
    long = pd.read_parquet(RESULTS / "features_long.parquet")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    fig.suptitle("Does protein-level normalisation change the phospho pathway answer?",
                 fontsize=13, fontweight="bold", x=0.01, ha="left")

    # ---- A: significant sites per cohort -------------------------------
    ax = axes[0, 0]
    s = site[site["scheme"] == SCHEME].sort_values("dataset")
    x = np.arange(len(s))
    ax.bar(x - 0.19, s["n_sig_raw"], 0.38, color=RAW, label="raw")
    ax.bar(x + 0.19, s["n_sig_adj"], 0.38, color=ADJ, label="protein-adjusted")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n({'paired' if v else 'subset only'})"
                        for d, v in zip(s["dataset"], s["paired_comparison_valid"])])
    ax.set_ylabel("significant phosphosites")
    ax.set_title("A  Significant sites per cohort", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)

    # ---- B: paired effect sizes, PG4 ------------------------------------
    ax = axes[0, 1]
    pg4 = long[(long["layer"] == "phospho") & (long["dataset"] == "PG4")
               & (long["scheme"] == SCHEME)]
    r = pg4[pg4["norm"] == "raw"].set_index("feature_id")["logFC"]
    a = pg4[pg4["norm"] == "prot_adj"].set_index("feature_id")["logFC"]
    common = r.index.intersection(a.index)
    ax.scatter(r.loc[common], a.loc[common], s=2, alpha=0.18, color="#4a5568", linewidths=0)
    lim = float(np.nanmax(np.abs([r.loc[common].min(), r.loc[common].max()]))) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="#999", lw=0.8, ls="--")
    ax.axhline(0, color="#ccc", lw=0.6)
    ax.axvline(0, color="#ccc", lw=0.6)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC, raw phospho")
    ax.set_ylabel("log2FC, protein-adjusted")
    rho = s.loc[s["dataset"] == "PG4", "spearman_logFC"].iloc[0]
    ax.set_title(f"B  PG4 paired sites (n={len(common):,}, Spearman {rho:.2f})",
                 loc="left", fontsize=11, fontweight="bold")

    # ---- C: pathway-set agreement ---------------------------------------
    ax = axes[1, 0]
    p = pw[(pw["scheme"] == SCHEME) & (pw["direction"] != "any")].copy()
    p["grp"] = p["dataset"] + "  " + p["collection"] + "  " + p["direction"]
    p = p.sort_values(["dataset", "collection", "direction"])
    ax.barh(np.arange(len(p)), p["jaccard"].fillna(0), color="#5b6472")
    ax.set_yticks(np.arange(len(p)))
    ax.set_yticklabels(p["grp"], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Jaccard of enriched pathway sets (raw vs protein-adjusted)")
    ax.set_xlim(0, 1)
    ax.set_title("C  How much of the pathway answer survives normalisation",
                 loc="left", fontsize=11, fontweight="bold")

    # ---- D: cross-cohort concentration -----------------------------------
    ax = axes[1, 1]
    signed = conc[conc["direction"] != "any"]
    c = signed[signed["conc_padj"] <= 0.05]
    grp = (c.groupby(["layer", "norm", "scheme"]).size()
           .rename("n").reset_index())
    # every group that was tested, so a normalisation with no hits shows as zero
    labels = sorted((signed["layer"] + " / " + signed["norm"]).unique().tolist())
    schemes = sorted(grp["scheme"].unique())
    width = 0.8 / len(schemes)
    for i, sc in enumerate(schemes):
        vals = [grp.loc[(grp["layer"] + " / " + grp["norm"] == lab)
                        & (grp["scheme"] == sc), "n"].sum() for lab in labels]
        ax.bar(np.arange(len(labels)) + i * width - 0.4 + width / 2, vals, width, label=sc)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("pathways concentrated across cohorts\n(FDR < 0.05)")
    ax.set_yscale("symlog", linthresh=10)
    ax.legend(frameon=False, fontsize=8.5, title="threshold scheme", title_fontsize=8.5)
    ax.set_title("D  Cross-cohort concentration by layer", loc="left",
                 fontsize=11, fontweight="bold")

    for a_ in axes.ravel():
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"norm_summary.{ext}", dpi=200, bbox_inches="tight")
    print(f"written -> {FIGURES / 'norm_summary.png'}")


if __name__ == "__main__":
    main()
