"""Which genes drive a pathway, and in which tumour type and which layer.

Two outputs:

  results/tissue_pathway_examples.csv   every pathway that appears or disappears
                                        when the cohorts are split by tumour
                                        type, with its FDR in each analysis
  figures/pathway_components.png        for a handful of those pathways, the
                                        member genes with their pooled evidence
                                        in each tumour type and each layer

The proteome columns are signed (red = higher in tumour, blue = lower); the
genome columns are not, because a truncating variant is loss of function and
has only one direction. The two therefore get separate colour scales, and a
gene can only be dark in the genome columns by being recurrently truncated.
"""
import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

from pmap.paths import FIGURES, RESULTS, ROOT

TOP_N = 1000
ALPHA = 0.05
MAX_GENES = 16          # rows per pathway panel

# One example of each kind, chosen from results/tissue_pathway_examples.csv.
PANELS = [
    ("Hallmark | Oxidative Phosphorylation", "gastric only"),
    ("Reactome | Citric Acid Cycle (TCA Cycle)", "gained by splitting (gastric)"),
    ("KEGG | RNA transport", "oesophageal only"),
    ("KEGG | Protein processing in endoplasmic reticulum", "gained by splitting (oesophageal)"),
    ("Reactome | mRNA Splicing", "both tumour types"),
    ("KEGG | Non-homologous end-joining", "lost by splitting"),
]


def load_pooled():
    out = {}
    for layer in ["proteome", "genome"]:
        for t in ["gastric", "esophageal", "all"]:
            f = RESULTS / f"pooled_{layer}_{t}_stouffer.csv"
            if f.exists():
                d = pd.read_csv(f).set_index("gene")
                out[(layer, t)] = d
    return out


def examples_table(po: pd.DataFrame) -> pd.DataFrame:
    """Every pathway classified by what the split does to it."""
    rows = []
    for layer, s in po.groupby("layer"):
        padj, sig = {}, {}
        for grp, t in s.groupby("group"):
            k = t["collection"] + " | " + t["pathway"]
            padj[grp] = dict(zip(k, t["padj"]))
            sig[grp] = set(k[t["padj"] <= ALPHA])
        g, e, a = sig.get("gastric", set()), sig.get("esophageal", set()), sig.get("all", set())
        for k in sorted(g | e | a):
            if k in g and k in e:
                cls = "both tumour types"
            elif k in g:
                cls = "gastric only"
            elif k in e:
                cls = "oesophageal only"
            else:
                cls = "pooled only"
            if k not in a and (k in g or k in e):
                cls += " (gained by splitting)"
            if k in a and k not in (g | e):
                cls = "lost by splitting"
            rows.append({
                "layer": layer, "pathway": k, "classification": cls,
                "padj_pooled": padj.get("all", {}).get(k, np.nan),
                "padj_gastric": padj.get("gastric", {}).get(k, np.nan),
                "padj_esophageal": padj.get("esophageal", {}).get(k, np.nan)})
    df = pd.DataFrame(rows)
    df["best_split_padj"] = df[["padj_gastric", "padj_esophageal"]].min(axis=1)
    return df.sort_values(["layer", "classification", "best_split_padj"])


def panel_matrix(pathway: str, members: set, pooled: dict):
    """Genes x (proteome gastric/eso, genome gastric/eso) of pooled Z."""
    cols = [("proteome", "gastric"), ("proteome", "esophageal"),
            ("genome", "gastric"), ("genome", "esophageal")]
    data = {}
    for key in cols:
        d = pooled.get(key)
        if d is None:
            data[key] = pd.Series(dtype=float)
            continue
        z = d["Z"].reindex(sorted(members))
        if key[0] == "genome":
            z = z.clip(lower=0)        # one-sided: only loss evidence counts
        data[key] = z
    M = pd.DataFrame(data)
    M.columns = ["PROT\ngastric", "PROT\noeso", "GEN\ngastric", "GEN\noeso"]
    prot = M.iloc[:, :2]
    keep = prot.abs().max(axis=1).sort_values(ascending=False).head(MAX_GENES).index
    return M.loc[keep]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    po = pd.read_csv(RESULTS / "pooled_ora.csv.gz")
    po = po[po.top_n == TOP_N]
    ex = examples_table(po)
    ex.to_csv(RESULTS / "tissue_pathway_examples.csv", index=False)
    print(f"{len(ex):,} pathway classifications -> results/tissue_pathway_examples.csv")
    print(ex.groupby(["layer", "classification"]).size().to_string())

    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    gs["key"] = gs["collection"] + " | " + gs["pathway"]
    members = {k: set(v) for k, v in gs.groupby("key")["gene"]}
    pooled = load_pooled()

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 10.5))
    fig.suptitle("Who drives each pathway: pooled evidence per gene, "
                 "by tumour type and layer", fontsize=13, fontweight="bold",
                 x=0.01, ha="left")

    pnorm = TwoSlopeNorm(vcenter=0, vmin=-12, vmax=12)
    gnorm = Normalize(vmin=0, vmax=6)
    for ax, (key, label) in zip(axes.ravel(), PANELS):
        if key not in members:
            ax.axis("off")
            continue
        M = panel_matrix(key, members[key], pooled)
        if M.empty:
            ax.axis("off")
            continue
        ax.imshow(M.iloc[:, :2].values, cmap="RdBu_r", norm=pnorm, aspect="auto",
                  extent=(-0.5, 1.5, len(M) - 0.5, -0.5))
        ax.imshow(M.iloc[:, 2:].values, cmap="Purples", norm=gnorm, aspect="auto",
                  extent=(1.5, 3.5, len(M) - 0.5, -0.5))
        ax.axvline(1.5, color="#333", lw=1.6)
        ax.set_xticks(range(4))
        ax.set_xticklabels(M.columns, fontsize=8)
        ax.set_yticks(range(len(M)))
        ax.set_yticklabels(M.index, fontsize=7.5)
        ax.set_xlim(-0.5, 3.5)
        name = key.split(" | ", 1)[1]
        ax.set_title(f"{name}\n{label}", loc="left", fontsize=9.5, fontweight="bold")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)

    fig.subplots_adjust(bottom=0.11)
    cax1 = fig.add_axes([0.13, 0.045, 0.3, 0.014])
    cax2 = fig.add_axes([0.57, 0.045, 0.3, 0.014])
    fig.colorbar(plt.cm.ScalarMappable(norm=pnorm, cmap="RdBu_r"), cax=cax1,
                 orientation="horizontal").set_label(
        "proteome pooled Z  (red = up in tumour, blue = down)", fontsize=8.5)
    fig.colorbar(plt.cm.ScalarMappable(norm=gnorm, cmap="Purples"), cax=cax2,
                 orientation="horizontal").set_label(
        "genome pooled Z  (recurrent truncation = loss)", fontsize=8.5)

    fig.tight_layout(rect=[0, 0.09, 1, 0.95])
    for e in ("png", "pdf"):
        fig.savefig(FIGURES / f"pathway_components.{e}", dpi=200, bbox_inches="tight")
    print(f"figure -> {(FIGURES / 'pathway_components.png').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
