"""One pathway, gene by gene: which cohort of which layer moved it, and how.

The cell map carries one row per layer and one cell per cohort, for a whole
pathway. This opens a single pathway up: one row per member gene, one column
per cohort, so "gene A changed in 6 of 8 proteomes, 3 of 7 scRNA cohorts and
no genome" is read off the row.

    red / blue cell   that cohort calls the gene up / down (shade = |log2FC|)
    pale grey cell    the cohort measured the gene and did not call it
    dashed outline    the cohort did not measure it
    right-hand text   k/n per layer, the row's summary

Two figures are written, the same genes in two orders:

    _by_support   layers agreeing, then combined z - the strongest rows first
    _by_reaction  grouped by Reactome reaction order, so the biochemistry is
                  kept; genes that no reaction places sit in a last block

Phospho is left out by default (--with-phospho adds its arms), because only 3
of the 8 proteome cohorts carry it and the columns would be mostly empty.

    python scripts/16_pathway_gene_matrix.py
    python scripts/16_pathway_gene_matrix.py --pathway "Cell Cycle Checkpoints"
    python scripts/16_pathway_gene_matrix.py --min-layers 1 --max-genes 120
"""
import _bootstrap  # noqa: F401

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle

from pmap import genesets, reactome
from pmap.paths import FIGURES, RESULTS, load_config

UP, DOWN = "#c0392b", "#2a6fb0"
NONE, GRID = "#e9ecef", "#cfd4db"
LAYERS = ["proteome", "scrna", "genome"]
LABEL = {"proteome": "proteome", "scrna": "scRNA-seq", "genome": "genome",
         "phospho": "phospho"}


def shade(base: str, w: float) -> tuple:
    """Blend towards white: a weak call must still be readable."""
    c = np.array(matplotlib.colors.to_rgb(base))
    return tuple(1 - (1 - c) * (0.35 + 0.65 * min(max(w, 0.0), 1.0)))


def cohort_columns(reg, layers, scheme_for, default_scheme):
    """[(layer, cohort, scheme)], in layer order, cohorts sorted by name."""
    cols = []
    for lay in layers:
        names = sorted({e["name"] for e in reg if e["layer"] == lay})
        for n in names:
            cols.append((lay, n, scheme_for.get(lay, default_scheme)))
    return cols


def build(calls, feats, members, cols):
    """state[gene][col] = ('up'|'down'|'none'|'untested', |log2FC|)."""
    state = {g: {} for g in members}
    for lay, ds, sch in cols:
        c = calls[(calls["layer"] == lay) & (calls["dataset"] == ds)
                  & (calls["scheme"] == sch)]
        c = c[c["gene"].isin(members)]
        called = dict(zip(c["gene"], zip(c["direction"], c["mean_logFC"].abs())))
        measured = set(feats[(feats["layer"] == lay)
                             & (feats["dataset"] == ds)]["gene"]) & members
        for g in members:
            if g in called:
                d, fc = called[g]
                state[g][(lay, ds)] = (d, 0.0 if pd.isna(fc) else float(fc))
            else:
                state[g][(lay, ds)] = ("none" if g in measured else "untested", 0.0)
    return state


def draw(genes, cols, state, counts, title, subtitle, sections, out, fc_ref):
    nrow, ncol = len(genes), len(cols)
    gap = {lay: i for i, lay in enumerate(dict.fromkeys(c[0] for c in cols))}
    xs = [j + 1.1 * gap[lay] for j, (lay, _, _) in enumerate(cols)]
    width = 6.4 + 0.26 * (max(xs) + 1)
    fig, ax = plt.subplots(figsize=(width, 1.9 + 0.17 * nrow))

    for i, g in enumerate(genes):
        y = nrow - 1 - i
        for x, (lay, ds, _) in zip(xs, cols):
            st, fc = state[g][(lay, ds)]
            if st == "untested":
                ax.add_patch(Rectangle((x, y), 0.86, 0.86, facecolor="white",
                                       edgecolor="#9aa3b0", lw=0.45, ls=(0, (1, 1))))
            else:
                col = (NONE if st == "none"
                       else shade(UP if st == "up" else DOWN, fc / fc_ref))
                ax.add_patch(Rectangle((x, y), 0.86, 0.86, facecolor=col,
                                       edgecolor=GRID, lw=0.3))
        txt = " ".join(f"{lay[0].upper()} {counts[g][lay][0]}/{counts[g][lay][1]:<4}"
                       for lay in LAYERS if counts[g][lay][1])
        ax.text(max(xs) + 1.6, y + 0.43, txt, va="center", fontsize=6.4,
                color="#4a5568", family="monospace")
        ax.text(-0.5, y + 0.43, g, va="center", ha="right", fontsize=6.6)

    # layer headers
    head = nrow + 5.2
    for lay in dict.fromkeys(c[0] for c in cols):
        xx = [x for x, c in zip(xs, cols) if c[0] == lay]
        ax.text(np.mean(xx) + 0.43, head + 0.5, LABEL[lay], ha="center",
                fontsize=8, fontweight="600", color="#374151")
        ax.plot([min(xx), max(xx) + 0.86], [head + 0.25] * 2, color="#9aa3b0", lw=0.8)
    for x, (lay, ds, _) in zip(xs, cols):
        ax.text(x + 0.43, nrow + 0.25, ds, ha="center", va="bottom", fontsize=5.6,
                rotation=90, color="#6b7280")

    # section rules, when the rows are grouped by reaction order
    for label, start, stop in sections:
        y = nrow - stop
        ax.plot([-0.2, max(xs) + 9.1], [y - 0.12] * 2, color="#c3c9d2", lw=0.6)
        # on the right, clear of the gene names
        ax.text(max(xs) + 9.4, nrow - start - 0.1, label, ha="left", va="top",
                fontsize=6.4, style="italic", color="#8b93a1")

    ax.set_xlim(-6.2, max(xs) + 14.5)
    ax.set_ylim(-0.6, nrow + 6.4)
    ax.axis("off")
    ax.set_title(title + "\n", fontsize=12, fontweight="bold", loc="left")
    ax.text(0, 1.0, subtitle, transform=ax.transAxes, fontsize=7.6,
            color="#6b7280", va="bottom")
    ax.legend(handles=[Patch(facecolor=UP, label="called up"),
                       Patch(facecolor=DOWN, label="called down"),
                       Patch(facecolor=NONE, label="measured, not called"),
                       Patch(facecolor="white", edgecolor="#9aa3b0", ls=(0, (1, 1)),
                             label="not measured"),
                       Patch(facecolor=shade(UP, 0.2),
                             label="pale to deep = that cohort's own |log2FC|")],
              loc="lower left", bbox_to_anchor=(0, -0.055 - 6.0 / max(nrow, 1)),
              ncol=5, frameon=False, fontsize=7)
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"figure -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pathway", default=None,
                    help="default: the drawn pathway with most 2-layer genes")
    ap.add_argument("--collection", default="Reactome")
    ap.add_argument("--min-layers", type=int, default=2,
                    help="keep genes supported by at least this many layers")
    ap.add_argument("--max-genes", type=int, default=0, help="0 = no cap")
    ap.add_argument("--with-phospho", action="store_true")
    ap.add_argument("--max-reactions", type=int, default=500,
                    help="a 358-member pathway runs past the drawing default")
    a = ap.parse_args()

    reg = load_config("datasets.yaml")["datasets"]
    scheme_for = load_config("analysis.yaml").get("cross_layer_scheme") or {}
    layers = LAYERS + (["phospho"] if a.with_phospho else [])

    ev = pd.read_csv(RESULTS / "gene_multiomic.csv")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)
    pathway = a.pathway
    if not pathway:
        sup = pd.read_csv(RESULTS / "pathway_layer_support.csv")
        pathway = sup.sort_values("n_2_layers", ascending=False).iloc[0]["pathway"]
    members = sets[(a.collection, pathway)]

    best = (ev[ev["gene"].isin(members)].sort_values("z_combined", ascending=False)
            .groupby("gene").first())
    keep = best[best["n_layers_supporting"] >= a.min_layers]
    keep = keep.sort_values(["n_layers_supporting", "z_combined"], ascending=False)
    if a.max_genes:
        keep = keep.head(a.max_genes)
    genes = list(keep.index)
    if not genes:
        raise SystemExit(f"no gene in {pathway} reaches {a.min_layers} layers")

    cols = cohort_columns(reg, layers, scheme_for, "fdr05")
    calls = pd.read_parquet(RESULTS / "gene_calls.parquet",
                            columns=["layer", "dataset", "scheme", "gene",
                                     "direction", "mean_logFC"])
    feats = pd.read_parquet(RESULTS / "features_long.parquet",
                            columns=["layer", "dataset", "gene"]).drop_duplicates()
    state = build(calls, feats, set(genes), cols)

    counts, rows = {}, []
    for g in genes:
        counts[g] = {}
        for lay in layers:
            cs = [state[g][(l, d)][0] for l, d, _ in cols if l == lay]
            counts[g][lay] = (sum(c in ("up", "down") for c in cs), len(cs))
        rows.append({"pathway": pathway, "gene": g,
                     "n_layers_supporting": int(best.loc[g, "n_layers_supporting"]),
                     "direction": best.loc[g, "direction"],
                     "z_combined": round(float(best.loc[g, "z_combined"]), 3),
                     **{f"{lay}_k": counts[g][lay][0] for lay in layers},
                     **{f"{lay}_n": counts[g][lay][1] for lay in layers},
                     **{f"{l}|{d}": state[g][(l, d)][0] for l, d, _ in cols}})
    tag = re.sub(r"[^\w]+", "_", pathway)[:60]
    pd.DataFrame(rows).to_csv(RESULTS / f"gene_matrix_{tag}.csv", index=False)

    fcs = [fc for g in genes for (st, fc) in state[g].values() if st in ("up", "down")]
    fc_ref = float(np.percentile(fcs, 90)) if fcs else 1.0
    sub = (f"{len(genes)} of {len(members)} members, supported by "
           f"{a.min_layers}+ layers | one column per cohort, "
           f"{', '.join(f'{LABEL[l]} {sum(1 for c in cols if c[0] == l)}' for l in layers)}"
           f" | shade = |log2FC| (90th pct = {fc_ref:.2f})")

    draw(genes, cols, state, counts, pathway, sub, [],
         FIGURES / f"gene_matrix_{tag}_by_support.png", fc_ref)

    # --- the same genes, in the order the pathway runs ----------------------
    RG, meta = reactome.pathway_graph(pathway, sorted(members),
                                  max_reactions=a.max_reactions)
    if RG is None:
        print("no Reactome graph; the reaction-ordered figure is skipped")
        return
    lay_of = reactome.layers(RG)
    placed = {g: lay_of[g] for g in genes if g in lay_of and RG.nodes[g].get("kind") == "gene"}
    order = sorted(genes, key=lambda g: (placed.get(g, 10 ** 6),
                                         -float(best.loc[g, "z_combined"])))
    sections = []
    cur, start = None, 0
    for i, g in enumerate(order):
        k = placed.get(g)
        if k != cur:
            if cur is not None:
                sections.append((f"step {cur + 1}" if cur is not None else "", start, i))
            cur, start = k, i
    sections.append(("unplaced" if cur is None else f"step {cur + 1}", start, len(order)))
    sections = [(lbl if lbl != "step 1000001" else "no reaction", s, e)
                for lbl, s, e in sections]
    draw(order, cols, state, counts, pathway,
         sub + f" | grouped by Reactome reaction order ({meta['n_reactions']} reactions)",
         sections, FIGURES / f"gene_matrix_{tag}_by_reaction.png", fc_ref)
    print(f"table  -> {RESULTS / f'gene_matrix_{tag}.csv'}")


if __name__ == "__main__":
    main()
