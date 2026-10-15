"""One diagram per Reactome pathway, drawn in the order the pathway runs.

Scoring uses all three gene-set collections (see 12_multiomic.py); only
Reactome is drawn, because only Reactome states what happens before what.
Every reaction lists its participants and the reactions that must precede it,
so the figure can put the sequence on the page:

    grey squares   reactions
    blue arrows    precedence - the running order, left to right
    grey arrows    gene -> reaction (input or catalyst), reaction -> gene (output)
    circles        member genes, coloured by direction and by how many omics
                   layers agree

Colour is a two-ended gradient, six steps in all:

    up    1 layer pale  ->  2 layers mid red   ->  3 layers deep red
    down  1 layer pale  ->  2 layers mid blue  ->  3 layers deep blue
    grey  measured, called by no layer

Up can reach at most two layers: a truncating variant is loss of function, so
the genome layer can only ever support "down".

A Reactome set with no reaction chain - an assembly step such as Complex I
Biogenesis is one reaction with 68 members - has no order to show and is
skipped rather than drawn some other way.

    python scripts/13_pathway_network.py --top 30
    python scripts/13_pathway_network.py --top 0      # every pathway that resolves
"""
import _bootstrap  # noqa: F401

import argparse
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

from pmap import genesets, reactome
from pmap.paths import FIGURES, RESULTS, ROOT

UP = ["#fcbba1", "#fb6a4a", "#a50f15"]       # 1, 2, 3 layers
DOWN = ["#c6dbef", "#4292c6", "#08306b"]
NONE = "#9e9e9e"
UP_CMAP = LinearSegmentedColormap.from_list("up", ["#fff5f0"] + UP)
DOWN_CMAP = LinearSegmentedColormap.from_list("down", ["#f7fbff"] + DOWN)
OUT = FIGURES / "pathway_networks"


def node_color(n, d):
    if n < 1:
        return NONE
    return (UP if d == "up" else DOWN)[min(int(n), 3) - 1]


def gene_state(ev: pd.DataFrame, genes):
    """Each gene's stronger direction, its layer support and its combined z."""
    e = ev[ev["gene"].isin(genes)]
    best = e.sort_values("z_combined", ascending=False).groupby("gene").first()
    idx = sorted(genes)
    return (best["n_layers_supporting"].reindex(idx).fillna(0).astype(int),
            best["direction"].reindex(idx).fillna(""),
            best["z_combined"].reindex(idx).fillna(0.0))


def reaction_layout(G, layer):
    """Columns in running order: reaction layer on x, members spread on y."""
    cols = {}
    for n in G.nodes:
        cols.setdefault(layer.get(n, 0), []).append(n)
    ncol = max(len(cols), 1)
    pos = {}
    for i, (lay, nodes) in enumerate(sorted(cols.items())):
        rx = sorted(n for n in nodes if G.nodes[n].get("kind") == "reaction")
        gn = sorted(n for n in nodes if G.nodes[n].get("kind") != "reaction")
        x = i / max(ncol - 1, 1)
        for j, n in enumerate(rx):
            pos[n] = (x, 0.93 - 0.86 * (j + 0.5) / max(len(rx), 1))
        for j, n in enumerate(gn):
            # genes sit just right of their column's reactions, in two files so
            # a 40-gene column does not become one unreadable line
            pos[n] = (x + (0.028 if j % 2 else 0.055),
                      0.93 - 0.86 * (j + 0.5) / max(len(gn), 1))
    return pos


def draw_reactions(ax, G, layer, n_sup, direct, z, title, subtitle, max_labels=26):
    pos = reaction_layout(G, layer)
    kinds = nx.get_edge_attributes(G, "kind")
    prec = [e for e, k in kinds.items() if k == "precedes"]
    io_ = [e for e, k in kinds.items() if k in ("input", "output")]
    # membership lines are not drawn: a gene already sits in its reaction's
    # column, and several hundred of them turn the order into a haze

    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=io_, edge_color="#9aa3b0",
                           width=0.7, alpha=0.7, arrows=True, arrowsize=6,
                           node_size=60)
    # the order itself
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=prec, edge_color="#2b6cb0",
                           width=1.5, alpha=0.85, arrows=True, arrowsize=12,
                           node_size=120, connectionstyle="arc3,rad=0.12")

    rx = [n for n, k in G.nodes(data="kind") if k == "reaction"]
    gn = [n for n, k in G.nodes(data="kind") if k != "reaction"]
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=rx, node_shape="s",
                           node_color="#4a5568", node_size=55,
                           edgecolors="white", linewidths=0.5)
    zz = np.abs(np.array([z.get(g, 0.0) for g in gn])) if gn else np.array([0.0])
    nx.draw_networkx_nodes(
        G, pos, ax=ax, nodelist=gn,
        node_color=[node_color(n_sup.get(g, 0), direct.get(g, "")) for g in gn],
        node_size=50 + 170 * (zz / max(zz.max(), 1e-9)),
        edgecolors="#3f3f3f", linewidths=0.5)

    cand = sorted([g for g in gn if n_sup.get(g, 0) >= 2],
                  key=lambda g: -abs(z.get(g, 0)))[:max_labels]
    nx.draw_networkx_labels(G, pos, labels={g: g for g in cand}, ax=ax, font_size=7)
    ax.set_title(f"{title}\n", fontsize=12.5, fontweight="bold", loc="left")
    ax.text(0, 1.005, subtitle, transform=ax.transAxes, fontsize=8.5,
            color="#6b7280", va="bottom")
    ax.text(0.5, 0.055, "left to right is the running order; blue arrows are "
            "Reactome precedence, grey squares are reactions",
            transform=ax.transAxes, fontsize=7.5, color="#9aa3b0", ha="center")
    ax.margins(0.07)
    ax.axis("off")


def legend(fig):
    """Two gradient bars: the more layers agree, the deeper the colour."""
    for cmap, label, x in [(DOWN_CMAP, "down  (layers agreeing)", 0.11),
                           (UP_CMAP, "up  (layers agreeing)", 0.55)]:
        cax = fig.add_axes([x, 0.035, 0.30, 0.016])
        cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 3), cmap=cmap),
                          cax=cax, orientation="horizontal")
        cb.set_ticks([0, 1, 2, 3])
        cb.set_ticklabels(["none", "1", "2", "3"])
        cb.set_label(label, fontsize=8.5)
        cb.outline.set_visible(False)
    fig.text(0.5, 0.082, "grey = no layer calls it   |   node size = combined z",
             fontsize=8, color="#6b7280", ha="center")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=30,
                    help="0 draws every pathway that passes the filter")
    ap.add_argument("--universe", default="three_layers")
    ap.add_argument("--method", default="wilcoxon")
    ap.add_argument("--min-analyses", type=int, default=6)
    ap.add_argument("--max-members", type=int, default=90)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    ev = pd.read_csv(RESULTS / "gene_multiomic.csv")
    res = pd.read_csv(RESULTS / "pathway_multiomic.csv")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)
    with open(RESULTS / "universes.json") as fh:
        universes = {k: set(v) for k, v in json.load(fh).items()}
    from pmap.multiomic import universes_by_coverage
    uni = universes_by_coverage(universes, 3 if a.universe == "three_layers" else 2)

    sig = res[res["padj"] <= 0.05]
    agree = (sig.groupby(["collection", "pathway", "direction"]).size()
             .rename("n_analyses").reset_index())
    sel = (res[(res["universe"] == a.universe) & (res["method"] == a.method)
               & (res["padj"] <= 0.05)]
           .merge(agree, on=["collection", "pathway", "direction"]))
    sel = sel[(sel["collection"] == "Reactome")
              & (sel["n_analyses"] >= a.min_analyses)
              & (sel["n_members"] <= a.max_members)]
    sel = sel.sort_values(["n_analyses", "padj"], ascending=[False, True])
    print(f"{len(sel)} Reactome pathways pass {a.min_analyses}/6 analyses with "
          f"<= {a.max_members} members ({a.universe}, {a.method})")

    made, skipped = [], []
    for r in sel.itertuples(index=False):
        if a.top and len(made) >= a.top:
            break
        members = sorted(sets[(r.collection, r.pathway)] & uni)
        RG, meta = reactome.pathway_graph(r.pathway, members)
        if RG is None or not reactome.usable(meta):
            skipped.append({"pathway": r.pathway,
                            "reason": meta.get("reason") or "no reaction chain",
                            "n_reactions": meta.get("n_reactions", 0),
                            "coverage": round(meta.get("coverage", 0), 3)})
            continue
        n_sup, direct, z = gene_state(ev, members)
        lay = reactome.layers(RG)
        sub = (f"{r.direction} | {r.n_members} members, "
               f"{int((n_sup >= 1).sum())} changed, {int((n_sup >= 2).sum())} in 2+ layers, "
               f"{int((n_sup >= 3).sum())} in all 3 | {r.method} FDR {r.padj:.1e}, "
               f"{r.n_analyses}/6 analyses | {meta['n_reactions']} reactions, "
               f"{meta['n_precedence']} precedence links, "
               f"{meta['coverage']:.0%} of members placed")
        fig, ax = plt.subplots(figsize=(12, 9))
        draw_reactions(ax, RG, lay, n_sup.to_dict(), direct.to_dict(), z.to_dict(),
                       r.pathway, sub)
        fig.subplots_adjust(bottom=0.11)
        legend(fig)
        fn = re.sub(r"[^\w]+", "_", f"{r.collection}_{r.pathway}")[:80] + ".png"
        fig.savefig(OUT / fn, dpi=170, bbox_inches="tight")
        plt.close(fig)
        made.append({"collection": r.collection, "pathway": r.pathway,
                     "direction": r.direction, "n_members": r.n_members,
                     "n_changed": int((n_sup >= 1).sum()),
                     "n_2plus": int((n_sup >= 2).sum()),
                     "n_3": int((n_sup >= 3).sum()),
                     "n_reactions": meta["n_reactions"],
                     "n_precedence": meta["n_precedence"],
                     "coverage": round(meta["coverage"], 3),
                     "n_analyses": r.n_analyses, "padj": r.padj, "file": fn})

    pd.DataFrame(skipped).to_csv(RESULTS / "pathway_networks_skipped.csv", index=False)
    print(f"drawn {len(made)}, skipped {len(skipped)} with no reaction chain")

    m = pd.DataFrame(made)
    m.to_csv(RESULTS / "pathway_networks.csv", index=False)
    print(m[["pathway", "direction", "n_members", "n_2plus", "n_3",
             "n_reactions", "n_precedence", "coverage"]].head(20).to_string(index=False))

    hero = m.sort_values(["n_3", "n_2plus", "n_changed"], ascending=False).iloc[0]
    import shutil
    shutil.copy(OUT / hero["file"], FIGURES / "pathway_network_top.png")
    print(f"\nREADME figure -> figures/pathway_network_top.png  ({hero['pathway']})")


if __name__ == "__main__":
    main()
