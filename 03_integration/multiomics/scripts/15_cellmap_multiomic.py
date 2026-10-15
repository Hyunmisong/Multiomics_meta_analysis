"""The closing figure: the integrated Reactome result on a schematic cell.

Script 05 draws one layer's cross-cohort concentration. This one draws the
*integrated* call - the pathways that survive all six three-layer analyses
(3 methods x 2 universes, scripts/12) - and keeps the cellular layout, so the
result can be read by location as well as by pathway:

    compartment boxes   where the pathway acts
    node card           one pathway; the badge is how many of the 6 analyses
                        called it, and how many of its measured members two or
                        more layers agree on
    evidence tracks     one cell per cohort, per layer: red enriched among
                        up-regulated genes, blue among down, grey measured but
                        not enriched, pale not measured
    grey lines          pathways sharing members (Jaccard >= 0.10)

Reactome only, to match the network figures: the diagrams next to this one are
drawn in Reactome reaction order, so the map has to name the same pathways.
Scoring itself still uses all three collections.

    python scripts/15_cellmap_multiomic.py
    python scripts/15_cellmap_multiomic.py --top-n 30 --direction down
    python scripts/15_cellmap_multiomic.py --layout results/cellmap_multiomic_layout.yaml
"""
import _bootstrap  # noqa: F401

import argparse

import pandas as pd
import yaml

from pmap import cellmap, genesets
from pmap.paths import FIGURES, RESULTS, load_config


def select(res: pd.DataFrame, collection: str, universe: str, method: str,
           direction: str, top_n: int, min_analyses: int,
           top_down: int = 10) -> pd.DataFrame:
    """Pathways called by `min_analyses` of the six analyses, best first.

    With direction "any" the up and down lists are taken separately: 323
    Reactome sets are called up against 23 down, so one ranked list would be
    all up and the figure would lose the metabolic half of the result.
    """
    sig = res[res["padj"] <= 0.05]
    agree = (sig.groupby(["collection", "pathway", "direction"]).size()
             .rename("n_analyses").reset_index())
    sel = (res[(res["universe"] == universe) & (res["method"] == method)
               & (res["padj"] <= 0.05) & (res["collection"] == collection)]
           .merge(agree, on=["collection", "pathway", "direction"]))
    if direction != "any":
        sel = sel[sel["direction"] == direction]
    sel = sel[sel["n_analyses"] >= min_analyses]
    # one row per pathway: keep the stronger direction if both are called
    sel = (sel.sort_values(["n_analyses", "padj"], ascending=[False, True])
           .drop_duplicates("pathway"))
    if direction != "any":
        return sel.head(top_n).reset_index(drop=True)
    down = sel[sel["direction"] == "down"].head(top_down)
    up = sel[sel["direction"] == "up"].head(max(top_n - len(down), 0))
    return pd.concat([up, down]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collection", default="Reactome")
    ap.add_argument("--universe", default="three_layers")
    ap.add_argument("--method", default="wilcoxon")
    ap.add_argument("--direction", default="any", choices=["any", "up", "down"])
    ap.add_argument("--top-n", type=int, default=26)
    ap.add_argument("--top-down", type=int, default=10,
                    help="how many of the top-n slots go to down pathways")
    ap.add_argument("--min-analyses", type=int, default=6)
    ap.add_argument("--max-members", type=int, default=400)
    ap.add_argument("--no-phospho", action="store_true",
                    help="drop the two phospho tracks; written as a separate file")
    ap.add_argument("--layout", help="YAML of manual {pathway: {x, y}} positions")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    cm = load_config("cellmap.yaml")
    if a.no_phospho:
        cm["tracks"] = [t for t in cm["tracks"] if t["layer"] != "phospho"]
    reg = load_config("datasets.yaml")["datasets"]
    res = pd.read_csv(RESULTS / "pathway_multiomic.csv")
    ev = pd.read_csv(RESULTS / "gene_multiomic.csv")
    ora = pd.read_parquet(RESULTS / "ora.parquet")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)

    picked = select(res, a.collection, a.universe, a.method, a.direction,
                    a.top_n, a.min_analyses, a.top_down)
    picked = picked[picked["n_members"] <= a.max_members]
    if picked.empty:
        raise SystemExit("no pathway passes the filter")
    print(f"{len(picked)} {a.collection} pathways, {a.min_analyses}/6 analyses "
          f"({a.universe}, {a.method}, {a.direction})")

    # how many of each pathway's members two or more layers agree on
    best = (ev.sort_values("z_combined", ascending=False)
            .groupby("gene").first())
    sup = best["n_layers_supporting"]
    called = {lay: best[f"{lay}_k"].fillna(0) >= 1
              for lay in ("proteome", "scrna", "genome")}
    support = []
    rules, default = cm["compartments"]["rules"], cm["compartments"]["default"]
    canvas = dict(cm["canvas"])
    nodes = []
    for r in picked.itertuples(index=False):
        members = sets.get((r.collection, r.pathway), set())
        s = sup.reindex(sorted(members)).fillna(0)
        nodes.append(cellmap.Node(
            pathway=r.pathway,
            compartment=cellmap.assign_compartment(r.pathway, rules, default),
            w=canvas["node_w"], h=canvas["node_h"],
            meta={"badge": f"{r.direction} {int(r.n_analyses)}/6",
                  "support": (int((s >= 2).sum()), int((s >= 3).sum())),
                  "tip": (f"{r.direction} | {int(r.n_members)} members, "
                          f"{int((s >= 1).sum())} changed, "
                          f"{int((s >= 2).sum())} in 2+ layers, "
                          f"{int((s >= 3).sum())} in all 3 | "
                          f"{a.method} FDR {r.padj:.2g}")}))
        support.append({"pathway": r.pathway, "direction": r.direction,
                        "n_members": int(r.n_members),
                        "n_changed": int((s >= 1).sum()),
                        "n_2_layers": int((s >= 2).sum()),
                        "n_3_layers": int((s >= 3).sum()),
                        **{f"called_{k}": int(v.reindex(sorted(members)).astype("boolean").fillna(False).sum())
                           for k, v in called.items()}})

    manual = dict(cm.get("manual_positions") or {})
    if a.layout:
        with open(a.layout, encoding="utf-8") as fh:
            manual.update(yaml.safe_load(fh) or {})

    boxes, canvas, cell_rect = cellmap.layout(nodes, canvas, manual)
    names = [n.pathway for n in nodes]
    edges = cellmap.pathway_edges(sets, names, a.collection)
    scheme_for = load_config("analysis.yaml").get("cross_layer_scheme") or {}
    grid = cellmap.evidence_grid(ora, names, a.collection, "fdr05",
                                 cm["tracks"], reg, scheme_for=scheme_for)

    svg = cellmap.render_svg(
        nodes, boxes, edges, grid, cm["tracks"], canvas, cell_rect,
        title="Integrated multi-omic pathway map (Reactome)",
        subtitle=(f"pathways called by {a.min_analyses} of 6 three-layer analyses "
                  f"(3 methods x 2 universes); badge = analyses agreeing, "
                  f"tracks = per-cohort enrichment in "
                  + ("proteome, genome and scRNA-seq" if a.no_phospho else
                     "proteome, phospho (raw and protein-adjusted), genome "
                     "and scRNA-seq")))

    out = a.out or FIGURES / (f"cellmap_{a.collection}_multiomic"
                              + ("_nophospho" if a.no_phospho else "") + ".svg")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(svg)
    with open(str(out).replace(".svg", ".html"), "w", encoding="utf-8") as fh:
        fh.write('<!doctype html><meta charset="utf-8">'
                 '<body style="margin:0;background:#fff;overflow-x:auto">' + svg + "</body>")
    lay = RESULTS / "cellmap_multiomic_layout.yaml"
    with open(lay, "w", encoding="utf-8") as fh:
        yaml.safe_dump({n.pathway: {"x": round(n.x, 1), "y": round(n.y, 1),
                                    "compartment": n.compartment} for n in nodes},
                       fh, sort_keys=True, allow_unicode=True)
    tab = pd.DataFrame(support)
    tab.to_csv(RESULTS / "pathway_layer_support.csv", index=False)
    print(f"figure -> {out}\nlayout -> {lay} (edit and pass --layout)")
    print("genes called, summed over the drawn pathways: "
          + ", ".join(f"{k.split('_')[1]} {int(tab[k].sum())}"
                      for k in ("called_proteome", "called_scrna", "called_genome"))
          + f" | members agreed by 3 layers: {int(tab['n_3_layers'].sum())}")
    for comp, group in cellmap._groupby(nodes, lambda n: n.compartment):
        print(f"  {comp:16s} {len(group):2d}  " + ", ".join(n.pathway for n in group)[:88])


if __name__ == "__main__":
    main()
