"""Draw the cellular map with the multi-omic evidence overlay.

Automatic layout on the first run; coordinates are written to
results/cellmap_layout.yaml, which you can edit by hand and feed back with
--layout to nudge nodes without rerunning anything upstream.

    python scripts/05_cellmap.py
    python scripts/05_cellmap.py --collection KEGG --top-n 30
    python scripts/05_cellmap.py --layout results/cellmap_layout.yaml
"""
import _bootstrap  # noqa: F401

import argparse

import pandas as pd
import yaml

from pmap import cellmap, genesets
from pmap.paths import FIGURES, RESULTS, load_config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collection")
    ap.add_argument("--scheme")
    ap.add_argument("--direction")
    ap.add_argument("--top-n", type=int)
    ap.add_argument("--layout", help="YAML of manual {pathway: {x, y}} positions")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cm = load_config("cellmap.yaml")
    reg = load_config("datasets.yaml")["datasets"]
    sel = dict(cm["selection"])
    for k, v in [("collection", args.collection), ("scheme", args.scheme),
                 ("direction", args.direction), ("top_n", args.top_n)]:
        if v is not None:
            sel[k] = v

    conc = pd.read_csv(RESULTS / "concentration.csv")
    ora = pd.read_parquet(RESULTS / "ora.parquet")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)

    picked = cellmap.select_pathways(conc, sel)
    if picked.empty:
        raise SystemExit(f"no pathways pass the selection filter: {sel}")
    print(f"{len(picked)} pathways selected from {sel['collection']} "
          f"({sel['scheme']}, {sel['direction']}, {sel['layer']})")

    rules, default = cm["compartments"]["rules"], cm["compartments"]["default"]
    canvas = cm["canvas"]
    nodes = [cellmap.Node(pathway=r.pathway,
                          compartment=cellmap.assign_compartment(r.pathway, rules, default),
                          w=canvas["node_w"], h=canvas["node_h"],
                          meta={"n_cohorts_enriched": int(r.n_cohorts_enriched),
                                "conc_padj": float(r.conc_padj)})
             for r in picked.itertuples()]

    manual = dict(cm.get("manual_positions") or {})
    if args.layout:
        with open(args.layout, encoding="utf-8") as fh:
            manual.update(yaml.safe_load(fh) or {})

    boxes, canvas, cell_rect = cellmap.layout(nodes, canvas, manual)
    names = [n.pathway for n in nodes]
    edges = cellmap.pathway_edges(sets, names, sel["collection"])
    scheme_for = load_config("analysis.yaml").get("cross_layer_scheme") or {}
    grid = cellmap.evidence_grid(ora, names, sel["collection"], sel["scheme"],
                                 cm["tracks"], reg, scheme_for=scheme_for)

    layers = {t["layer"] for t in cm["tracks"]}
    present = {e["layer"] for e in reg}
    subtitle = (f"{sel['collection']} pathways concentrated across cohorts "
                f"({sel['scheme']}); layers present: {', '.join(sorted(layers & present))}"
                f"; awaiting: {', '.join(sorted(layers - present)) or 'none'}")
    svg = cellmap.render_svg(nodes, boxes, edges, grid, cm["tracks"], canvas,
                             cell_rect,
                             title="Multi-omic pathway concentration map",
                             subtitle=subtitle)

    out = args.out or FIGURES / f"cellmap_{sel['collection']}_{sel['scheme']}.svg"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(svg)

    with open(RESULTS / "cellmap_layout.yaml", "w", encoding="utf-8") as fh:
        yaml.safe_dump({n.pathway: {"x": round(n.x, 1), "y": round(n.y, 1),
                                    "compartment": n.compartment} for n in nodes},
                       fh, sort_keys=True, allow_unicode=True)

    # a scaled HTML wrapper, so the SVG can be eyeballed at a glance in a browser
    html_out = str(out).replace(".svg", ".html")
    with open(html_out, "w", encoding="utf-8") as fh:
        fh.write('<!doctype html><meta charset="utf-8">'
                 '<body style="margin:0;background:#fff;overflow-x:auto">' + svg + "</body>")

    print(f"figure  -> {out}")
    print(f"preview -> {html_out}")
    print(f"layout  -> {RESULTS / 'cellmap_layout.yaml'} (edit and pass --layout)")
    for comp, group in cellmap._groupby(nodes, lambda n: n.compartment):
        print(f"  {comp:16s} {len(group):2d}  " + ", ".join(n.pathway for n in group)[:90])


if __name__ == "__main__":
    main()
