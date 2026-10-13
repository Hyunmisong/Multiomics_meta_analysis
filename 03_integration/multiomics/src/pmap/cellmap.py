"""Lay pathways out on a schematic cell and render the multi-omic overlay.

The figure is an SVG built directly (no plotting library) so that the layout
stays editable: every node's position is dumped to results/cellmap_layout.yaml
and read back on the next run, so an automatic first pass can be nudged by hand
without rerunning the analysis.

Each pathway is drawn as a card carrying one evidence track per omics layer.
A track is a row of cells, one per dataset in that layer, filled red when the
dataset shows the pathway enriched among up-regulated genes, blue when among
down-regulated, purple when both, and left empty when the dataset measured the
pathway but found nothing. Phosphoproteomics gets two tracks, raw and
protein-adjusted, because the whole point of the prototype is that they differ.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ palette
PALETTE = {
    "bg": "#ffffff",
    "cell_fill": "#fbfbfd", "cell_stroke": "#9aa3b0",
    "nucleus": "#eef1f8", "nucleus_stroke": "#6f7da3",
    "mito": "#f3eef8", "mito_stroke": "#8d7aa8",
    "er": "#eef6f2", "er_stroke": "#6f9f88",
    "ecm": "#f8f3ee", "ecm_stroke": "#a8907a",
    "membrane": "#eef5f8", "membrane_stroke": "#7a98a8",
    "cytosol": "#fcfcfd", "cytosol_stroke": "#b7bdc7",
    "node_fill": "#ffffff", "node_stroke": "#5b6472",
    "text": "#1d2129", "muted": "#6b7280",
    "up": "#c0392b", "down": "#2a6fb0", "both": "#7d4a9e",
    "none": "#e9ecef", "untested": "#f7f8fa", "edge": "#c3c9d2",
}
COMPARTMENT_STYLE = {
    "extracellular": ("ecm", "ecm_stroke", "Extracellular / secreted"),
    "membrane": ("membrane", "membrane_stroke", "Plasma membrane & signalling"),
    "cytosol": ("cytosol", "cytosol_stroke", "Cytosol"),
    "er_golgi": ("er", "er_stroke", "ER / Golgi / vesicles"),
    "mitochondrion": ("mito", "mito_stroke", "Mitochondrion"),
    "nucleus": ("nucleus", "nucleus_stroke", "Nucleus"),
}


# --------------------------------------------------------------- assignment
def assign_compartment(pathway: str, rules: dict, default: str) -> str:
    name = pathway.lower()
    for comp, patterns in rules.items():
        for pat in patterns:
            if re.search(pat, name):
                return comp
    return default


# ---------------------------------------------------------------- selection
def select_pathways(conc: pd.DataFrame, sel: dict) -> pd.DataFrame:
    df = conc[(conc["collection"] == sel["collection"])
              & (conc["scheme"] == sel["scheme"])
              & (conc["direction"] == sel["direction"])
              & (conc["layer"] == sel["layer"])].copy()
    df = df[~df["pathway"].isin(sel.get("drop") or [])]
    df = df[df["conc_padj"] <= sel.get("max_padj", 1.0)]
    df = df.sort_values(["n_cohorts_enriched", "stouffer_z"], ascending=False)
    top = df.head(sel.get("top_n", 24))
    pinned = conc[(conc["collection"] == sel["collection"])
                  & (conc["scheme"] == sel["scheme"])
                  & (conc["direction"] == sel["direction"])
                  & (conc["layer"] == sel["layer"])
                  & (conc["pathway"].isin(sel.get("pin") or []))]
    return pd.concat([top, pinned]).drop_duplicates("pathway").reset_index(drop=True)


# ------------------------------------------------------------ evidence grid
def evidence_grid(ora: pd.DataFrame, pathways: List[str], collection: str,
                  scheme: str, tracks: List[dict], registry: dict,
                  alpha: float = 0.05,
                  scheme_for: Dict[str, str] | None = None) -> Dict[str, Dict[str, List[dict]]]:
    """{pathway: {track_key: [{dataset, state}, ...]}}.

    state is one of up / down / both / none (tested, not enriched) /
    untested (the dataset does not measure enough of the pathway).
    """
    # A layer may be scored by its own schemes (the genome layer is), so each
    # track reads the scheme that layer was actually run under.
    scheme_for = scheme_for or {}
    sub = ora[(ora["collection"] == collection) & (ora["pathway"].isin(pathways))]
    grid = {p: {} for p in pathways}
    for tr in tracks:
        arms = [e for e in registry
                if e["layer"] == tr["layer"]
                and (tr["norm"] is None or e.get("norm") == tr["norm"])]
        datasets = sorted({e["name"] for e in arms})
        t = sub[(sub["layer"] == tr["layer"])
                & (sub["scheme"] == scheme_for.get(tr["layer"], scheme))]
        if tr["norm"] is not None:
            t = t[t["norm"] == tr["norm"]]
        for p in pathways:
            tp = t[t["pathway"] == p]
            cells = []
            for ds in datasets:
                d = tp[tp["dataset"] == ds]
                if d.empty:
                    cells.append({"dataset": ds, "state": "untested"})
                    continue
                up = bool((d[(d["direction"] == "up") & (d["padj"] <= alpha)]).shape[0])
                dn = bool((d[(d["direction"] == "down") & (d["padj"] <= alpha)]).shape[0])
                state = "both" if up and dn else "up" if up else "down" if dn else "none"
                cells.append({"dataset": ds, "state": state})
            grid[p][tr["key"]] = cells
    return grid


# ------------------------------------------------------------------- layout
@dataclass
class Node:
    pathway: str
    compartment: str
    x: float = 0.0
    y: float = 0.0
    w: float = 182
    h: float = 62
    meta: dict = field(default_factory=dict)


COLS = {"extracellular": 1, "membrane": 3, "cytosol": 2,
        "er_golgi": 2, "nucleus": 2, "mitochondrion": 2}


def _box_height(n: int, cols: int, nh: float, gy: float, label: float = 26) -> float:
    rows = max(1, -(-n // cols))
    return label + gy + rows * (nh + gy)


def _compartment_boxes(canvas: dict, counts: Dict[str, int]) -> tuple:
    """Schematic cell sized to its contents.

    Layout, mirroring how these figures are usually drawn: a secreted/ECM strip
    outside the cell on the left; inside the cell a plasma-membrane band across
    the top, then two columns -- cytosol and ER/Golgi on the left, nucleus and
    mitochondrion on the right. Every box shrinks to the rows it actually needs,
    so the canvas height follows the number of pathways instead of being fixed.
    """
    M = canvas["margin"]
    nw, nh = canvas["node_w"], canvas["node_h"]
    gx, gy = canvas["node_gap_x"], canvas["node_gap_y"]
    top = M + 56

    def width(cols):
        return cols * nw + (cols + 1) * gx

    boxes = {}
    # outside the cell
    ecm_w = width(COLS["extracellular"])
    ecm_h = _box_height(counts.get("extracellular", 0), COLS["extracellular"], nh, gy)
    cell_x = M + ecm_w + 26
    # the cell must hold the two inner columns side by side
    cell_w = 2 * gx + width(COLS["cytosol"]) + 8 + width(COLS["nucleus"]) + 12

    y = top + 22
    if counts.get("membrane"):
        h = _box_height(counts["membrane"], COLS["membrane"], nh, gy)
        boxes["membrane"] = (cell_x + gx, y, cell_w - 2 * gx, h)
        y += h + gy + 6
    left_y = right_y = y
    left_x = cell_x + gx
    right_x = cell_x + gx + width(COLS["cytosol"]) + 8
    for comp, col in (("cytosol", "L"), ("er_golgi", "L"),
                      ("nucleus", "R"), ("mitochondrion", "R")):
        if not counts.get(comp):
            continue
        h = _box_height(counts[comp], COLS[comp], nh, gy)
        if col == "L":
            boxes[comp] = (left_x, left_y, width(COLS[comp]), h)
            left_y += h + gy + 6
        else:
            boxes[comp] = (right_x, right_y, width(COLS[comp]), h)
            right_y += h + gy + 6

    cell_h = max(left_y, right_y) - (top + 22) + 20
    if counts.get("extracellular"):
        boxes["extracellular"] = (M, top + 22, ecm_w, ecm_h)
    canvas = dict(canvas)
    canvas["width"] = int(cell_x + cell_w + M)
    canvas["height"] = int(top + 22 + cell_h + 54)
    cell_rect = (cell_x, top, cell_w, cell_h + 44)
    return boxes, canvas, cell_rect


def layout(nodes: List[Node], canvas: dict, manual: dict) -> tuple:
    """Grid-pack nodes inside their compartment box, then apply overrides.

    A packed grid beats a force layout here: compartments are small, the cards
    are wide, and a reader scans a figure like this row by row.
    """
    counts = {c: len(g) for c, g in _groupby(nodes, lambda n: n.compartment)}
    boxes, canvas, cell_rect = _compartment_boxes(canvas, counts)
    gx, gy = canvas["node_gap_x"], canvas["node_gap_y"]
    for comp, group in _groupby(nodes, lambda n: n.compartment):
        bx, by, bw, bh = boxes[comp]
        cols = COLS.get(comp, 2)
        for i, n in enumerate(group):
            r, c = divmod(i, cols)
            n.x = bx + gx + c * (n.w + gx)
            n.y = by + 26 + gy + r * (n.h + gy)
    for n in nodes:
        if n.pathway in manual:
            n.x, n.y = manual[n.pathway]["x"], manual[n.pathway]["y"]
    return boxes, canvas, cell_rect


def _groupby(items, key):
    out = {}
    for it in items:
        out.setdefault(key(it), []).append(it)
    return out.items()


def pathway_edges(sets: dict, pathways: List[str], collection: str,
                  min_jaccard: float = 0.10) -> List[tuple]:
    """Edges between pathways sharing genes -- the 'relationship' layer."""
    edges = []
    for i, a in enumerate(pathways):
        ga = sets.get((collection, a), set())
        for b in pathways[i + 1:]:
            gb = sets.get((collection, b), set())
            if not ga or not gb:
                continue
            j = len(ga & gb) / len(ga | gb)
            if j >= min_jaccard:
                edges.append((a, b, j))
    return edges


# ------------------------------------------------------------------- render
def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _wrap(text: str, width: int = 26, max_lines: int = 2) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 <= width:
            cur = f"{cur} {w}".strip()
        else:
            lines.append(cur)
            cur = w
            if len(lines) == max_lines:
                break
    if cur and len(lines) < max_lines:
        lines.append(cur)
    if len(lines) == max_lines and len(" ".join(lines)) < len(text):
        lines[-1] = lines[-1][:width - 1] + "…"
    return lines


def render_svg(nodes: List[Node], boxes: dict, edges: List[tuple],
               grid: dict, tracks: List[dict], canvas: dict, cell_rect: tuple,
               title: str, subtitle: str) -> str:
    W, H, M = canvas["width"], canvas["height"], canvas["margin"]
    pos = {n.pathway: n for n in nodes}
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}" font-family="Inter, Helvetica, Arial, sans-serif">',
           f'<rect width="{W}" height="{H}" fill="{PALETTE["bg"]}"/>',
           f'<text x="{M}" y="32" font-size="21" font-weight="700" '
           f'fill="{PALETTE["text"]}">{_esc(title)}</text>',
           f'<text x="{M}" y="52" font-size="11.5" '
           f'fill="{PALETTE["muted"]}">{_esc(subtitle)}</text>']

    cx, cy, cw, ch = cell_rect
    out.append(f'<rect x="{cx:.1f}" y="{cy:.1f}" width="{cw:.1f}" height="{ch:.1f}" '
               f'rx="40" fill="{PALETTE["cell_fill"]}" stroke="{PALETTE["cell_stroke"]}" '
               f'stroke-width="2.5"/>')
    out.append(f'<text x="{cx + 16:.1f}" y="{cy + 17:.1f}" font-size="10.5" '
               f'font-weight="700" letter-spacing="0.8" fill="{PALETTE["cell_stroke"]}">'
               f'TUMOUR CELL</text>')

    for comp, (x, y, w, h) in boxes.items():
        fill, stroke, label = COMPARTMENT_STYLE[comp]
        rx = 34 if comp in ("nucleus", "mitochondrion") else 12
        dash = ' stroke-dasharray="6 4"' if comp == "cytosol" else ""
        out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
                   f'rx="{rx}" fill="{PALETTE[fill]}" stroke="{PALETTE[stroke]}" '
                   f'stroke-width="1.6"{dash}/>')
        out.append(f'<text x="{x + 14:.1f}" y="{y + 18:.1f}" font-size="10.5" '
                   f'font-weight="700" letter-spacing="0.7" '
                   f'fill="{PALETTE[stroke]}">{_esc(label.upper())}</text>')

    for a, b, j in edges:
        na, nb = pos[a], pos[b]
        out.append(f'<line x1="{na.x + na.w / 2:.1f}" y1="{na.y + na.h / 2:.1f}" '
                   f'x2="{nb.x + nb.w / 2:.1f}" y2="{nb.y + nb.h / 2:.1f}" '
                   f'stroke="{PALETTE["edge"]}" stroke-width="{0.6 + 4 * j:.2f}" '
                   f'opacity="0.5"/>')

    for n in nodes:
        out.append(_node_svg(n, grid.get(n.pathway, {}), tracks))

    out.append(_legend_svg(W, H, M))
    out.append("</svg>")
    return "\n".join(out)


def _node_svg(n: Node, ev: dict, tracks: List[dict]) -> str:
    x, y, w, h = n.x, n.y, n.w, n.h
    k = n.meta.get("n_cohorts_enriched", 0)
    s = [f'<g><rect x="{x:.1f}" y="{y:.1f}" width="{w}" height="{h}" rx="7" '
         f'fill="{PALETTE["node_fill"]}" stroke="{PALETTE["node_stroke"]}" '
         f'stroke-width="1.1"/>',
         f'<title>{_esc(n.pathway)} — enriched in {k} cohorts '
         f'(concentration FDR {n.meta.get("conc_padj", float("nan")):.3g})</title>']
    lines = _wrap(n.pathway, 27, 2)
    for i, line in enumerate(lines):
        s.append(f'<text x="{x + 9:.1f}" y="{y + 15 + i * 11:.1f}" font-size="9.8" '
                 f'font-weight="600" fill="{PALETTE["text"]}">{_esc(line)}</text>')
    s.append(f'<text x="{x + w - 9:.1f}" y="{y + 15:.1f}" font-size="9" '
             f'text-anchor="end" fill="{PALETTE["muted"]}">{k}/8</text>')

    ty = y + 24 + len(lines) * 11 - 11
    for tr in tracks:
        cells = ev.get(tr["key"], [])
        s.append(f'<text x="{x + 9:.1f}" y="{ty + 6.6:.1f}" font-size="6.8" '
                 f'fill="{PALETTE["muted"]}">{_esc(tr["label"])}</text>')
        cxx = x + 46
        if not cells:
            s.append(f'<text x="{cxx:.1f}" y="{ty + 6.6:.1f}" font-size="6.4" '
                     f'fill="#c9ced6">no data yet</text>')
        for c in cells:
            col = {"up": PALETTE["up"], "down": PALETTE["down"], "both": PALETTE["both"],
                   "none": PALETTE["none"], "untested": PALETTE["untested"]}[c["state"]]
            s.append(f'<rect x="{cxx:.1f}" y="{ty:.1f}" width="8" height="7" rx="1.3" '
                     f'fill="{col}" stroke="#d2d6dd" stroke-width="0.4">'
                     f'<title>{_esc(tr["label"])} {_esc(c["dataset"])}: {c["state"]}</title></rect>')
            cxx += 9.4
        ty += 8.6
    s.append("</g>")
    return "".join(s)


def _legend_svg(W: int, H: int, M: int) -> str:
    x, y = M, H - 18
    items = [("up", "enriched among up"), ("down", "enriched among down"),
             ("both", "both directions"), ("none", "measured, not enriched"),
             ("untested", "not measured")]
    s = []
    for key, label in items:
        s.append(f'<rect x="{x}" y="{y - 8}" width="9" height="9" rx="1.5" '
                 f'fill="{PALETTE[key]}" stroke="#d2d6dd" stroke-width="0.5"/>')
        s.append(f'<text x="{x + 13}" y="{y}" font-size="9.5" '
                 f'fill="{PALETTE["muted"]}">{label}</text>')
        x += 22 + len(label) * 5.4
    s.append(f'<text x="{W - M}" y="{y}" font-size="9.5" '
             f'text-anchor="end" fill="{PALETTE["muted"]}">'
             f'one cell = one cohort; tracks = PROT / phospho raw / phospho protein-adjusted '
             f'/ genomics / scRNA-seq</text>')
    return "".join(s)
