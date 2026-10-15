"""Reaction order inside a Reactome pathway, as a directed graph.

STRING says which proteins bind; it cannot say what happens before what.
Reactome can: every reaction lists its participants and the reactions that must
precede it, so a pathway can be drawn in the order it runs.

The graph is bipartite on purpose. Reactions are nodes too, because the order
lives on them: `precedingEvent` is a reaction-to-reaction relation, and
collapsing it onto genes loses it (a binding reaction has the same genes on
both sides, so a gene-only graph of this pathway has edges but no direction to
speak of).

    gene -> reaction      the gene is an input or a catalyst
    reaction -> gene      the gene is in the output
    reaction -> reaction  the first must happen before the second
    gene -- reaction      participation, when the service does not resolve the
                          gene to one side (complexes are flattened by the
                          participants endpoint, which does not say which side)

Not every pathway suits this. An assembly pathway such as Complex I Biogenesis
is one reaction with 68 members, so it has no order to show; `usable()` says
when to fall back to an interaction network instead.

Everything is cached under data/genesets/reactome/: a pathway costs two calls
per reaction and the figures are redrawn often.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from typing import Dict, Iterable, Set, Tuple

import networkx as nx

from .paths import GENESETS

API = "https://reactome.org/ContentService"
UA = {"User-Agent": "Mozilla/5.0 (multiomics meta-analysis pipeline)"}
CACHE = GENESETS / "reactome"
COMPARTMENT = re.compile(r"\s*\[[^\]]*\]\s*$")
TOKEN = re.compile(r"[A-Z0-9][A-Z0-9\-]{1,14}")
TAIL = re.compile(r"\s(\S+)$")


def _get(path: str, retries: int = 3):
    CACHE.mkdir(parents=True, exist_ok=True)
    key = CACHE / (re.sub(r"[^\w.-]+", "_", path)[:120] + ".json")
    if key.exists():
        return json.loads(key.read_text(encoding="utf-8"))
    for i in range(retries):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(API + path, headers=UA), timeout=120) as r:
                d = json.loads(r.read())
            key.write_text(json.dumps(d), encoding="utf-8")
            return d
        except Exception:
            if i == retries - 1:
                key.write_text("null", encoding="utf-8")
                return None
            time.sleep(1.5 * (i + 1))


def pathway_id(name: str) -> str | None:
    """Reactome stable id for a pathway name as the gene-set file spells it.

    Only an exact name match is accepted: the search endpoint is fuzzy enough
    that anything looser would quietly draw a different pathway.
    """
    q = urllib.parse.quote(name)
    d = _get(f"/search/query?query={q}&species=Homo%20sapiens&types=Pathway&cluster=true")
    if not d:
        return None
    for res in d.get("results", []):
        for e in res.get("entries", []):
            plain = re.sub(r"<[^>]+>", "", e.get("name", "")).strip().lower()
            if plain == name.strip().lower():
                return e.get("stId")
    return None


def _from_display(display: str, allowed: Set[str]) -> Set[str]:
    """Symbols written into a physical entity's own name.

    'MCM4:MCM6:MCM7 heterotrimer [nucleoplasm]' -> {MCM4, MCM6, MCM7}. Tokens
    survive only if the pathway's gene set contains them, which discards words
    like 'heterotrimer' without needing a dictionary.
    """
    s = COMPARTMENT.sub("", display or "").upper()
    return {t for t in TOKEN.findall(s) if t in allowed}


def _participants(stid: str, allowed: Set[str]) -> Set[str]:
    """Every gene in a reaction, with complexes and sets flattened."""
    out = set()
    for p in _get(f"/data/participants/{stid}") or []:
        if not isinstance(p, dict):
            continue
        for r in p.get("refEntities", []) or []:
            m = TAIL.search(r.get("displayName", "") if isinstance(r, dict) else "")
            if m and m.group(1).upper() in allowed:
                out.add(m.group(1).upper())
    return out


def pathway_graph(name: str, members: Iterable[str],
                  max_reactions: int = 120) -> Tuple[nx.DiGraph | None, dict]:
    """Bipartite gene/reaction graph laid out in running order."""
    allowed = {g.upper() for g in members}
    sid = pathway_id(name)
    if not sid:
        return None, {"reason": "name not found in Reactome"}
    ev = _get(f"/data/pathway/{sid}/containedEvents") or []
    rxn = [e for e in ev if isinstance(e, dict)
           and "Reaction" in str(e.get("className", ""))][:max_reactions]
    if not rxn:
        return None, {"reason": "no reactions", "stId": sid}

    def side(d, key, sub=None):
        out = set()
        for p in d.get(key) or []:
            if not isinstance(p, dict):
                continue
            if sub:
                p = p.get(sub) or {}
                if not isinstance(p, dict):
                    continue
            out |= _from_display(p.get("displayName", ""), allowed)
        return out

    info = {}
    for e in rxn:
        st = e["stId"]
        d = _get(f"/data/query/{st}")
        d = d if isinstance(d, dict) else {}
        info[st] = {
            "name": d.get("displayName", e.get("displayName", "")),
            "genes": _participants(st, allowed),
            "in": side(d, "input"),
            "out": side(d, "output"),
            "cat": side(d, "catalystActivity", "physicalEntity"),
            "pre": [p.get("stId") for p in (d.get("precedingEvent") or [])
                    if isinstance(p, dict)],
        }

    G = nx.DiGraph()
    for st, r in info.items():
        G.add_node(st, kind="reaction", label=r["name"])
    covered = set()
    for st, r in info.items():
        covered |= r["genes"]
        for g in r["genes"]:
            G.add_node(g, kind="gene")
        for g in r["in"] | r["cat"]:
            G.add_edge(g, st, kind="input")
        for g in r["out"]:
            G.add_edge(st, g, kind="output")
        for g in r["genes"] - r["in"] - r["out"] - r["cat"]:
            G.add_edge(st, g, kind="member")      # side unresolved
        for p in r["pre"]:
            if p in info:
                G.add_edge(p, st, kind="precedes")

    meta = {"stId": sid, "n_reactions": len(info),
            "n_precedence": sum(1 for _, _, k in G.edges(data="kind") if k == "precedes"),
            "coverage": len(covered) / max(len(allowed), 1),
            "n_genes": len(covered), "n_members": len(allowed)}
    return G, meta


def usable(meta: dict, min_reactions: int = 4, min_coverage: float = 0.35,
           min_precedence: int = 2) -> bool:
    """Is there an order worth drawing, or should this fall back to STRING?"""
    return bool(meta) and meta.get("n_reactions", 0) >= min_reactions \
        and meta.get("coverage", 0) >= min_coverage \
        and meta.get("n_precedence", 0) >= min_precedence


def layers(G: nx.DiGraph) -> Dict[str, int]:
    """Running order: reactions get their depth in the precedence DAG, genes
    inherit the earliest reaction they take part in.

    Reactome pathways contain cycles, so the precedence graph is condensed into
    its strongly connected components before the order is taken.
    """
    rx = [n for n, k in G.nodes(data="kind") if k == "reaction"]
    R = nx.DiGraph()
    R.add_nodes_from(rx)
    R.add_edges_from([(a, b) for a, b, k in G.edges(data="kind") if k == "precedes"])
    depth = {}
    if len(R):
        C = nx.condensation(R)
        d = {}
        for comp in nx.topological_sort(C):
            pred = list(C.predecessors(comp))
            d[comp] = 0 if not pred else 1 + max(d[p] for p in pred)
        for comp, mem in C.nodes(data="members"):
            for n in mem:
                depth[n] = d[comp]
    out = dict(depth)
    for g, k in G.nodes(data="kind"):
        if k != "gene":
            continue
        near = sorted(depth.get(n, 0) for n in
                      list(G.predecessors(g)) + list(G.successors(g))
                      if G.nodes[n].get("kind") == "reaction")
        # the median, not the earliest: a scaffold present in every step would
        # otherwise pile every gene into the first column
        out[g] = near[len(near) // 2] if near else 0
    return out
