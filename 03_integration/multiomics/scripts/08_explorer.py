"""Interactive pathway explorer: tick the layers, see which pathways light up.

Writes one self-contained HTML file (data embedded, no server, no network) so
it can be opened from the repo or dropped into a GitHub Pages branch.

The page answers: given a chosen set of omics layers and a required level of
agreement, which pathways are supported by all of them at once? Layers with no
data yet (genomics, scRNA-seq) appear as disabled checkboxes, so the panel
already shows the shape of the finished figure.

    python scripts/08_explorer.py
"""
import _bootstrap  # noqa: F401

import json

import pandas as pd

from pmap import cellmap
from pmap.paths import FIGURES, RESULTS, load_config

OUT = FIGURES / "pathway_explorer.html"


def build_records(ora: pd.DataFrame, rules, default) -> tuple:
    """Per (scheme, collection, pathway, arm, direction): which cohorts hit it.

    Both bars are kept -- the per-cohort FDR and the nominal p -- because which
    one you use changes the answer a lot when a layer has only three cohorts,
    and the page lets the reader switch.
    """
    keep = ora[(ora["pval"] <= 0.05)].copy()
    keep["arm"] = keep["layer"] + ("|" + keep["norm"]).where(keep["norm"] != "na", "")

    recs = {}
    for (scheme, coll, pw, arm, d), sub in keep.groupby(
            ["scheme", "collection", "pathway", "arm", "direction"]):
        key = (scheme, coll, pw)
        r = recs.setdefault(key, {})
        slot = r.setdefault(arm, {})
        slot[d] = {
            "f": sorted(sub.loc[sub["padj"] <= 0.05, "dataset"].unique().tolist()),
            "n": sorted(sub["dataset"].unique().tolist()),
        }

    out = []
    for (scheme, coll, pw), arms in recs.items():
        out.append({
            "s": scheme, "c": coll, "p": pw,
            "comp": cellmap.assign_compartment(pw, rules, default),
            "a": arms,
        })
    totals = (ora.groupby(ora["layer"] + ("|" + ora["norm"]).where(ora["norm"] != "na", ""))
              ["dataset"].nunique().to_dict())
    return out, totals


HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pathway Evidence Explorer</title>
<style>
  :root {
    --bg:#ffffff; --panel:#f7f8fa; --line:#dfe3e9; --ink:#1d2129; --muted:#6b7280;
    --up:#c0392b; --down:#2a6fb0; --both:#7d4a9e; --none:#e9ecef;
    --nucleus:#eef1f8; --nucleus-s:#6f7da3; --mito:#f3eef8; --mito-s:#8d7aa8;
    --er:#eef6f2; --er-s:#6f9f88; --ecm:#f8f3ee; --ecm-s:#a8907a;
    --mem:#eef5f8; --mem-s:#7a98a8; --cyt:#fcfcfd; --cyt-s:#b7bdc7;
  }
  :root:not([data-theme="light"]) { }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --bg:#14171c; --panel:#1b1f26; --line:#2b313b; --ink:#e8eaed; --muted:#9aa3b0;
      --none:#2b313b;
      --nucleus:#1c2130; --nucleus-s:#7b89ad; --mito:#241d2e; --mito-s:#9d8ab8;
      --er:#17241f; --er-s:#7fae98; --ecm:#241f19; --ecm-s:#b8a08a;
      --mem:#172028; --mem-s:#8aa8b8; --cyt:#181c22; --cyt-s:#4a5058;
    }
  }
  :root[data-theme="dark"] {
    --bg:#14171c; --panel:#1b1f26; --line:#2b313b; --ink:#e8eaed; --muted:#9aa3b0;
    --none:#2b313b;
    --nucleus:#1c2130; --nucleus-s:#7b89ad; --mito:#241d2e; --mito-s:#9d8ab8;
    --er:#17241f; --er-s:#7fae98; --ecm:#241f19; --ecm-s:#b8a08a;
    --mem:#172028; --mem-s:#8aa8b8; --cyt:#181c22; --cyt-s:#4a5058;
  }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--bg); color:var(--ink);
         font:14px/1.45 Inter,-apple-system,Segoe UI,Helvetica,Arial,sans-serif; }
  .wrap { display:flex; gap:0; min-height:100vh; align-items:stretch; }
  aside { width:290px; flex:0 0 290px; background:var(--panel);
          border-right:1px solid var(--line); padding:18px 16px; }
  main { flex:1; padding:18px 16px 40px; min-width:0; }
  h1 { font-size:19px; margin:0 0 2px; }
  .sub { color:var(--muted); font-size:12px; margin:0 0 16px; }
  fieldset { border:0; padding:0; margin:0 0 16px; }
  legend { font-size:11px; font-weight:700; letter-spacing:.6px;
           text-transform:uppercase; color:var(--muted); margin-bottom:7px; }
  label.row { display:flex; align-items:center; gap:8px; padding:4px 0; cursor:pointer; }
  label.row.off { opacity:.45; cursor:not-allowed; }
  label.row .tag { margin-left:auto; font-size:11px; color:var(--muted); }
  select, input[type=range] { width:100%; }
  select { padding:5px 7px; border:1px solid var(--line); border-radius:6px;
           background:var(--bg); color:var(--ink); font-size:13px; }
  .count { font-variant-numeric:tabular-nums; }
  .comp { border:1px solid; border-radius:12px; padding:10px 12px 12px; margin:0 0 14px; }
  .comp > h2 { font-size:11px; font-weight:700; letter-spacing:.7px;
               text-transform:uppercase; margin:0 0 9px; }
  .chips { display:flex; flex-wrap:wrap; gap:8px; }
  .chip { border:1px solid var(--line); background:var(--bg); border-radius:8px;
          padding:7px 9px; min-width:190px; max-width:260px; }
  .chip .nm { font-size:12px; font-weight:600; line-height:1.25; }
  .chip .tracks { margin-top:6px; display:grid; gap:3px; }
  .tr { display:flex; align-items:center; gap:5px; }
  .tr .lb { font-size:9px; color:var(--muted); width:40px; flex:0 0 40px; }
  .dot { width:9px; height:9px; border-radius:2px; background:var(--none);
         border:1px solid var(--line); }
  .dot.up { background:var(--up); border-color:var(--up); }
  .dot.down { background:var(--down); border-color:var(--down); }
  .dot.both { background:var(--both); border-color:var(--both); }
  .empty { color:var(--muted); font-size:13px; padding:28px 4px; }
  .legend { display:flex; flex-wrap:wrap; gap:12px; font-size:11px;
            color:var(--muted); margin:4px 0 14px; align-items:center; }
  .legend span.k { display:inline-block; width:10px; height:10px; border-radius:2px;
                   margin-right:4px; vertical-align:-1px; }
  @media (max-width:860px) {
    .wrap { flex-direction:column; }
    aside { width:auto; flex:none; border-right:0; border-bottom:1px solid var(--line); }
  }
</style></head><body>
<div class="wrap">
<aside>
  <h1>Pathway evidence</h1>
  <p class="sub">Tick the layers. A pathway is shown when every ticked layer
  supports it in at least the required number of its own cohorts.</p>

  <fieldset><legend>Layers</legend><div id="layers"></div></fieldset>

  <fieldset><legend>Agreement</legend>
    <label class="row"><input type="checkbox" id="fdr" checked>
      per-cohort FDR &le; 0.05</label>
    <label class="row" style="display:block">
      min cohorts per layer: <b id="minlab" class="count">1</b>
      <input type="range" id="minc" min="1" max="4" value="1"></label>
    <label class="row"><input type="checkbox" id="signed">
      require one direction</label>
  </fieldset>

  <fieldset><legend>Gene sets</legend>
    <select id="coll"></select></fieldset>
  <fieldset><legend>Threshold scheme</legend>
    <select id="scheme"></select>
    <p class="sub" style="margin:6px 0 0">The genome layer has no effect size,
    so it is scored by its own scheme and is unaffected by this choice.</p>
  </fieldset>

  <p class="sub" id="status"></p>
</aside>
<main>
  <div class="legend">
    <span><span class="k" style="background:var(--up)"></span>up</span>
    <span><span class="k" style="background:var(--down)"></span>down</span>
    <span><span class="k" style="background:var(--both)"></span>both</span>
    <span><span class="k" style="background:var(--none)"></span>not enriched</span>
    <span>one dot = one cohort</span>
  </div>
  <div id="out"></div>
</main>
</div>
<script>
const DATA = __DATA__, TOTALS = __TOTALS__, ARMS = __ARMS__, COMPS = __COMPS__;
const SCHEME_FOR = __SCHEME_FOR__;   // layers scored by their own schemes
// index by scheme -> collection -> pathway, so each layer can be read under
// the scheme it was actually run with
const IDX = {};
DATA.forEach(d => {
  ((IDX[d.s] = IDX[d.s] || {})[d.c] = IDX[d.s][d.c] || {})[d.p] = d;
});
const schemeOf = arm => SCHEME_FOR[arm.split("|")[0]] || null;
const PATHWAYS = {};  // collection -> pathway -> compartment
DATA.forEach(d => { (PATHWAYS[d.c] = PATHWAYS[d.c] || {})[d.p] = d.comp; });
const $ = s => document.querySelector(s);

$("#layers").innerHTML = ARMS.map(a => {
  const n = TOTALS[a.key] || 0, off = n === 0;
  return `<label class="row ${off?'off':''}">
    <input type="checkbox" data-arm="${a.key}" ${off?'disabled':''}>
    ${a.label}<span class="tag">${off?'no data yet':n+' cohorts'}</span></label>`;
}).join("");

const ownSchemes = new Set(Object.values(SCHEME_FOR));
const schemes = [...new Set(DATA.map(d => d.s))].filter(s => !ownSchemes.has(s)).sort();
$("#scheme").innerHTML = schemes.map(s => `<option${s==="fdr05"?" selected":""}>${s}</option>`).join("");
const colls = [...new Set(DATA.map(d => d.c))].sort();
$("#coll").innerHTML = colls.map(c => `<option${c==="Hallmark"?" selected":""}>${c}</option>`).join("");

function cohorts(coll, pw, scheme, arm, dir, useFdr) {
  const s = schemeOf(arm) || scheme;
  const rec = (IDX[s] || {})[coll] && IDX[s][coll][pw];
  if (!rec) return [];
  const a = rec.a[arm]; if (!a) return [];
  const slot = a[dir]; if (!slot) return [];
  return useFdr ? slot.f : slot.n;
}

function render() {
  const picked = [...document.querySelectorAll("#layers input:checked")]
                   .map(i => i.dataset.arm);
  const useFdr = $("#fdr").checked, signed = $("#signed").checked;
  const minc = +$("#minc").value, scheme = $("#scheme").value, coll = $("#coll").value;
  $("#minlab").textContent = minc;

  const rows = Object.entries(PATHWAYS[coll] || {})
    .map(([p, comp]) => ({ p, comp }))
    .filter(rec => {
      if (!picked.length) return false;
      return picked.every(arm => {
        if (signed) {
          return cohorts(coll, rec.p, scheme, arm, "up", useFdr).length >= minc
              || cohorts(coll, rec.p, scheme, arm, "down", useFdr).length >= minc;
        }
        const u = cohorts(coll, rec.p, scheme, arm, "up", useFdr);
        const dn = cohorts(coll, rec.p, scheme, arm, "down", useFdr);
        const an = cohorts(coll, rec.p, scheme, arm, "any", useFdr);
        return new Set([...u, ...dn, ...an]).size >= minc;
      });
    });

  if (!picked.length) {
    $("#out").innerHTML = `<p class="empty">Tick at least one layer.</p>`;
    $("#status").textContent = ""; return;
  }
  $("#status").textContent =
    `${rows.length} pathway${rows.length===1?"":"s"} supported by `
    + `${picked.length} layer${picked.length===1?"":"s"}`
    + ` at ${useFdr ? "FDR\\u22640.05" : "nominal p<0.05"}, ${minc}+ cohort${minc===1?"":"s"} each`;

  if (!rows.length) {
    $("#out").innerHTML = `<p class="empty">Nothing clears this combination.
      Loosen the bar, lower the cohort requirement, or try a larger collection
      such as Reactome.</p>`;
    return;
  }

  const byComp = {};
  rows.forEach(r => (byComp[r.comp] = byComp[r.comp] || []).push(r));
  $("#out").innerHTML = COMPS.filter(c => byComp[c.key]).map(c => `
    <section class="comp" style="background:var(--${c.v});border-color:var(--${c.v}-s)">
      <h2 style="color:var(--${c.v}-s)">${c.label} <span class="count">(${byComp[c.key].length})</span></h2>
      <div class="chips">${byComp[c.key]
        .sort((a,b)=>a.p.localeCompare(b.p))
        .map(r => chip(coll, r, scheme, useFdr)).join("")}</div>
    </section>`).join("");
}

function chip(coll, rec, scheme, useFdr) {
  const tracks = ARMS.map(a => {
    const n = TOTALS[a.key] || 0;
    if (!n) return `<div class="tr"><span class="lb">${a.short}</span>
        <span style="font-size:9px;color:var(--muted)">no data yet</span></div>`;
    const up = new Set(cohorts(coll, rec.p, scheme, a.key, "up", useFdr));
    const dn = new Set(cohorts(coll, rec.p, scheme, a.key, "down", useFdr));
    const an = new Set(cohorts(coll, rec.p, scheme, a.key, "any", useFdr));
    const all = [...new Set([...up, ...dn, ...an])].sort();
    const dots = (a.cohorts || []).map(c => {
      const cls = up.has(c) && dn.has(c) ? "both" : up.has(c) ? "up"
                : dn.has(c) ? "down" : an.has(c) ? "both" : "";
      return `<span class="dot ${cls}" title="${a.short} ${c}"></span>`;
    }).join("");
    return `<div class="tr"><span class="lb">${a.short}</span>${dots}</div>`;
  }).join("");
  return `<div class="chip"><div class="nm">${rec.p}</div>
          <div class="tracks">${tracks}</div></div>`;
}

document.addEventListener("change", render);
document.addEventListener("input", e => { if (e.target.id === "minc") render(); });
["proteome"].forEach(k => {
  const el = document.querySelector(`#layers input[data-arm="${k}"]`);
  if (el) el.checked = true;
});
render();
</script></body></html>
"""


def main():
    ora = pd.read_parquet(RESULTS / "ora.parquet")
    cm = load_config("cellmap.yaml")
    rules, default = cm["compartments"]["rules"], cm["compartments"]["default"]

    scheme_for = load_config("analysis.yaml").get("cross_layer_scheme") or {}
    recs, totals = build_records(ora, rules, default)

    cohort_by_arm = (ora.assign(arm=ora["layer"] + ("|" + ora["norm"]).where(ora["norm"] != "na", ""))
                     .groupby("arm")["dataset"].unique().apply(lambda v: sorted(v)).to_dict())
    arms = [
        {"key": "proteome", "label": "Proteome", "short": "PROT"},
        {"key": "phospho|raw", "label": "Phospho, raw", "short": "P-raw"},
        {"key": "phospho|prot_adj", "label": "Phospho, protein-adjusted", "short": "P-adj"},
        {"key": "genome", "label": "Genomics", "short": "GEN"},
        {"key": "scrna", "label": "scRNA-seq", "short": "scRNA"},
    ]
    for a in arms:
        a["cohorts"] = cohort_by_arm.get(a["key"], [])

    comps = [
        {"key": "extracellular", "v": "ecm", "label": "Extracellular / secreted"},
        {"key": "membrane", "v": "mem", "label": "Plasma membrane & signalling"},
        {"key": "cytosol", "v": "cyt", "label": "Cytosol"},
        {"key": "er_golgi", "v": "er", "label": "ER / Golgi / vesicles"},
        {"key": "mitochondrion", "v": "mito", "label": "Mitochondrion"},
        {"key": "nucleus", "v": "nucleus", "label": "Nucleus"},
    ]

    html = (HTML
            .replace("__DATA__", json.dumps(recs, separators=(",", ":")))
            .replace("__TOTALS__", json.dumps(totals, separators=(",", ":")))
            .replace("__ARMS__", json.dumps(arms, separators=(",", ":")))
            .replace("__COMPS__", json.dumps(comps, separators=(",", ":")))
            .replace("__SCHEME_FOR__", json.dumps(scheme_for, separators=(",", ":"))))
    OUT.write_text(html, encoding="utf-8")
    print(f"{len(recs):,} pathway records -> {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)")
    print("cohorts per layer:", {k: len(v) for k, v in cohort_by_arm.items()})


if __name__ == "__main__":
    main()
