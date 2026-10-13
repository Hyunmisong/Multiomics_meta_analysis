"""Write the threshold-scheme comparison as tracked Markdown.

The point of running every cutoff instead of picking one is that the choice is
then evidence, not taste. That evidence has to be visible on GitHub without
cloning the repo and running anything, so it is written to a tracked file and
spliced into the README between markers, regenerated on every run.

Outputs:
  docs/scheme_comparison.md        full tables, every layer and collection
  README.md                        headline tables, between the SCHEME markers
"""
import _bootstrap  # noqa: F401

import pandas as pd

from pmap.paths import ROOT, RESULTS, load_config

DOCS = ROOT / "docs"
OUT = DOCS / "scheme_comparison.md"
START, END = "<!-- SCHEME_TABLES:START -->", "<!-- SCHEME_TABLES:END -->"


def md_table(df: pd.DataFrame, floatfmt="{:.3f}") -> str:
    def cell(v):
        if isinstance(v, float):
            return "-" if pd.isna(v) else floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(cell(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, rule, *rows])


def scheme_definitions(cfg: dict) -> str:
    lines = ["| scheme | layer | rule |", "|---|---|---|"]
    for name, sc in cfg["threshold_schemes"].items():
        bits = []
        if "padj_max" in sc:
            bits.append(f"BH-adjusted p ≤ {sc['padj_max']}")
        if "pval_max" in sc:
            bits.append(f"raw p ≤ {sc['pval_max']}")
        if sc.get("abs_logfc_min"):
            bits.append(f"\\|log2FC\\| ≥ {sc['abs_logfc_min']}")
        if sc.get("mode") == "rank":
            f = sc["tail_frac"]
            bits.append(f"log2FC in the bottom {f:.0%} (down) or top {f:.0%} (up) of that arm")
        lines.append(f"| `{name}` | proteome, phospho | {' and '.join(bits)} |")
    for name, sc in (cfg.get("genome_schemes") or {}).items():
        bits = []
        if "q_max" in sc:
            bits.append(f"MutSig2CV q ≤ {sc['q_max']}")
        if "max_rank" in sc:
            bits.append(f"MutSig2CV rank ≤ {sc['max_rank']}")
        bits.append(f"≥ {sc.get('min_variants', 1)} truncating variants")
        if sc.get("min_strong"):
            bits.append(f"≥ {sc['min_strong']} of them destroying ≥ 50% of the protein")
        lines.append(f"| `{name}` | genome | {' and '.join(bits)} |")
    return "\n".join(lines)


def stability_table(stab: pd.DataFrame, layer: str, norm: str | None = None) -> pd.DataFrame:
    s = stab[stab["layer"] == layer]
    if s.empty:
        return pd.DataFrame()
    piv = s[s["direction"] != "any"].pivot_table(
        index="scheme", columns=["collection", "direction"],
        values="mean_pairwise_jaccard")
    if piv.empty:  # genome has only the down direction
        piv = s.pivot_table(index="scheme", columns="collection",
                            values="mean_pairwise_jaccard")
    piv.columns = [" ".join(c) if isinstance(c, tuple) else str(c) for c in piv.columns]
    med = s.groupby("scheme")["median_n_enriched"].median().rename("median sets/arm")
    out = piv.join(med).reset_index().rename(columns={"scheme": "scheme"})
    return out


def concentration_counts(conc: pd.DataFrame, layer: str) -> pd.DataFrame:
    c = conc[(conc["layer"] == layer) & (conc["conc_padj"] <= 0.05)]
    if c.empty:
        return pd.DataFrame()
    t = (c.groupby(["scheme", "norm", "direction"]).size().rename("pathways")
         .reset_index())
    return t.pivot_table(index=["scheme", "norm"], columns="direction",
                         values="pathways", fill_value=0).reset_index()


def main():
    DOCS.mkdir(exist_ok=True)
    cfg = load_config("analysis.yaml")
    stab = pd.read_csv(RESULTS / "scheme_stability.csv")
    conc = pd.read_csv(RESULTS / "concentration.csv")
    arms = pd.read_csv(RESULTS / "arm_summary.csv")

    chosen = cfg.get("cross_layer_scheme") or {}
    headline = []

    headline.append("#### Continuous layers (proteome, phosphoproteome)\n")
    cont = stability_table(stab, "proteome")
    headline.append(md_table(cont))
    headline.append("\n#### Genome layer\n")
    gen = stability_table(stab, "genome")
    gconc = conc[(conc["layer"] == "genome") & (conc["conc_padj"] <= 0.05)]
    gcount = gconc[gconc["direction"] == "down"].groupby("scheme").size()
    if not gen.empty:
        gen["concentrated pathways"] = gen["scheme"].map(gcount).fillna(0).astype(int)
        headline.append(md_table(gen))
    headline.append("")
    headline.append("Higher mean pairwise Jaccard means the cohorts agree more. "
                    "Read it next to `median sets/arm`: a scheme that calls every "
                    "pathway enriched scores well for the wrong reason.")
    headline.append("")
    for layer, name in chosen.items():
        headline.append(f"Chosen for cross-layer views - **{layer}: `{name}`**.")
    headline_md = "\n".join(headline)

    # ---------------------------------------------------------------- docs
    doc = ["# Threshold-scheme comparison", "",
           "Regenerated by `scripts/09_report_schemes.py`; do not edit by hand.", "",
           "Every cutoff is run on every cohort, and the one that reproduces best "
           "across independent cohorts is chosen. Nothing here is a preference.",
           "", "## Scheme definitions", "", scheme_definitions(cfg), "",
           "## Cross-cohort reproducibility", "",
           "For one (layer, collection, scheme, direction): take the set of "
           "pathways enriched at BH-adjusted p <= 0.05 in each cohort, then average "
           "the Jaccard index over every pair of cohorts.", ""]
    for layer in ["proteome", "phospho", "genome"]:
        t = stability_table(stab, layer)
        if t.empty:
            continue
        doc += [f"### {layer}", "", md_table(t), ""]
    doc += ["## Cross-cohort concentration (pathways at FDR <= 0.05)", ""]
    for layer in ["proteome", "phospho", "genome"]:
        t = concentration_counts(conc, layer)
        if t.empty:
            continue
        doc += [f"### {layer}", "", md_table(t, floatfmt="{:.0f}"), ""]
    doc += ["## Genes called per cohort", "", md_table(
        arms.sort_values(["scheme", "layer", "dataset"])[
            ["scheme", "layer", "dataset", "norm", "n_genes_called", "n_up", "n_down"]],
        floatfmt="{:.0f}"), ""]
    OUT.write_text("\n".join(doc) + "\n", encoding="utf-8")

    # -------------------------------------------------------------- README
    rp = ROOT / "README.md"
    txt = rp.read_text(encoding="utf-8")
    if START in txt and END in txt:
        pre, rest = txt.split(START, 1)
        _, post = rest.split(END, 1)
        rp.write_text(pre + START + "\n" + headline_md + "\n" + END + post,
                      encoding="utf-8")
        print("README scheme tables updated")
    else:
        print("README markers not found; only docs/ written")
    print(f"written -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
