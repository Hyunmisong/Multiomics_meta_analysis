"""Gastric vs oesophageal: pool the cohorts, or analyse them apart?

Four analyses, in the order they constrain the decision:

  1  Proteome cohort clustering -- ANOSIM / PERMANOVA on cohort profiles.
     Fully powered (4 + 4 cohorts, 35 labellings, min p = 0.029).
  2  Pooled results per tumour type per layer -- Stouffer on p-values, plus a
     direction vote so the cohort that ships no p-values still counts.
  3  Layer concordance 2x2 -- does the genome of one tumour type agree more
     with the proteome of the same type than with the other's?
  4  Power-matched discovery -- is a tumour-type-pure set of cohorts more
     consistent than a mixed set of the same size?

Outputs land in results/ and docs/tissue_split.md.
"""
import _bootstrap  # noqa: F401

import json
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pmap import enrich, genesets, pool, tissue

from pmap.paths import FIGURES, RESULTS, ROOT, load_config

DOCS = ROOT / "docs"
PROT_SCHEME = "fdr05"
ALPHA = 0.05
MIN_COHORTS = 2          # a gene must be measured in >= 2 cohorts of the group
# Pooled gene sets are size-matched rather than FDR-thresholded: see
# pool.pooled_gene_set. TOP_N is the primary size; the others are the
# sensitivity sweep that shows whether the conclusion depends on it.
TOP_N = 1000
TOP_N_SWEEP = [250, 500, 1000, 2000]


def tissue_map(long: pd.DataFrame) -> pd.DataFrame:
    return (long[["layer", "dataset", "tissue"]].drop_duplicates()
            .reset_index(drop=True))


# ----------------------------------------------------------------- analysis 1
def cohort_clustering(ora, tmap, layer, scheme, out):
    prof = tissue.enrichment_profile(ora, layer, scheme)
    if prof.shape[1] < 4:
        print(f"  {layer}: only {prof.shape[1]} cohorts with a profile, skipped")
        return None
    D = tissue.profile_distance(prof)
    lab = (tmap[tmap.layer == layer].set_index("dataset")["tissue"]
           .reindex(D.index))
    if lab.nunique() < 2 or lab.value_counts().min() < 2:
        print(f"  {layer}: one tumour type has <2 cohorts, no test possible")
        return None

    res = {"layer": layer, "scheme": scheme, "n_cohorts": len(D),
           "n_pathways": len(prof)}
    res.update({f"anosim_{k}": v for k, v in tissue.anosim(D, lab).items()})
    res.update({f"permanova_{k}": v for k, v in tissue.permanova(D, lab).items()})
    sil = tissue.silhouette(D, lab)
    res["mean_silhouette"] = float(sil["silhouette"].mean())

    D.to_csv(out / f"tissue_distance_{layer}.csv")
    sil.to_csv(out / f"tissue_silhouette_{layer}.csv", index=False)
    co = tissue.pcoa(D)
    co.insert(0, "tissue", lab.values)
    co.to_csv(out / f"tissue_pcoa_{layer}.csv")

    print(f"  {layer}: ANOSIM R={res['anosim_R']:+.3f} p={res['anosim_p']:.3f} | "
          f"PERMANOVA R2={res['permanova_R2']:.3f} p={res['permanova_p']:.3f} | "
          f"silhouette {res['mean_silhouette']:+.3f} "
          f"(min attainable p = {res['anosim_min_attainable_p']:.3f})")
    return res


# ----------------------------------------------------------------- analysis 2
def pooled_layer(long, gcalls, universes, layer, cohorts, one_sided, scheme):
    """Stouffer + vote for one group of cohorts of one layer."""
    per = pool.gene_stats(long[long["scheme"] == scheme], layer, cohorts)
    st = pool.stouffer(per, one_sided=one_sided, min_cohorts=MIN_COHORTS)
    calls = gcalls[(gcalls["layer"] == layer) & (gcalls["scheme"] == scheme)
                   & (gcalls["dataset"].isin(cohorts))]
    vt = pool.vote(per, calls, universes, layer, min_cohorts=MIN_COHORTS)
    return st, vt


def main():
    DOCS.mkdir(exist_ok=True)
    long = pd.read_parquet(RESULTS / "features_long.parquet")
    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    ora = pd.read_parquet(RESULTS / "ora.parquet")
    gs = pd.read_csv(RESULTS / "genesets_long.csv.gz")
    sets = genesets.geneset_dict(gs)
    with open(RESULTS / "universes.json") as fh:
        universes = {k: set(v) for k, v in json.load(fh).items()}
    cfg = load_config("analysis.yaml")
    gen_scheme = (cfg.get("cross_layer_scheme") or {}).get("genome", "g_top200_n2")

    tmap = tissue_map(long)
    groups = {}
    for layer in ["proteome", "genome"]:
        for t in ["gastric", "esophageal"]:
            groups[(layer, t)] = sorted(
                tmap[(tmap.layer == layer) & (tmap.tissue == t)]["dataset"])

    print("cohorts per group:")
    for k, v in groups.items():
        print(f"  {k[0]:9s} {k[1]:11s} {v}")

    # ------------------------------------------------------------- analysis 1
    print("\n== 1. cohort clustering by tumour type ==")
    clus = [r for r in [cohort_clustering(ora, tmap, "proteome", PROT_SCHEME, RESULTS),
                        cohort_clustering(ora, tmap, "genome", gen_scheme, RESULTS)]
            if r]
    pd.DataFrame(clus).to_csv(RESULTS / "tissue_clustering.csv", index=False)

    # ------------------------------------------------------------- analysis 2
    print("\n== 2. pooled result per tumour type ==")
    pooled, shared_universe = {}, {}
    for layer, scheme, one_sided in [("proteome", PROT_SCHEME, False),
                                     ("genome", gen_scheme, True)]:
        unis = {}
        for t in ["gastric", "esophageal"]:
            arms = [f"{layer}|{d}|na" for d in groups[(layer, t)]]
            unis[t] = pool.pooled_universe(universes, arms, MIN_COHORTS)
        shared = unis["gastric"] & unis["esophageal"]
        shared_universe[layer] = shared
        print(f"  {layer}: gastric universe {len(unis['gastric']):,}, "
              f"esophageal {len(unis['esophageal']):,}, shared {len(shared):,}")
        for t in ["gastric", "esophageal"]:
            st, vt = pooled_layer(long, gcalls, universes, layer,
                                  groups[(layer, t)], one_sided, scheme)
            pooled[(layer, t)] = st
            st.to_csv(RESULTS / f"pooled_{layer}_{t}_stouffer.csv", index=False)
            if len(vt):
                vt.to_csv(RESULTS / f"pooled_{layer}_{t}_vote.csv", index=False)
            n_st = int((st["padj_pooled"] <= ALPHA).sum())
            n_vt = int((vt["padj_vote"] <= ALPHA).sum()) if len(vt) else 0
            print(f"    {t:11s} cohorts={len(groups[(layer, t)])} "
                  f"Stouffer genes={n_st:,} | vote genes={n_vt:,}")
        # all cohorts together, for the pooled-vs-split comparison
        allc = groups[(layer, "gastric")] + groups[(layer, "esophageal")]
        st, _ = pooled_layer(long, gcalls, universes, layer, allc, one_sided, scheme)
        pooled[(layer, "all")] = st
        st.to_csv(RESULTS / f"pooled_{layer}_all_stouffer.csv", index=False)
        print(f"    {'all':11s} cohorts={len(allc)} "
              f"Stouffer genes={int((st['padj_pooled'] <= ALPHA).sum()):,}")

    # ---- enrichment of each pooled gene set, on the shared background -------
    print("\n== 3. enrichment of the pooled sets (size-matched) ==")
    sweeps = {}
    for n_top in TOP_N_SWEEP:
        rows = []
        for (layer, t), st in pooled.items():
            uni = shared_universe[layer]
            hits = pool.pooled_gene_set(st, uni, top_n=n_top)
            if len(hits) < 50:
                continue
            res = enrich.run_ora(
                pd.DataFrame({"arm": f"{layer}|{t}", "layer": layer, "dataset": t,
                              "norm": "pooled", "scheme": "pooled",
                              "gene": sorted(hits), "direction": "any"}),
                {f"{layer}|{t}": uni}, sets, {"any": ["any"]}, min_overlap=3)
            res["group"] = t
            res["top_n"] = n_top
            rows.append(res)
            if n_top == TOP_N:
                print(f"  {layer:9s} {t:11s} {len(hits):,} genes -> "
                      f"{int((res['padj'] <= ALPHA).sum()):,} enriched pathways")
        sweeps[n_top] = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    pooled_ora = sweeps[TOP_N]
    pd.concat(sweeps.values(), ignore_index=True).to_csv(
        RESULTS / "pooled_ora.csv.gz", index=False)

    # ------------------------------------------- analysis 3: layer concordance
    print("\n== 4. layer concordance (same tumour type vs crossed) ==")
    allc = []
    for n_top, po in sweeps.items():
        c = concordance(po)
        if c is None:
            continue
        c["top_n"] = n_top
        c["delta"] = c.attrs["delta"]
        allc.append(c)
    conc = pd.concat(allc, ignore_index=True) if allc else None
    if conc is not None:
        conc.to_csv(RESULTS / "tissue_concordance.csv", index=False)
        print(conc[["top_n", "genome_group", "proteome_group", "n_pathways",
                    "spearman", "matched"]].to_string(index=False))
        print("\n  delta (matched - crossed) by set size:")
        for n_top, d in conc.groupby("top_n")["delta"].first().items():
            print(f"    top {n_top:5d}: {d:+.3f}")

    over = overlap_table(pd.concat(sweeps.values(), ignore_index=True), TOP_N)
    over.to_csv(RESULTS / "tissue_pathway_overlap.csv", index=False)
    print("\n== 5. pooled vs split, enriched pathway sets ==")
    print(over.to_string(index=False))
    make_figure(clus, conc, pooled_ora, over)

    write_report(clus, conc, groups, shared_universe, pooled, over)
    print(f"\nreport -> {(DOCS / 'tissue_split.md').relative_to(ROOT)}")


def overlap_table(pooled_ora: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Enriched pathway sets: pooled, per tumour type, and what splitting costs."""
    po = pooled_ora[pooled_ora["top_n"] == top_n]
    rows = []
    for layer, s in po.groupby("layer"):
        sets = {}
        for grp, t in s.groupby("group"):
            sig = t[t["padj"] <= ALPHA]
            sets[grp] = set(sig["collection"] + " | " + sig["pathway"])
        a, g, e = sets.get("all", set()), sets.get("gastric", set()), sets.get("esophageal", set())
        rows.append({"layer": layer, "pooled": len(a), "gastric": len(g),
                     "esophageal": len(e), "shared_by_both_types": len(g & e),
                     "gastric_only": len(g - e), "esophageal_only": len(e - g),
                     "jaccard": len(g & e) / max(len(g | e), 1),
                     "gained_by_splitting": len((g | e) - a),
                     "lost_by_splitting": len(a - (g | e))})
    return pd.DataFrame(rows)


def make_figure(clus, conc, pooled_ora, over):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    COL = {"gastric": "#c0392b", "esophageal": "#2a6fb0"}
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    fig.suptitle("Gastric vs oesophageal: pool the cohorts, or split them?",
                 fontsize=13, fontweight="bold", x=0.01, ha="left")

    # A: PCoA of proteome cohorts
    ax = axes[0, 0]
    try:
        co = pd.read_csv(RESULTS / "tissue_pcoa_proteome.csv", index_col=0)
        for t, s in co.groupby("tissue"):
            ax.scatter(s["PCo1"], s["PCo2"], s=120, color=COL[t], label=t,
                       edgecolor="white", linewidth=1.5, zorder=3)
            for c, r in s.iterrows():
                ax.annotate(c, (r["PCo1"], r["PCo2"]), fontsize=8,
                            xytext=(6, 4), textcoords="offset points")
        r = [x for x in clus if x["layer"] == "proteome"]
        if r:
            r = r[0]
            ax.set_title(f"A  Proteome cohorts in pathway space\n"
                         f"PERMANOVA R$^2$={r['permanova_R2']:.2f}, p={r['permanova_p']:.3f} "
                         f"(floor {r['permanova_min_attainable_p']:.3f})",
                         loc="left", fontsize=10.5, fontweight="bold")
        ax.set_xlabel("PCo1")
        ax.set_ylabel("PCo2")
        ax.legend(frameon=False, fontsize=9)
    except FileNotFoundError:
        ax.axis("off")

    # B: silhouette per cohort
    ax = axes[0, 1]
    try:
        sil = pd.read_csv(RESULTS / "tissue_silhouette_proteome.csv").sort_values("silhouette")
        ax.barh(range(len(sil)), sil["silhouette"],
                color=[COL[t] for t in sil["tissue"]])
        ax.set_yticks(range(len(sil)))
        ax.set_yticklabels(sil["cohort"])
        ax.axvline(0, color="#555", lw=0.8)
        ax.set_xlabel("silhouette (negative = closer to the other tumour type)")
        ax.set_title("B  How cleanly each cohort belongs", loc="left",
                     fontsize=10.5, fontweight="bold")
    except FileNotFoundError:
        ax.axis("off")

    # C: enriched pathway sets
    ax = axes[1, 0]
    if len(over):
        y = range(len(over))
        ax.barh(y, over["gastric_only"], color=COL["gastric"], label="gastric only")
        ax.barh(y, over["shared_by_both_types"], left=over["gastric_only"],
                color="#7d4a9e", label="both types")
        ax.barh(y, over["esophageal_only"],
                left=over["gastric_only"] + over["shared_by_both_types"],
                color=COL["esophageal"], label="oesophageal only")
        for i, r in enumerate(over.itertuples(index=False)):
            ax.text(r.gastric_only + r.shared_by_both_types + r.esophageal_only + 4, i,
                    f"Jaccard {r.jaccard:.3f}", va="center", fontsize=8.5)
        ax.set_yticks(list(y))
        ax.set_yticklabels(over["layer"])
        ax.set_xlabel("enriched pathways (size-matched gene sets)")
        ax.legend(frameon=False, fontsize=8.5, ncol=3)
        ax.set_title("C  The two tumour types barely share a pathway answer",
                     loc="left", fontsize=10.5, fontweight="bold")

    # D: concordance delta across set sizes
    ax = axes[1, 1]
    if conc is not None:
        d = conc.groupby("top_n")["delta"].first()
        ax.plot(d.index, d.values, "o-", color="#5b6472")
        ax.axhline(0, color="#999", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_xticks(list(d.index))
        ax.set_xticklabels([str(x) for x in d.index])
        ax.set_xlabel("pooled gene-set size (top N by |Z|)")
        ax.set_ylabel("Spearman delta\n(matched - crossed)")
        ax.set_ylim(-0.25, 0.25)
        ax.set_title("D  Genome and proteome of the same type do NOT\n"
                     "agree more than crossed pairs",
                     loc="left", fontsize=10.5, fontweight="bold")

    for a in axes.ravel():
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"tissue_split.{ext}", dpi=200, bbox_inches="tight")
    print(f"figure -> {(FIGURES / 'tissue_split.png').relative_to(ROOT)}")


def concordance(pooled_ora: pd.DataFrame):
    """2x2 Spearman between each layer's pooled enrichment profiles."""
    if pooled_ora.empty:
        return None
    p = pooled_ora.copy()
    p["key"] = p["collection"] + " | " + p["pathway"]
    p["score"] = -np.log10(p["pval"].clip(lower=1e-300))
    wide = p.pivot_table(index="key", columns=["layer", "group"], values="score")
    cols = [c for c in wide.columns if c[1] in ("gastric", "esophageal")]
    wide = wide[cols].dropna(how="all")
    rows = []
    for gl in [c for c in cols if c[0] == "genome"]:
        for pl in [c for c in cols if c[0] == "proteome"]:
            both = wide[[gl, pl]].fillna(0.0)
            if len(both) < 10:
                continue
            r = spearmanr(both[gl], both[pl]).statistic
            rows.append({"genome_group": gl[1], "proteome_group": pl[1],
                         "n_pathways": len(both), "spearman": float(r),
                         "matched": gl[1] == pl[1]})
    if not rows:
        return None
    df = pd.DataFrame(rows)
    m = df[df.matched]["spearman"].mean()
    x = df[~df.matched]["spearman"].mean()
    df.attrs["delta"] = m - x
    return df


def write_report(clus, conc, groups, shared_universe, pooled, over):
    L = ["# Gastric vs oesophageal: pool or split?", "",
         "Regenerated by `scripts/10_tissue_split.py`; do not edit by hand.", "",
         "## Cohorts", "", "| layer | tumour type | cohorts |", "|---|---|---|"]
    for (layer, t), v in groups.items():
        L.append(f"| {layer} | {t} | {', '.join(v)} ({len(v)}) |")
    L += ["", "## 1. Does tumour type explain how the cohorts cluster?", ""]
    if clus:
        c = pd.DataFrame(clus)
        L += ["| layer | cohorts | pathways | ANOSIM R | p | PERMANOVA R2 | p | "
              "mean silhouette | min attainable p |", "|---|---|---|---|---|---|---|---|---|"]
        for r in clus:
            L.append(f"| {r['layer']} | {r['n_cohorts']} | {r['n_pathways']} | "
                     f"{r['anosim_R']:+.3f} | {r['anosim_p']:.3f} | "
                     f"{r['permanova_R2']:.3f} | {r['permanova_p']:.3f} | "
                     f"{r['mean_silhouette']:+.3f} | {r['anosim_min_attainable_p']:.3f} |")
        L += ["", "`R2` is the share of variation the tumour-type label explains. "
                  "A small `R2` with a significant `p` means the split is real but "
                  "not worth the cohorts it costs. Negative silhouette for a cohort "
                  "means it sits closer to the other tumour type than to its own "
                  "(see `results/tissue_silhouette_*.csv`)."]
    L += ["", "## 2. Shared background used for the pooled comparison", "",
          "| layer | genes measured in >= 2 cohorts of BOTH types |", "|---|---|"]
    for layer, u in shared_universe.items():
        L.append(f"| {layer} | {len(u):,} |")
    L += ["", "## 3. Layer concordance", ""]
    if conc is not None:
        L += ["Pooled gene sets are size-matched (top N by |Z|), because "
              "meta-analysing four proteome cohorts puts ~70% of the universe "
              "under FDR 0.05 and a set that large cannot be tested for "
              "over-representation.", "",
              "| set size | genome | proteome | pathways | Spearman | matched |",
              "|---|---|---|---|---|---|"]
        for r in conc.itertuples(index=False):
            L.append(f"| {r.top_n} | {r.genome_group} | {r.proteome_group} | "
                     f"{r.n_pathways} | {r.spearman:+.3f} | "
                     f"{'yes' if r.matched else 'no'} |")
        L += ["", "| set size | delta (matched - crossed) |", "|---|---|"]
        for n_top, d in conc.groupby("top_n")["delta"].first().items():
            L.append(f"| {n_top} | {d:+.3f} |")
        L += ["", "Positive delta means the genome and the proteome of the same "
                  "tumour type agree more than crossed pairs, which is the "
                  "evidence for splitting. A delta that flips sign across set "
                  "sizes is noise, not evidence."]
    L += ["", "## 4. What splitting costs and buys", "",
          "Enriched pathway sets from size-matched pooled gene sets "
          f"(top {TOP_N} genes by |Z|).", "",
          "| layer | pooled | gastric | oesophageal | both types | gastric only | "
          "oesophageal only | Jaccard | gained by splitting | lost by splitting |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in over.itertuples(index=False):
        L.append(f"| {r.layer} | {r.pooled} | {r.gastric} | {r.esophageal} | "
                 f"{r.shared_by_both_types} | {r.gastric_only} | {r.esophageal_only} | "
                 f"{r.jaccard:.3f} | {r.gained_by_splitting} | {r.lost_by_splitting} |")
    L += ["", "## Caveat that bounds every number above", "",
          "Tumour type is completely confounded with study: the four gastric "
          "proteome cohorts are four different studies, and so are the four "
          "oesophageal ones. Nothing in this design separates *gastric vs "
          "oesophageal* from *these four laboratories vs those four*. Sequencing "
          "depth is **not** confounded -- the measured-gene range overlaps between "
          "the types (gastric 3,903-11,340, oesophageal 1,714-9,325) -- but batch "
          "is. The split is justified on reproducibility grounds; calling it "
          "biology needs cohorts of both types from the same pipeline."]
    (DOCS / "tissue_split.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
