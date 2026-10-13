"""The subproject question: how much does protein-level normalisation of the
phosphoproteome change the pathway answer?

Outputs (results/):
  norm_site_level.csv       paired per-site comparison (PG4 only, full tables)
  norm_gene_level.csv       regulated-gene-set overlap per cohort
  norm_pathway_level.csv    enriched-pathway-set overlap per cohort
  norm_pathway_shifts.csv   pathways gained / lost on normalisation, listed
  scheme_stability.csv      which threshold scheme reproduces across cohorts
"""
import _bootstrap  # noqa: F401

import pandas as pd

from pmap import compare_norm
from pmap.paths import RESULTS


def pathway_shifts(ora: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """One row per pathway that changes enrichment status on normalisation."""
    ph = ora[ora["layer"] == "phospho"]
    keys = ["dataset", "collection", "scheme", "direction", "pathway"]
    wide = ph.pivot_table(index=keys, columns="norm",
                          values=["padj", "n_overlap", "odds_ratio"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    for c in ["padj_raw", "padj_prot_adj"]:
        if c not in wide:
            wide[c] = 1.0
    wide = wide.fillna({"padj_raw": 1.0, "padj_prot_adj": 1.0})
    sig_raw = wide["padj_raw"] <= alpha
    sig_adj = wide["padj_prot_adj"] <= alpha
    wide["status"] = "not enriched"
    wide.loc[sig_raw & sig_adj, "status"] = "robust"
    wide.loc[sig_raw & ~sig_adj, "status"] = "lost on adjustment"
    wide.loc[~sig_raw & sig_adj, "status"] = "gained on adjustment"
    return wide[wide["status"] != "not enriched"].sort_values(
        ["collection", "status", "padj_raw"]).reset_index(drop=True)


def main():
    long = pd.read_parquet(RESULTS / "features_long.parquet")
    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    ora = pd.read_parquet(RESULTS / "ora.parquet")

    site = compare_norm.site_level(long)
    site.to_csv(RESULTS / "norm_site_level.csv", index=False)
    print("\n== site level (paired) ==")
    print(site[["dataset", "scheme", "paired_comparison_valid", "n_sites_common",
                "n_sig_raw", "n_sig_adj", "n_sig_raw_only", "n_sig_adj_only",
                "n_direction_flips", "median_absFC_raw", "median_absFC_adj",
                "spearman_logFC"]].to_string(index=False))

    gene = compare_norm.gene_level(gcalls)
    gene.to_csv(RESULTS / "norm_gene_level.csv", index=False)
    print("\n== gene level ==")
    print(gene.to_string(index=False))

    pw = compare_norm.pathway_level(ora)
    pw.to_csv(RESULTS / "norm_pathway_level.csv", index=False)
    print("\n== pathway level (Jaccard of enriched sets, raw vs protein-adjusted) ==")
    print(pw[["dataset", "collection", "scheme", "direction", "n_enriched_raw",
              "n_enriched_adj", "n_shared", "jaccard", "lost_on_adjust",
              "gained_on_adjust", "spearman_neglog10padj"]].to_string(index=False))

    shifts = pathway_shifts(ora)
    shifts.to_csv(RESULTS / "norm_pathway_shifts.csv", index=False)
    print(f"\npathway shifts: {len(shifts):,} rows")
    print(shifts.groupby(["scheme", "status"]).size().to_string())

    stab = compare_norm.scheme_stability(ora)
    stab.to_csv(RESULTS / "scheme_stability.csv", index=False)
    print("\n== threshold-scheme stability (mean pairwise Jaccard across arms) ==")
    print(stab.sort_values(["layer", "collection", "direction", "scheme"]).to_string(index=False))


if __name__ == "__main__":
    main()
