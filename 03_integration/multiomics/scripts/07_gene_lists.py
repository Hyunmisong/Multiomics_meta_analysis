"""Gene- and pathway-level lists behind the phospho conclusion.

Who actually changes after protein adjustment, and in how many of the three
cohorts, listed so the claim can be checked rather than taken on trust.

Outputs (results/):
  phospho_genes_<norm>.csv        every gene called in that arm, with the
                                  per-cohort direction and how many cohorts
  phospho_genes_shared.csv        the same, restricted to >= 2 cohorts, both
                                  arms side by side
  pathway_cohort_counts.csv       per pathway: in how many cohorts enriched,
                                  counted at the per-cohort FDR and at the
                                  nominal p, for both arms
"""
import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

from pmap.paths import RESULTS, load_config

SCHEME = "fdr05"
COHORTS = ["PE2", "PE5", "PG4"]


def gene_list(gcalls: pd.DataFrame, norm: str) -> pd.DataFrame:
    g = gcalls[(gcalls["layer"] == "phospho") & (gcalls["norm"] == norm)
               & (gcalls["scheme"] == SCHEME)]
    wide = g.pivot_table(index="gene", columns="dataset", values="direction",
                         aggfunc="first").reindex(columns=COHORTS)
    sites = g.pivot_table(index="gene", columns="dataset", values="n_sig_features",
                          aggfunc="sum").reindex(columns=COHORTS)
    fc = g.pivot_table(index="gene", columns="dataset", values="top_absFC",
                       aggfunc="max").reindex(columns=COHORTS)
    out = pd.DataFrame(index=wide.index)
    for c in COHORTS:
        out[f"{c}_direction"] = wide[c].fillna("")
        out[f"{c}_n_sites"] = sites[c].fillna(0).astype(int)
        out[f"{c}_max_absLog2FC"] = fc[c]
    out["n_cohorts"] = (wide.notna()).sum(axis=1)
    up = (wide == "up").sum(axis=1)
    dn = (wide == "down").sum(axis=1)
    out["n_up"], out["n_down"] = up, dn
    out["consistent_direction"] = np.where(
        (up > 0) & (dn == 0), "Up",
        np.where((dn > 0) & (up == 0), "Down", "conflict"))
    out.loc[out["n_cohorts"] == 0, "consistent_direction"] = ""
    out["norm"] = norm
    return out.reset_index().sort_values(
        ["n_cohorts", "consistent_direction", "gene"], ascending=[False, True, True])


def pathway_counts(ora: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    ph = ora[(ora["layer"] == "phospho") & (ora["scheme"] == SCHEME)]
    rows = []
    for (coll, pw, norm, d), sub in ph.groupby(
            ["collection", "pathway", "norm", "direction"]):
        rows.append({
            "collection": coll, "pathway": pw, "norm": norm, "direction": d,
            "n_cohorts_tested": sub["dataset"].nunique(),
            "n_cohorts_fdr": sub.loc[sub["padj"] <= alpha, "dataset"].nunique(),
            "n_cohorts_nominal": sub.loc[sub["pval"] <= alpha, "dataset"].nunique(),
            "cohorts_fdr": ";".join(sorted(sub.loc[sub["padj"] <= alpha, "dataset"])),
            "min_padj": sub["padj"].min(),
        })
    df = pd.DataFrame(rows)
    return df.sort_values(["n_cohorts_fdr", "n_cohorts_nominal", "min_padj"],
                          ascending=[False, False, True]).reset_index(drop=True)


def main():
    gcalls = pd.read_parquet(RESULTS / "gene_calls.parquet")
    ora = pd.read_parquet(RESULTS / "ora.parquet")

    lists = {}
    for norm in ["raw", "prot_adj"]:
        gl = gene_list(gcalls, norm)
        gl.to_csv(RESULTS / f"phospho_genes_{norm}.csv", index=False)
        lists[norm] = gl
        n = gl["n_cohorts"]
        print(f"{norm:9s}: {len(gl):5,d} genes called | "
              f"3/3 cohorts {int((n == 3).sum()):4d} | 2/3 {int((n == 2).sum()):4d}")
        for k in (3, 2):
            sub = gl[gl["n_cohorts"] == k]
            cons = sub[sub["consistent_direction"].isin(["Up", "Down"])]
            print(f"             of the {k}/3: {len(cons)} with a consistent direction "
                  f"({(cons.consistent_direction == 'Up').sum()} up, "
                  f"{(cons.consistent_direction == 'Down').sum()} down)")

    shared = []
    for norm, gl in lists.items():
        s = gl[gl["n_cohorts"] >= 2].copy()
        shared.append(s)
    pd.concat(shared).to_csv(RESULTS / "phospho_genes_shared.csv", index=False)

    pc = pathway_counts(ora)
    pc.to_csv(RESULTS / "pathway_cohort_counts.csv", index=False)
    print("\npathways by cohorts reaching per-cohort FDR<=0.05:")
    print(pc.pivot_table(index=["norm", "direction"], columns="n_cohorts_fdr",
                         values="pathway", aggfunc="count").fillna(0).astype(int).to_string())

    print("\nprotein-adjusted, enriched in 2 of 3 cohorts at FDR<=0.05:")
    top = pc[(pc["norm"] == "prot_adj") & (pc["n_cohorts_fdr"] >= 2)]
    print(top[["collection", "pathway", "direction", "cohorts_fdr", "min_padj"]]
          .head(30).to_string(index=False))


if __name__ == "__main__":
    main()
