"""Consolidate the scattered cohort files into one master input.

The cohort tables arrive in three shapes -- a master workbook, a per-dataset
result workbook and a per-dataset CSV -- with different column names, fold
changes on different scales, some sheets carrying no FDR, and some carrying
direction calls instead of statistics. `config/sources.yaml` describes each of
those; this script reads them through the same adapter the pipeline uses and
writes the result as two uncut, harmonised tables:

    data/master/GEC_master_proteome.csv   one row per protein per cohort
    data/master/GEC_master_phospho.csv    one row per site per cohort,
                                          raw and protein-adjusted side by side
    data/master/GEC_master_input.xlsx     the same two, as sheets, for reading
    data/master/MANIFEST.md               provenance and row counts

Nothing is filtered by significance. From then on `config/datasets.yaml` reads
only the two CSVs, so the analysis has a single input.

    python scripts/00_build_master.py
"""
import _bootstrap  # noqa: F401

import hashlib
from datetime import date

import pandas as pd

from pmap import harmonize
from pmap.idmap import load_uniprot_map
from pmap.paths import DATA, ROOT, load_config, resolve_data_path

MASTER = DATA / "master"
PROTEOME_CSV = MASTER / "GEC_master_proteome.csv"
PHOSPHO_CSV = MASTER / "GEC_master_phospho.csv"
GENOME_CSV = MASTER / "GEC_master_genome.csv"
XLSX = MASTER / "GEC_master_input.xlsx"
MANIFEST = MASTER / "MANIFEST.md"


def _sha256(path, limit=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(limit):
            h.update(chunk)
    return h.hexdigest()[:16]


def main():
    MASTER.mkdir(parents=True, exist_ok=True)
    cfg = load_config("sources.yaml")
    uni = load_uniprot_map()

    frames, prov = [], []
    for e in cfg["datasets"]:
        if e.get("reader") == "genomics":
            df = harmonize.load_genomics(e, cfg["paths"])
        else:
            df = harmonize.load_dataset(e, cfg["paths"], uni)
        frames.append(df)
        src = resolve_data_path(cfg["paths"][e["file"]])
        prov.append({
            "layer": e["layer"], "dataset": e["name"], "norm": e.get("norm", "na"),
            "tissue": e.get("tissue", ""), "n_patients": e.get("n_patients", ""),
            "disease": e.get("disease", ""), "ethnicity": e.get("ethnicity", ""),
            "source_file": src.name, "sheet": e.get("sheet", "-"),
            "rows": len(df), "genes": df["gene"].nunique(),
            "calls": e.get("calls", "stats"),
            "fc_scale_in": e.get("fc_scale", "log2"),
        })
        print(f"  {e['layer']:9s} {e['name']:5s} {e.get('norm','na'):9s} "
              f"{len(df):7,d} rows  <- {src.name} [{e.get('sheet','-')}]")

    long = pd.concat(frames, ignore_index=True)

    # ------------------------------------------------------------ proteome
    prot = long[long["layer"] == "proteome"].copy()
    prot = prot[["dataset", "tissue", "feature_id", "gene", "logFC", "pval", "padj",
                 "provided_direction", "calls"]]
    prot.to_csv(PROTEOME_CSV, index=False)

    # ------------------------------------------------------------- phospho
    # raw and protein-adjusted go side by side on the site's row; a site tested
    # in only one arm keeps the row with NaN for the other, which is what lets
    # each arm declare its own background downstream.
    ph = long[long["layer"] == "phospho"]
    arms = []
    for norm, tag in [("raw", "raw"), ("prot_adj", "adj")]:
        a = ph[ph["norm"] == norm].set_index(["dataset", "feature_id"])
        a = a[["gene", "tissue", "site", "logFC", "pval", "padj"]].rename(
            columns={"logFC": f"{tag}_logFC", "pval": f"{tag}_P", "padj": f"{tag}_adj.P"})
        arms.append(a)
    raw_a, adj_a = arms
    phospho = raw_a.join(adj_a[[c for c in adj_a.columns if c.startswith("adj_")]],
                         how="outer")
    # gene/site come from whichever arm has them
    missing = phospho["gene"].isna()
    if missing.any():
        phospho.loc[missing, ["gene", "tissue", "site"]] = adj_a.loc[
            phospho.index[missing], ["gene", "tissue", "site"]].values
    phospho["raw_tested"] = phospho["raw_P"].notna()
    phospho["adj_tested"] = phospho["adj_P"].notna()
    phospho = phospho.reset_index()
    phospho.to_csv(PHOSPHO_CSV, index=False)

    # -------------------------------------------------------------- genome
    gen = long[long["layer"] == "genome"].copy()
    gen = gen[["dataset", "tissue", "gene", "n_alt", "n_alt_strong",
               "pval", "padj", "alt_rank"]].rename(
        columns={"pval": "mutsig_p", "padj": "mutsig_q"})
    gen[["n_alt", "n_alt_strong", "alt_rank"]] = gen[
        ["n_alt", "n_alt_strong", "alt_rank"]].astype("Int64")
    gen.to_csv(GENOME_CSV, index=False)

    with pd.ExcelWriter(XLSX) as xw:
        prot.to_excel(xw, sheet_name="Proteome", index=False)
        phospho.to_excel(xw, sheet_name="Phospho", index=False)
        gen.to_excel(xw, sheet_name="Genome", index=False)
        pd.DataFrame(prov).to_excel(xw, sheet_name="Provenance", index=False)

    # ------------------------------------------------------------ manifest
    pv = pd.DataFrame(prov)
    lines = [
        "# GEC master input",
        "",
        f"Built {date.today().isoformat()} by `scripts/00_build_master.py` from the",
        "sources registered in `config/sources.yaml`. Harmonised, **not filtered by",
        "significance**: linear ratios are log2-transformed, BH FDR is computed where",
        "the source ships none, and gene symbols are resolved (UniProt accessions via",
        "the reviewed human ID map).",
        "",
        "| file | rows | columns |",
        "|---|---|---|",
        f"| `{PROTEOME_CSV.name}` | {len(prot):,} | dataset, feature_id, gene, logFC, pval, padj, provided_direction, calls |",
        f"| `{PHOSPHO_CSV.name}` | {len(phospho):,} | dataset, tissue, feature_id, gene, site, raw_*, adj_*, raw_tested, adj_tested |",
        f"| `{GENOME_CSV.name}` | {len(gen):,} | dataset, tissue, gene, n_alt, n_alt_strong, mutsig_p, mutsig_q, alt_rank |",
        "",
        "`logFC` is log2(tumour / normal) throughout; PE5 phospho writes this as",
        "log2(T2+T3 / T0), T0 being normal mucosa and T2+T3 carcinoma.",
        "",
        "## Provenance",
        "",
        "| layer | cohort | tissue | patients | norm | source file | sheet | rows | genes | calls |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in prov:
        lines.append(f"| {r['layer']} | {r['dataset']} | {r['tissue']} | {r['n_patients']} | "
                     f"{r['norm']} | `{r['source_file']}` | {r['sheet']} | {r['rows']:,} | "
                     f"{r['genes']:,} | {r['calls']} |")
    lines += [
        "",
        "`calls` says how significance is decided for that cohort: `stats` from the",
        "p-values, `provided` when the source ships only a direction call (PE3).",
        "",
        "## Source checksums (sha256, first 16 hex)",
        "",
        "| file | sha256 |",
        "|---|---|",
    ]
    for name, path in cfg["paths"].items():
        p = resolve_data_path(path)
        lines.append(f"| `{p.name}` | `{_sha256(p) if p.exists() else 'missing'}` |")
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\nproteome  {len(prot):8,d} rows -> {PROTEOME_CSV.relative_to(ROOT)}")
    print(f"phospho   {len(phospho):8,d} rows -> {PHOSPHO_CSV.relative_to(ROOT)}")
    print(f"          raw tested {int(phospho['raw_tested'].sum()):,}, "
          f"adjusted tested {int(phospho['adj_tested'].sum()):,}")
    print(f"genome    {len(gen):8,d} rows -> {GENOME_CSV.relative_to(ROOT)}")
    print(f"workbook  -> {XLSX.relative_to(ROOT)}")
    print(f"manifest  -> {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
