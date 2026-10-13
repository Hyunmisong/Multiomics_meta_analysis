"""Load every registered dataset into one standard long table.

Standard schema (one row per measured feature per dataset arm):

    layer       proteome | phospho | transcriptome | genome | scrna
    dataset     cohort id (PG1 ... PE5)
    norm        raw | prot_adj | na
    arm         "{layer}|{dataset}|{norm}" -- the unit everything else keys on
    feature_id  native id (protein, phosphosite, peptide, ENSG)
    gene        HGNC symbol, upper case
    site        residue/peptide for phospho, else ''
    logFC       log2 tumour / normal (linear ratios are log2-transformed)
    pval, padj  BH-adjusted p is computed when the source ships none
    calls       stats | provided | prefiltered  (how significance is decided)
    provided_direction  Up/Down call shipped by the source, when there is one

Threshold schemes are applied downstream (`apply_scheme`), not here, so one
load serves all three.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .idmap import clean_symbol, load_uniprot_map, strip_isoform
from .paths import load_config, resolve_data_path

STD_COLS = [
    "layer", "dataset", "norm", "tissue", "arm", "feature_id", "gene", "site",
    "logFC", "pval", "padj", "calls", "provided_direction",
    # genome layer only: truncating-variant counts and the MutSig rank. NaN
    # elsewhere, so every layer still shares one schema.
    "n_alt", "n_alt_strong", "alt_rank",
]

# Sequence Ontology terms that are unambiguously loss of function. A row is
# kept when any of its terms (the sources join them with ";" or "&") is one of
# these, so "frameshift_variant&splice_region_variant" counts and a bare
# "splice_region_variant" would not.
LOF_TERMS = {"frameshift_variant", "stop_gained"}


def arm_id(layer: str, dataset: str, norm: str) -> str:
    return f"{layer}|{dataset}|{norm}"


def _get(df: pd.DataFrame, col, default=np.nan):
    if col is None or col not in df.columns:
        return pd.Series(default, index=df.index)
    return df[col]


def load_dataset(entry: dict, paths: dict, uniprot_map: dict) -> pd.DataFrame:
    src = resolve_data_path(paths[entry["file"]])
    if src.suffix.lower() in {".csv", ".tsv", ".txt"}:
        sep = "	" if src.suffix.lower() in {".tsv", ".txt"} else ","
        raw = pd.read_csv(src, sep=sep, low_memory=False)
    else:
        raw = pd.read_excel(src, sheet_name=entry["sheet"], header=entry.get("header", 0))
    raw = _apply_where(raw, entry.get("where")).reset_index(drop=True)

    out = pd.DataFrame(index=raw.index)
    out["layer"] = entry["layer"]
    out["dataset"] = entry["name"]
    out["norm"] = entry.get("norm", "na")
    out["arm"] = arm_id(entry["layer"], entry["name"], entry.get("norm", "na"))
    out["feature_id"] = _get(raw, entry.get("feature_id")).astype(str)
    out["site"] = _get(raw, entry.get("site"), "").astype(str).replace("nan", "")

    # --- gene symbol ------------------------------------------------------
    gene = _get(raw, entry.get("gene"), "").map(clean_symbol)
    acc_col = entry.get("uniprot")
    if (gene == "").any() and acc_col and acc_col in raw.columns:
        acc = raw[acc_col].map(strip_isoform)
        mapped = acc.map(lambda a: uniprot_map.get(a, ""))
        gene = gene.where(gene != "", mapped)
    out["gene"] = gene

    # --- effect size ------------------------------------------------------
    fc = pd.to_numeric(_get(raw, entry.get("logFC")), errors="coerce")
    if entry.get("fc_scale") == "linear":
        fc = np.log2(fc.where(fc > 0))
    out["logFC"] = fc

    # --- p-values ---------------------------------------------------------
    out["pval"] = pd.to_numeric(_get(raw, entry.get("pval")), errors="coerce")
    padj = pd.to_numeric(_get(raw, entry.get("padj")), errors="coerce")
    if padj.isna().all() and out["pval"].notna().any():
        padj = pd.Series(np.nan, index=raw.index)
        ok = out["pval"].notna()
        padj.loc[ok] = multipletests(out.loc[ok, "pval"], method="fdr_bh")[1]
    out["padj"] = padj

    out["tissue"] = entry.get("tissue", "")
    for c in ("n_alt", "n_alt_strong", "alt_rank"):
        out[c] = pd.to_numeric(_get(raw, entry.get(c)), errors="coerce")
    out["calls"] = entry.get("calls", "stats")
    pd_col = entry.get("direction")
    if pd_col and pd_col in raw.columns:
        out["provided_direction"] = raw[pd_col].astype(str).str.strip().str.title()
        out.loc[~out["provided_direction"].isin(["Up", "Down"]), "provided_direction"] = ""
    else:
        out["provided_direction"] = ""

    out = out[out["gene"] != ""].copy()

    # A table may carry rows that this arm never actually tested -- a full site
    # table lists every feature, but a feature can fail the min-patients filter
    # in one arm and not the other. `require_notna` names the column that must
    # be present for the row to belong to this arm, which keeps each arm's
    # background to what it really measured.
    req = entry.get("require_notna")
    if req:
        out = out[pd.to_numeric(raw.loc[out.index, req], errors="coerce").notna()].copy()
    return out[STD_COLS]


def _apply_where(raw: pd.DataFrame, where: dict | None) -> pd.DataFrame:
    """Row filter, so one consolidated master table can serve several arms."""
    if not where:
        return raw
    mask = pd.Series(True, index=raw.index)
    for col, val in where.items():
        mask &= raw[col].astype(str) == str(val)
    return raw[mask]


def load_genomics(entry: dict, paths: dict) -> pd.DataFrame:
    """Gene-level evidence for one genomics cohort.

    Two sheets are combined: the truncating-variant list gives how often a gene
    is hit and how badly, and MutSig2CV gives whether that is more than the
    background mutation rate explains. The MutSig gene list is also the
    cohort's universe, so genes with no variant stay in the table with a count
    of zero -- an enrichment background needs them.
    """
    src = resolve_data_path(paths[entry["file"]])
    var = pd.read_excel(src, sheet_name=str(entry.get("variant_sheet", "2")))
    sig = pd.read_excel(src, sheet_name=str(entry.get("mutsig_sheet", "3")))

    var = var.rename(columns={c: c.replace("positon", "position") for c in var.columns})
    terms = var["Sequence ontology"].fillna("").map(
        lambda v: {t.strip() for t in re.split(r"[;&]", str(v)) if t.strip()})
    var = var[terms.map(lambda t: bool(t & LOF_TERMS))].copy()
    var["gene"] = var["Gene symbol"].map(clean_symbol)
    var = var[var["gene"] != ""]
    # fraction of the protein lost downstream of the variant
    ratio = pd.to_numeric(var["Protein truncation ratio"], errors="coerce")
    agg = var.groupby("gene").agg(n_alt=("gene", "size"))
    agg["n_alt_strong"] = var.assign(strong=ratio >= 50).groupby("gene")["strong"].sum()

    sig = sig.rename(columns={c: c.lower() for c in sig.columns})
    sig["gene"] = sig["gene"].map(clean_symbol)
    sig = sig[sig["gene"] != ""].drop_duplicates("gene").set_index("gene")
    sig["alt_rank"] = sig["p"].rank(method="min")

    out = pd.DataFrame(index=sig.index)
    out["layer"] = "genome"
    out["dataset"] = entry["name"]
    out["norm"] = entry.get("norm", "na")
    out["tissue"] = entry.get("tissue", "")
    out["arm"] = arm_id("genome", entry["name"], entry.get("norm", "na"))
    out["feature_id"] = sig.index
    out["gene"] = sig.index
    out["site"] = ""
    out["logFC"] = np.nan                      # LOF has no effect size
    out["pval"] = pd.to_numeric(sig["p"], errors="coerce")
    out["padj"] = pd.to_numeric(sig["q"], errors="coerce")
    out["calls"] = "genome"
    out["provided_direction"] = ""
    out["n_alt"] = agg["n_alt"].reindex(sig.index).fillna(0)
    out["n_alt_strong"] = agg["n_alt_strong"].reindex(sig.index).fillna(0)
    out["alt_rank"] = sig["alt_rank"]
    return out.reset_index(drop=True)[STD_COLS]


def load_all(config: str = "datasets.yaml", force_idmap: bool = False) -> pd.DataFrame:
    cfg = load_config(config)
    uniprot_map = load_uniprot_map(force=force_idmap)
    frames = [load_genomics(e, cfg["paths"]) if e.get("reader") == "genomics"
              else load_dataset(e, cfg["paths"], uniprot_map)
              for e in cfg["datasets"]]
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------- universes
def universes(long: pd.DataFrame, config: str = "datasets.yaml") -> dict:
    """{arm: set(genes)} honouring the `universe:` key (borrowed backgrounds).

    An arm whose table is a pre-filtered subset cannot supply its own
    background, so it borrows the background of the arm named in the config.
    """
    cfg = load_config(config)
    own = {arm: set(g) for arm, g in long.groupby("arm")["gene"]}
    uni = {}
    for e in cfg["datasets"]:
        arm = arm_id(e["layer"], e["name"], e.get("norm", "na"))
        if arm not in own:
            continue
        spec = e.get("universe", "self")
        if isinstance(spec, str) and spec.startswith("ref:"):
            ref = spec[4:]
            uni[arm] = own.get(ref, own[arm]) | own[arm]
        else:
            uni[arm] = own[arm]
    return uni


# ------------------------------------------------------------------ schemes
def apply_scheme(long: pd.DataFrame, scheme_name: str, scheme: dict) -> pd.DataFrame:
    """Add is_sig / direction columns for one threshold scheme.

    Arms with calls == 'provided'    use their shipped Up/Down call verbatim.
    Arms with calls == 'prefiltered' treat every shipped row as significant and
                                     take the direction from the sign of logFC.
    Every other arm is thresholded by the scheme.
    """
    df = long.copy()
    df["scheme"] = scheme_name
    df["is_sig"] = False
    df["direction"] = ""

    # The genome layer is scored by its own schemes (see apply_genome_scheme);
    # under a continuous-layer scheme its rows are simply not called.
    gen = df["calls"] == "genome"
    df.loc[gen, "is_sig"] = False

    prov = df["calls"] == "provided"
    df.loc[prov, "is_sig"] = df.loc[prov, "provided_direction"].isin(["Up", "Down"])
    df.loc[prov, "direction"] = df.loc[prov, "provided_direction"].str.lower()

    pre = df["calls"] == "prefiltered"
    df.loc[pre, "is_sig"] = True
    df.loc[pre, "direction"] = np.where(df.loc[pre, "logFC"] >= 0, "up", "down")

    stat = df["calls"] == "stats"
    if scheme.get("mode") == "rank":
        frac = scheme["tail_frac"]
        for arm, idx in df[stat].groupby("arm").groups.items():
            fc = df.loc[idx, "logFC"]
            if fc.notna().sum() < 20:
                continue
            lo, hi = fc.quantile(frac), fc.quantile(1 - frac)
            df.loc[idx, "is_sig"] = (fc <= lo) | (fc >= hi)
            df.loc[idx, "direction"] = np.where(fc >= hi, "up", np.where(fc <= lo, "down", ""))
    else:
        if "padj_max" in scheme:
            sig = stat & (df["padj"] <= scheme["padj_max"])
        else:
            sig = stat & (df["pval"] <= scheme["pval_max"])
        if scheme.get("abs_logfc_min", 0) > 0:
            sig &= df["logFC"].abs() >= scheme["abs_logfc_min"]
        df.loc[sig, "is_sig"] = True
        df.loc[sig, "direction"] = np.where(df.loc[sig, "logFC"] >= 0, "up", "down")

    df.loc[~df["is_sig"], "direction"] = ""
    return df


def gene_level(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse features to genes within an arm.

    Phospho arms have many sites per gene; a gene counts as regulated in that
    arm if any site is significant. When sites disagree the gene is called
    'mixed' and is counted in neither the up nor the down set, but still counts
    as a hit under direction-agnostic (`any`) enrichment.
    """
    sig = df[df["is_sig"]]
    g = sig.groupby(["arm", "layer", "dataset", "norm", "tissue", "scheme", "gene"]).agg(
        n_sig_features=("direction", "size"),
        n_up=("direction", lambda s: int((s == "up").sum())),
        n_down=("direction", lambda s: int((s == "down").sum())),
        top_absFC=("logFC", lambda s: s.abs().max()),
        mean_logFC=("logFC", "mean"),
    ).reset_index()
    g["direction"] = np.where(g["n_down"] == 0, "up",
                       np.where(g["n_up"] == 0, "down", "mixed"))
    return g


def apply_genome_scheme(long: pd.DataFrame, scheme_name: str, scheme: dict) -> pd.DataFrame:
    """Call genome-layer genes under one genome scheme.

    A truncating variant is loss of function, so there is only one direction:
    every called gene is "down". That is what makes the layer combine cleanly
    with the others -- genome down + proteome down + phospho down is three
    layers agreeing, not three layers that happen to be significant.
    """
    df = long[long["calls"] == "genome"].copy()
    df["scheme"] = scheme_name
    ok = pd.Series(True, index=df.index)
    if "q_max" in scheme:
        ok &= df["padj"] <= scheme["q_max"]
    if "max_rank" in scheme:
        ok &= df["alt_rank"] <= scheme["max_rank"]
    ok &= df["n_alt"] >= scheme.get("min_variants", 1)
    if scheme.get("min_strong"):
        ok &= df["n_alt_strong"] >= scheme["min_strong"]
    df["is_sig"] = ok.fillna(False)
    df["direction"] = np.where(df["is_sig"], "down", "")
    return df
