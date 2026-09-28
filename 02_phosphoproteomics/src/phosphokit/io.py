"""Dataset-specific readers.

Each function returns tidy objects for the rest of the package: a numeric matrix
(features x samples) plus an annotation frame, and where applicable a sample table with
Case / Group columns.
"""
import re

import numpy as np
import pandas as pd

from .prep import fix_excel_dates, numeric

SAMPLE_RE_PE2 = re.compile(r"^ESCC-\d+[NT]$")


# ---------------------------------------------------------------- CDAP ratio tables (PG datasets)
def read_cdap_ratios(path, index_col, value_prefix="Log ", sheet=0, annot_cols=()):
    """
    CPTAC CDAP tables: one column per TMT/iTRAQ channel pair, already log2 tumor/normal.

    `value_prefix` selects "Log " (all peptides) or "Unshared Log " (peptides unique to the
    protein, the usual choice for protein-level quantification). Column names such as
    "Log T112.1/N111.1" carry the patient and the technical replicate.
    """
    df = pd.read_excel(path, sheet_name=sheet) if str(path).endswith((".xlsx", ".xls")) \
        else pd.read_csv(path, sep="\t", low_memory=False)
    df = df[~df[index_col].isin(["Mean", "Median", "StdDev"])]        # CDAP summary rows
    value_cols = [c for c in df.columns if str(c).startswith(value_prefix)]
    mat = numeric(df.set_index(index_col)[value_cols])
    mat.columns = [c.replace(value_prefix, "") for c in value_cols]
    annot = df.set_index(index_col)[[c for c in annot_cols if c in df.columns]]
    return mat, annot


def cdap_sample_table(columns):
    """'T112.1/N111.1' -> patient T112/N111, replicate 1."""
    out = pd.DataFrame({"sample": list(columns)})
    out["Case"] = out["sample"].str.replace(r"\.\d+", "", regex=True)
    out["replicate"] = out["sample"].str.extract(r"\.(\d+)/")[0].fillna("1")
    out["Group"] = "ratio"
    return out.set_index("sample")


# ---------------------------------------------------------------- PE2 (ESCC, per-sample values)
def read_pe2_sheet(path, sheet):
    """Row 0 is a group banner and row 1 the header; sample columns look like ESCC-001N / ESCC-001T."""
    df = pd.read_excel(path, sheet_name=sheet, header=1)
    sample_cols = [c for c in df.columns if SAMPLE_RE_PE2.match(str(c))]
    annot = df[[c for c in df.columns if c not in sample_cols]].copy()
    if "Gene name" in annot:
        annot["Gene name"] = annot["Gene name"].map(fix_excel_dates)
    return numeric(df[sample_cols]), annot


def pe2_sample_table(columns):
    out = pd.DataFrame({"sample": list(columns)})
    out["Case"] = out["sample"].str[:-1]
    out["Group"] = out["sample"].str[-1]          # N or T
    return out.set_index("sample")


# ---------------------------------------------------------------- PE5 (ESCC, Firmiana profile)
def norm_experiment_id(x):
    """Firmiana experiment IDs are 6-digit codes; Excel drops their leading zeros."""
    return re.sub(r"\.0$", "", str(x).strip()).zfill(6)


def read_pe5_profile(path, sheet, case_col, first_data_col, layer, substage_to_group):
    """
    Two rows per case: experiment IDs on the first, substage labels on the second.
    Substages use '-' in the proteome sheet and '_' in the phospho sheet.
    """
    df = pd.read_excel(path, sheet_name=sheet, dtype=str)
    data_cols = df.columns[df.columns.get_loc(first_data_col):]
    ids, subs = df.iloc[0::2].reset_index(drop=True), df.iloc[1::2].reset_index(drop=True)
    rows = []
    for i in range(len(ids)):
        for col in data_cols:
            exp_id, sub = ids.at[i, col], subs.at[i, col]
            if pd.isna(exp_id) or pd.isna(sub):
                continue
            rows.append({"Case": str(ids.at[i, case_col]).strip(),
                         "Substage": re.sub(r"[_\s]", "-", str(sub).strip()),
                         f"ExpID_{layer}": norm_experiment_id(exp_id)})
    out = pd.DataFrame(rows)
    out["Group"] = out["Substage"].map(substage_to_group)
    return out


def pair_by_substage(phospho_profile, proteome_profile):
    """
    Match the two layers on (Case, Substage), not on Case alone: one case contributes several
    substages, and only the same tissue piece may be paired.
    """
    for name, df in [("phospho", phospho_profile), ("proteome", proteome_profile)]:
        assert not df.duplicated(["Case", "Substage"]).any(), f"{name}: duplicated (Case, Substage)"
    merged = phospho_profile.merge(proteome_profile[["Case", "Substage", "ExpID_proteome"]],
                                   on=["Case", "Substage"], how="left", validate="one_to_one")
    merged["paired"] = merged["ExpID_proteome"].notna()
    return merged.set_index("ExpID_phospho")
