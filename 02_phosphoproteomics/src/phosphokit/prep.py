"""Cleaning steps that recur across datasets."""
import re

import numpy as np
import pandas as pd

# Excel silently converts some gene symbols to dates (SEPT9 -> 2021-09-09, MARCH5 -> 2021-03-05).
# The old HGNC spellings are used here because both proteome and phospho tables are affected the
# same way, so the two layers still match each other.
_MONTH_PREFIX = {3: "MARCH", 9: "SEPT", 10: "OCT", 12: "DEC"}


def fix_excel_dates(symbol):
    if pd.isna(symbol):
        return symbol
    m = re.match(r"^\d{4}-(\d{2})-(\d{2})", str(symbol))
    if not m:
        return str(symbol)
    prefix = _MONTH_PREFIX.get(int(m.group(1)))
    return f"{prefix}{int(m.group(2))}" if prefix else str(symbol)


def zero_to_nan(df, placeholder=0.0):
    """A 0 (or a placeholder such as 1e-5) means 'not detected', not 'measured as zero'."""
    return df.mask(df <= placeholder)


def phospho_form(seq):
    """
    Upper-case every modified residue except phosphorylated S/T/Y.

    Sequences mark modifications in lower case, so the same phosphopeptide appears several times
    when it also carries an oxidised Met ('m') or an N-terminal acetylation ('a'). Those forms are
    the same phosphorylation event and are merged.
    """
    return "".join(ch if ch in "sty" else ch.upper() for ch in str(seq))


def count_phosphosites(seq):
    return sum(ch in "sty" for ch in str(seq))


def collapse_duplicates(df, feature, value_cols, rule="sum"):
    """
    Merge rows that describe the same feature.

    "sum"     : add the intensities (min_count=1 keeps all-missing rows as NaN rather than 0)
    "max_obs" : keep the row with the most non-missing values
    """
    if rule == "sum":
        return df.groupby(feature)[value_cols].sum(min_count=1)
    if rule == "max_obs":
        order = df[value_cols].notna().sum(axis=1)
        return (df.assign(_n=order).sort_values("_n", ascending=False)
                .drop_duplicates(feature).set_index(feature)[value_cols])
    raise ValueError(rule)


def sequence_window(seq, flank=7):
    """
    +/- `flank` residues around each phosphosite, used to match sites across datasets.
    Sequence windows travel better than RefSeq position numbers, which shift between
    isoforms and database versions.
    """
    plain = str(seq).upper()
    out = []
    for i, ch in enumerate(str(seq)):
        if ch in "sty":
            start, end = max(0, i - flank), min(len(plain), i + flank + 1)
            out.append(plain[start:end].rjust(flank + 1 + (i - start == flank) * flank, "_"))
    return ";".join(out)


def numeric(df):
    return df.apply(pd.to_numeric, errors="coerce").astype(float).replace([np.inf, -np.inf], np.nan)
