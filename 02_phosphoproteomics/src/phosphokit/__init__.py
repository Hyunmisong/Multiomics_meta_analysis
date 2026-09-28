"""Shared code for the phosphoproteomics analyses.

The per-dataset notebooks in 02_phosphoproteomics/per_dataset differ only in how the files are
read; everything after that (normalization, testing, protein adjustment, plots) is the same and
lives here.
"""
from .config import Settings, DEFAULTS
from .prep import fix_excel_dates, phospho_form, collapse_duplicates, zero_to_nan
from .normalize import normalize, quantile_normalize, normalization_bias
from .impute import impute
from .stats import (limma_test, one_sample_test, two_group_test, collapse_patients,
                    bh_adjust, direction, is_significant)
from .adjust import (adjusted_matrix, method_b, categorize, category_table,
                     CATEGORY_ORDER, CATEGORY_COLORS)
from .plots import (volcano, save_volcanoes, pvalue_histograms, sample_boxplot,
                    phospho_protein_correlation, method_comparison, embedding_plot)
from . import io

__all__ = [n for n in dir() if not n.startswith("_")]
