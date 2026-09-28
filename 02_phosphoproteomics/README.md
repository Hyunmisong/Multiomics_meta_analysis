# Phosphoproteomics

Per-dataset analysis of phosphosites together with the matching proteome, followed by the
cross-dataset integration.

```
per_dataset/   PG4.ipynb, PE2.ipynb, PE5.ipynb
integration/   meta-analysis across datasets
src/phosphokit shared code used by the notebooks
```

## What each notebook does

1. Read the dataset and build the sample table (patient, group, and the proteome sample paired with
   each phospho sample).
2. log2 transform and normalize (`quantile_common` by default; the notebook prints the bias table
   that justifies the choice).
3. Build the protein-adjusted matrix: `log2(phospho) - log2(protein)` in every sample.
4. Test protein, phospho and adjusted phospho with limma.
5. Write the site table, the unadjusted DEP list and the adjusted DEP list, plus QC figures.

## Why phosphosites are adjusted, and why both lists are kept

A phosphosite signal rises either because the site is phosphorylated more often or because there is
simply more protein. Subtracting the protein value of the same sample separates the two.

- **Unadjusted** answers "did this phosphorylation signal change?" - the total signalling output.
- **Adjusted** answers "did phosphorylation change beyond the protein?" - closer to kinase and
  phosphatase activity.

Sites whose protein was not quantified cannot be adjusted and appear only in the unadjusted list, so
reporting the adjusted list alone silently drops them.

## Method A vs Method B

| | Method A (primary) | Method B (sensitivity) |
|---|---|---|
| computation | subtract per sample, then test | test separately, then subtract |
| SE | from the differences | `sqrt(SE_p^2 + SE_prot^2)` |
| covariance | included | ignored |

`Var(X - Y) = Var(X) + Var(Y) - 2 Cov(X, Y)`. Phospho and protein measured in the same sample are
positively correlated (median r about 0.4-0.5 here), so Method B overstates the SE and is
conservative. It remains useful when the layers cannot be paired, and it uses every sample where the
protein was measured, which helps sites seen in few samples. Agreement between A and B is reported
as a robustness flag.

## Package

```python
import sys; sys.path.append("../src")
import phosphokit as pk
```

| Module | Contents |
|---|---|
| `config` | `Settings`: normalization, missing values, cutoffs |
| `io` | dataset readers and sample tables |
| `prep` | Excel-date gene symbols, peptide forms, duplicate merging |
| `normalize` | quantile / median variants and the bias check |
| `impute` | none / downshift / half_min / knn |
| `stats` | limma wrappers, one-sample and two-group tests, BH, DEP calls |
| `adjust` | adjusted matrix, Method B, site categories |
| `plots` | boxplots, volcano, p-value histograms, PCA/UMAP/t-SNE |
