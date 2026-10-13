# Phosphoproteomics

Per-dataset analysis of phosphosites together with the matching proteome,
followed by the cross-dataset integration. Three cohorts: PG4 (gastric), PE2 and
PE5 (oesophageal), each with a paired proteome.

```
per_dataset/   PG4.ipynb, PE2.ipynb, PE5.ipynb
integration/   meta-analysis across datasets
src/phosphokit shared code used by the notebooks
```

## What each notebook does

1. Read the dataset and build the sample table (patient, group, and the proteome
   sample paired with each phospho sample).
2. log2 transform and normalize (`quantile_common` by default; the notebook
   prints the bias table that justifies the choice).
3. Build the protein-adjusted matrix: `log2(phospho) - log2(protein)` in every
   sample.
4. Test protein, phospho and adjusted phospho with limma.
5. Write the site table, the unadjusted DEP list and the adjusted DEP list, plus
   QC figures.

## Why phosphosites are adjusted, and why both lists are kept

A phosphosite signal rises either because the site is phosphorylated more often
or because there is simply more protein. Subtracting the protein value of the
same sample separates the two.

- **Unadjusted** answers "did this phosphorylation signal change?" - the total
  signalling output.
- **Adjusted** answers "did phosphorylation change beyond the protein?" - closer
  to kinase and phosphatase activity.

Sites whose protein was not quantified cannot be adjusted and appear only in the
unadjusted list, so reporting the adjusted list alone silently drops them.

## Method A vs Method B

| | Method A (primary) | Method B (sensitivity) |
|---|---|---|
| computation | subtract per sample, then test | test separately, then subtract |
| SE | from the differences | `sqrt(SE_p^2 + SE_prot^2)` |
| covariance | included | ignored |

`Var(X - Y) = Var(X) + Var(Y) - 2 Cov(X, Y)`. Phospho and protein measured in
the same sample are positively correlated (median r about 0.4-0.5 here), so
Method B overstates the SE and is conservative. It remains useful when the
layers cannot be paired, and it uses every sample where the protein was
measured, which helps sites seen in few samples. Agreement between A and B is
reported as a robustness flag.

---

## What the two lists actually do to the answer

The two lists were built because they answer different questions. Putting the
three cohorts side by side shows how far apart those answers are.

### Site level

Restricted to sites tested in **both** arms, called at FDR < 0.05:

| | PE2 | PE5 | PG4 |
|---|---|---|---|
| paired sites | 11,002 | 2,340 | 14,659 |
| significant, unadjusted | 4,714 | 755 | 5,742 |
| significant, adjusted | 2,688 | 621 | 3,914 |
| unadjusted only | 2,719 | 441 | 2,943 |
| adjusted only | 693 | 307 | 1,115 |
| **direction flips** | **126** | **43** | **77** |
| median \|log2FC\| unadjusted -> adjusted | 0.290 -> 0.229 | 0.968 -> 1.234 | 0.177 -> 0.132 |
| Spearman of log2FC | 0.65 | 0.44 | 0.69 |

Adjustment removes 18-43% of the significant sites and adds a smaller,
different set back. A few dozen sites per cohort change direction outright.
Effect sizes shrink in PE2 and PG4 and grow in PE5, which also has the weakest
agreement between arms.

### Pathway level

Median Jaccard of the enriched pathway sets, unadjusted vs adjusted:

| cutoff | PE2 | PE5 | PG4 |
|---|---|---|---|
| FDR < 0.05 | 0.09 | 0.00 | 0.08 |
| FDR < 0.05 & \|log2FC\| > 0.58 | 0.28 | 0.00 | 0.17 |
| p < 0.05 & \|log2FC\| > 1 | 0.11 | 0.00 | 0.13 |
| top/bottom 10% by log2FC | 0.00 | 0.00 | 0.11 |

Under no cutoff are these two views of one result.

### The difference is directional

Pathways concentrated across all three cohorts at FDR < 0.05:

| | direction-agnostic | up | down |
|---|---|---|---|
| unadjusted | 10 | **20** | **8** |
| adjusted | 13 | **0** | **0** |

Both arms concentrate - the adjusted arm slightly more when direction is
ignored (Mitotic Spindle and Apical Junction in 3/3 cohorts). What disappears
is agreement on **which way** a pathway moves. The unadjusted arm's directional
signal is Myc Targets V1 and V2 up in 3/3 and Myogenesis down in 3/3 - the same
proliferation signature the proteome layer produces on its own, i.e. a
protein-abundance echo.

### What survives adjustment

At two of three cohorts and per-cohort FDR < 0.05 the adjusted arm gives a
coherent signature, always in the same two cohorts (PE2 and PG4): **24 of 25
pathways down** - RHO and RAC1 GTPase cycles, regulation of actin cytoskeleton,
focal adhesion, tight junction, Netrin-1 signalling, and signalling by
non-receptor tyrosine kinases, PTK6, FGFR1/3/4 and ErbB.

Gene level, adjusted arm, called in all three cohorts: **174 genes**, of which
98 agree on direction (31 up, 67 down) and 76 conflict.

**Practical reading:** the unadjusted list should not be reported as signalling
evidence on its own - its directional pathway signal is largely the proteome's.
The adjusted list is not empty either; it carries a smaller, cytoskeleton- and
tyrosine-kinase-centric signature that reproduces across the two cohorts deep
enough to show it.

---

## What the comparison needs from the output format

The comparison above was initially impossible, and the reason is worth keeping
in mind when writing results out.

The notebooks produce three things per contrast: the **full site table** (every
tested site, both arms on the same row) and two **DEP lists** (the sites
significant in one arm or the other). A DEP list cannot support this
comparison:

- a paired raw-vs-adjusted count conditions on significance if either side is
  pre-filtered - everything in the adjusted DEP list is adjusted-significant by
  construction;
- an enrichment test needs a background of everything measured, which a DEP list
  does not carry.

Only the full site table answers both. When PE2 and PE5 first entered the joint
analysis as DEP lists, their paired statistics and enrichment backgrounds were
unusable; swapping in `Adjusted_site_table` (PE2) and the full
`site_table_<contrast>` (PE5) fixed it, and all three cohorts became comparable.

**Reproduction check.** Run under this part's own DEP rule (raw p < 0.05 and
|log2FC| > 1), the joint pipeline reproduces these notebooks' DEP counts
exactly: PE5 1,534 unadjusted and 903 adjusted, PE2 359 adjusted. (The paired
table above shows smaller unadjusted counts because it is restricted to sites
tested in both arms.)

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

`rbridge.py` converts numpy inputs to R explicitly and calls limma **outside**
the converter context, rather than wrapping the call in `with NP2R.context()`.
That is deliberate: the context converts the return value too, and the Python
type rpy2 gives an R named list has changed across versions (`ListVector`, then
`NamedList` whose `.items()` yields `NamedItem` objects). Converting only the
inputs leaves a plain R object that `.rx2()` can address on every rpy2 version.
A standalone notebook that inlined the older pattern broke on an rpy2 upgrade
with `cannot unpack non-iterable NamedItem object`.

---

## Where each piece lives

| | |
|---|---|
| per-dataset analysis, both lists, QC | `per_dataset/` and `src/phosphokit` (this folder) |
| the cross-cohort comparison above | [`03_integration/multiomics`](../03_integration/multiomics/README.md) |

The comparison is run by the joint pipeline, which reads the full site tables
from this part, re-calls every cohort under several cutoffs from one uncut
master input, and combines the result with proteomics and genomics. Full cutoff
comparison:
[scheme_comparison.md](../03_integration/multiomics/docs/scheme_comparison.md).
