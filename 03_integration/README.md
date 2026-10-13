# Integration

Combines the proteomics and phosphoproteomics results with genomics and
scRNA-seq data.

```
genomics/     processing of the collaborator's variant / CNV data
scrna/        scRNA-seq processing and cell-type signatures
multiomics/   joint analysis across layers
```

Inputs come from `results/` of the other two parts plus the external data
described in [`data/README.md`](../data/README.md).

## Status

| | state |
|---|---|
| **[multiomics](multiomics/README.md)** | **the joint pathway meta-analysis - 8 proteome, 3 phospho x 2 normalisations, 5 genome cohorts. Findings, statistics and figures are there.** |
| genomics | the 5 variant cohorts are read directly by `multiomics` (truncating LOF only; no copy number yet). This folder is for the collaborator's CNV data when it arrives. |
| scrna | not started. The joint pipeline already reserves a cell-type axis and counts cohorts rather than arms, so the layer slots in without changing the analysis. |

## What the joint analysis found

Details, figures and limitations in
**[multiomics/README.md](multiomics/README.md)**.

- **Phospho quantification changes the pathway answer.** Normalising
  phosphosites to their parent protein removes 18-43% of the significant sites,
  and the enriched pathway sets barely overlap under any cutoff. The unadjusted
  arm's directional signal is largely an echo of protein abundance.
- **Gastric and oesophageal cohorts should be analysed apart.** Tumour type
  explains 34% of the variation between proteome cohorts. The member genes show
  why: one pathway can contain genes moving in opposite directions in the two
  types, so pooling cancels rather than dilutes.
- **The genome layer covers truncating loss of function only** - no copy
  number, missense or structural variants - so an absent genomics signal means
  "not recurrently truncated", never "not altered".

Generated reports, regenerated on every run:
[scheme comparison](multiomics/docs/scheme_comparison.md),
[gastric vs oesophageal](multiomics/docs/tissue_split.md),
[master input manifest](multiomics/data/master/MANIFEST.md).
