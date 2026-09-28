# Multiomics_meta_analysis

Meta-analysis of published proteomics and phosphoproteomics datasets of gastric and esophageal
cancer, and integration with genomics and scRNA-seq data.

## Layout

```
01_proteomics/          differential analysis per dataset, then the cross-dataset DEP table
02_phosphoproteomics/   phospho vs protein per dataset (shared code in src/phosphokit), then integration
03_integration/         genomics, scRNA-seq and the multi-omics layer
data/                   inputs (not tracked; see data/README.md for sources)
results/                tables and figures written by the notebooks
docs/                   methods and decisions
```

`01_proteomics` and `02_phosphoproteomics` are independent: the shared package lives inside
`02_phosphoproteomics/src` and is used only by the phospho notebooks.

## Setup

```bash
conda env create -f environment.yml
conda activate multiomics
```

The phosphoproteomics notebooks call limma through `rpy2`, so R with `limma` and `statmod` is
required; both are in `environment.yml`.

## Running

1. Download the source files into `data/raw/<DATASET>/` following `data/README.md`.
2. Run the notebook of the dataset from its own folder, e.g.
   `02_phosphoproteomics/per_dataset/PE5.ipynb`. Paths in the notebooks are relative to that folder.
3. Results are written to `results/<DATASET>/`.

## Datasets

| ID | Tissue | Layers | Comparison |
|---|---|---|---|
| PG1-PG4 | gastric | proteome (PG4 also phospho) | tumor vs normal |
| PE1-PE5 | esophageal | proteome (PE2, PE5 also phospho) | tumor vs normal, PE5 by tissue stage |

## Method summary

- Differential analysis with limma; `duplicateCorrelation(block = patient)` whenever a patient
  contributes to both arms of a comparison.
- Normalization restricted to features detected in most samples, because tumor samples detect more
  low-abundance features and aligning all observed values shifts them upward.
- Phosphosites are reported twice: unadjusted, and adjusted for the abundance of their protein
  (per-sample subtraction). The two answer different questions, and sites whose protein was not
  quantified exist only in the unadjusted list.

Details are in `docs/methods.md`.
