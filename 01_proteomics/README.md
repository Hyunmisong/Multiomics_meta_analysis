# Proteomics

Differential analysis of each proteome dataset, then the cross-dataset DEP
table. 8 cohorts: PG1-PG4 (gastric) and PE1, PE2, PE3, PE5 (oesophageal), all
tumour vs normal.

```
per_dataset/   one notebook per dataset
integration/   make_dep_summary.py -> gene x dataset table of Up/Down calls
```

This part does not share code with `02_phosphoproteomics`.

---

## 1. The cross-dataset DEP table

`make_dep_summary.py` takes the master table (one sheet per dataset), calls
DEPs with the rule the source studies used, and writes a gene-by-dataset
summary workbook:

```bash
python integration/make_dep_summary.py ../data/processed/GEC_tissue_proteomics_DEP_Master_table.xlsx \
    -o ../results/proteomics/GEC_DEP_summary.xlsx --p 0.05 --lfc 1
```

Each cohort's sheet has its own column names and its own fold-change scale, so
the script normalises those first; PG2 ships UniProt accessions only and borrows
gene symbols from the other sheets.

### What it gave

| cohorts agreeing | genes |
|---|---|
| 7 of 8 | **1** (EPCAM) |
| 6 | 9 |
| 5 | 35 |
| 4 | 149 |
| 3 | 451 |
| 2 | 1,468 |
| 1 | 3,488 |

5,601 genes are called in at least one cohort. **One gene is called in seven,
none in all eight.** The workbook also flags 19 genes whose direction disagrees
between cohorts (`Conflicts` sheet: NDUFA13, PBLD, SERPINB5, SLC25A4, FAM3D,
ADH1C and others).

A one-gene consensus is not something to interpret. The question is why it
collapsed.

---

## 2. Why the merge collapsed: the cutoff, not the biology

The DEP rule is `raw p < 0.05 and |log2FC| > 1`. The fold-change half of that
rule means completely different things in different cohorts, because the
cohorts have completely different dynamic range:

| cohort | measured | DEPs at p<0.05 & \|log2FC\|>1 | share of measured |
|---|---|---|---|
| PG4 | 9,745 | **92** | **1%** |
| PE2 | 9,325 | 261 | 3% |
| PG3 | 11,340 | 460 | 4% |
| PG1 | 5,458 | 293 | 5% |
| PG2 | 3,903 | 553 | 14% |
| PE1 | 9,048 | 3,033 | 34% |
| PE3 | 5,653 | 2,774 | 49% |
| PE5 | 1,714 | 1,143 | **67%** |

**A 70-fold spread in how much each cohort calls.** PG4's proteome tops out at
|log2FC| 1.98, so demanding > 1 keeps almost nothing; PE5's reaches 7.4, so the
same rule keeps two thirds of what it measured. With PG4 contributing 92 genes,
no gene can reach 8 of 8 whatever the biology is.

The merge was measuring cutoff behaviour, not consensus.

---

## 3. The same cohorts at FDR < 0.05

Re-calling every cohort on its own FDR, with no fold-change filter:

| cohort | measured | p<0.05 & \|log2FC\|>1 | FDR<0.05 |
|---|---|---|---|
| PG1 | 5,458 | 293 (5%) | 763 (14%) |
| PG2 | 3,903 | 553 (14%) | 2,884 (74%) |
| PG3 | 11,340 | 460 (4%) | 4,420 (39%) |
| PG4 | 9,745 | 92 (1%) | 6,961 (71%) |
| PE1 | 9,048 | 3,033 (34%) | 4,059 (45%) |
| PE2 | 9,325 | 261 (3%) | 4,048 (43%) |
| PE3 | 5,653 | 2,774 (49%) | 2,774 (49%) |
| PE5 | 1,714 | 1,143 (67%) | 1,502 (88%) |

The spread narrows from 70-fold to 6-fold, and the consensus list stops being
empty:

| cohorts agreeing | p<0.05 & \|log2FC\|>1 | FDR<0.05 |
|---|---|---|
| 8 of 8 | 0 | **70** |
| 7 | 1 | 264 |
| 6 | 9 | 593 |
| 5 | 35 | 1,001 |
| 4 | 149 | 1,222 |

But the gene level does not become interpretable just by being larger: **21 of
those 70 eight-of-eight genes disagree on direction between cohorts.** A
consensus list that is 30% self-contradictory at its core is still not an
answer.

(PE3 is unchanged between the two rules because it ships direction calls only,
with no effect size or p-value; the pipeline reuses its shipped call under every
rule and flags it.)

---

## 4. Proteome-only pathway analysis

Moving up a level is what made the 8 cohorts speak together. Enrichment is
computed per cohort against **that cohort's own measured background**, and the
question becomes how many independent cohorts implicate the same pathway,
against a null where each cohort flags pathways at its own observed rate.

Hallmark, FDR < 0.05, up and down tested separately:

| up-regulated | cohorts | down-regulated | cohorts |
|---|---|---|---|
| E2F Targets | **8 / 8** | Oxidative Phosphorylation | 7 / 8 |
| Myc Targets V1 | **8 / 8** | Fatty Acid Metabolism | 7 / 8 |
| G2-M Checkpoint | **8 / 8** | Estrogen Response Late | 7 / 8 |
| Interferon Gamma Response | **8 / 8** | Adipogenesis | 6 / 8 |
| Myc Targets V2 | 7 / 8 | Xenobiotic Metabolism | 6 / 8 |
| Interferon Alpha Response | 7 / 8 | Myogenesis | 5 / 8 |
| DNA Repair | 6 / 8 | Cholesterol Homeostasis | 5 / 8 |
| UV Response Up | 6 / 8 | Glycolysis | 4 / 8 |

Four Hallmark sets are enriched in **every one of the eight cohorts**, against a
null expectation near one. Proliferation and interferon up, oxidative and lipid
metabolism down - a coherent picture that eight independent studies agree on,
from the same data where the gene-level consensus was a single gene.

Pathways work here for a specific reason: they pool weak per-gene evidence, and
they survive the fact that any two cohorts measure only about half a universe in
common, on a *different* half each time.

Across all three collections, proteome-only concentrated pathways at FDR < 0.05:

| rule | direction-agnostic | up | down |
|---|---|---|---|
| **FDR < 0.05** | **100** | **200** | **114** |
| FDR < 0.05 & \|log2FC\| > 0.58 | 65 | 64 | 80 |
| top/bottom 10% by log2FC | 24 | 78 | 69 |
| p < 0.05 & \|log2FC\| > 1 | 13 | 15 | 57 |

---

## 5. Which rule to use

Not a matter of taste: every rule is run on every cohort and compared on how
reproducibly its pathway calls hold across independent cohorts (mean pairwise
Jaccard of the enriched sets, 8 cohorts, 28 pairs):

| rule | Hallmark up | Hallmark down | KEGG down | Reactome up |
|---|---|---|---|---|
| **FDR < 0.05** | **0.419** | **0.295** | **0.276** | **0.272** |
| FDR < 0.05 & \|log2FC\| > 0.58 | 0.284 | 0.264 | 0.272 | 0.117 |
| top/bottom 10% by log2FC | 0.323 | 0.193 | 0.253 | 0.116 |
| p < 0.05 & \|log2FC\| > 1 | 0.213 | 0.196 | 0.171 | 0.067 |

FDR alone wins in every collection and direction; the `|log2FC| > 1` rule is
last in every one, for the dynamic-range reason in section 2. Splitting up and
down always beats the direction-agnostic call.

The original rule is kept available so results can be read against the source
studies' own DEP lists - it just should not be the setting for a meta-analysis.

---

## Where each piece lives

| | |
|---|---|
| per-dataset differential analysis | `per_dataset/` (this folder) |
| the cross-dataset DEP table (section 1) | `integration/make_dep_summary.py` (this folder) |
| sections 2-5 | [`03_integration/multiomics`](../03_integration/multiomics/README.md) |

Sections 2 to 5 are produced by the joint pipeline, which re-calls every cohort
under several rules from one uncut master input and runs the pathway analysis.
Full cutoff comparison:
[scheme_comparison.md](../03_integration/multiomics/docs/scheme_comparison.md).
The proteome results are also combined there with phosphoproteomics and
genomics, and that is where the tumour-type question - whether the gastric and
oesophageal cohorts should be analysed together at all - is answered.
