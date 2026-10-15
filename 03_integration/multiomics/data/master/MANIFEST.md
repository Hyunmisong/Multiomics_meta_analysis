# GEC master input

Built 2026-10-14 by `scripts/00_build_master.py` from the
sources registered in `config/sources.yaml`. Harmonised, **not filtered by
significance**: linear ratios are log2-transformed, BH FDR is computed where
the source ships none, and gene symbols are resolved (UniProt accessions via
the reviewed human ID map).

| file | rows | columns |
|---|---|---|
| `GEC_master_proteome.csv` | 56,511 | dataset, feature_id, gene, logFC, pval, padj, provided_direction, calls |
| `GEC_master_phospho.csv` | 31,163 | dataset, tissue, feature_id, gene, site, raw_*, adj_*, raw_tested, adj_tested |
| `GEC_master_scrna.csv` | 873,393 | dataset, tissue, lineage, gene, logFC, pval, padj |
| `GEC_master_genome.csv` | 96,290 | dataset, tissue, gene, n_alt, n_alt_strong, mutsig_p, mutsig_q, alt_rank |

`logFC` is log2(tumour / normal) throughout; PE5 phospho writes this as
log2(T2+T3 / T0), T0 being normal mucosa and T2+T3 carcinoma.

## Provenance

| layer | cohort | tissue | patients | norm | source file | sheet | rows | genes | calls |
|---|---|---|---|---|---|---|---|---|---|
| proteome | PG1 | gastric |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PG1 | 5,458 | 5,458 | stats |
| proteome | PG2 | gastric |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PG2 | 3,904 | 3,903 | stats |
| proteome | PG3 | gastric |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PG3 | 11,340 | 11,340 | stats |
| proteome | PG4 | gastric |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PG4 | 9,745 | 9,745 | stats |
| proteome | PE1 | esophageal |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PE1 | 9,084 | 9,048 | stats |
| proteome | PE2 | esophageal |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PE2 | 9,402 | 9,325 | stats |
| proteome | PE3 | esophageal |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PE3_merged | 5,864 | 5,653 | provided |
| proteome | PE5 | esophageal |  | na | `GEC_tissue_proteomics_DEP_Master table.xlsx` | PE5 | 1,714 | 1,714 | stats |
| phospho | PG4 | gastric |  | raw | `GEC_tissue_phosphoproteome_DEP_Master_table.xlsx` | PG4 | 15,950 | 3,611 | stats |
| phospho | PG4 | gastric |  | prot_adj | `GEC_tissue_phosphoproteome_DEP_Master_table.xlsx` | PG4_adj | 14,659 | 3,155 | stats |
| phospho | PE2 | esophageal |  | raw | `PE2_phospho_adjusted_results.xlsx` | Adjusted_site_table | 11,002 | 2,848 | stats |
| phospho | PE2 | esophageal |  | prot_adj | `PE2_phospho_adjusted_results.xlsx` | Adjusted_site_table | 11,002 | 2,848 | stats |
| phospho | PE5 | esophageal |  | raw | `PE5_site_table_T2T3_vs_T0_none.csv` | - | 4,211 | 1,642 | stats |
| phospho | PE5 | esophageal |  | prot_adj | `PE5_site_table_T2T3_vs_T0_none.csv` | - | 2,340 | 726 | stats |
| genome | GG1 | gastric | 272 | na | `GG1_white_Master_Table_20260910.xlsx` | - | 19,258 | 19,258 | stats |
| genome | GG2 | gastric | 25 | na | `GG2_Master_Table_20260907.xlsx` | - | 19,258 | 19,258 | stats |
| genome | GE1 | esophageal | 508 | na | `GE1_Master_Table_20260905.xlsx` | - | 19,258 | 19,258 | stats |
| genome | GE2 | esophageal | 114 | na | `GE2_white_Master_Table_20260914.xlsx` | - | 19,258 | 19,258 | stats |
| genome | GE3 | esophageal | 89 | na | `GE3_Master_Table_20260905.xlsx` | - | 19,258 | 19,258 | stats |
| scrna | SG1 | gastric |  | na | `GEC_master_scrna.csv` | - | 69,293 | 8,692 | stats |
| scrna | SE1 | esophageal |  | na | `GEC_master_scrna.csv` | - | 272,148 | 16,012 | stats |
| scrna | SG2 | gastric |  | na | `GEC_master_scrna.csv` | - | 92,822 | 14,604 | stats |
| scrna | SM1 | mixed |  | na | `GEC_master_scrna.csv` | - | 337,103 | 31,217 | stats |
| scrna | SG3 | gastric |  | na | `GEC_master_scrna.csv` | - | 5,854 | 5,854 | stats |
| scrna | SG4 | gastric |  | na | `GEC_master_scrna.csv` | - | 8,433 | 8,433 | stats |
| scrna | SG5 | gastric |  | na | `GEC_master_scrna.csv` | - | 87,740 | 17,386 | stats |

`calls` says how significance is decided for that cohort: `stats` from the
p-values, `provided` when the source ships only a direction call (PE3).

## Source checksums (sha256, first 16 hex)

| file | sha256 |
|---|---|
| `GEC_tissue_proteomics_DEP_Master table.xlsx` | `ae87564a8c27afe3` |
| `GEC_tissue_phosphoproteome_DEP_Master_table.xlsx` | `1241706ddf120f49` |
| `PE2_phospho_adjusted_results.xlsx` | `302a1f379207fe9d` |
| `PE5_site_table_T2T3_vs_T0_none.csv` | `d8547dc3cf7faaa5` |
| `GE3_Master_Table_20260905.xlsx` | `c6841ffc28b3b1b9` |
| `GE1_Master_Table_20260905.xlsx` | `31ec5d0ab25ce042` |
| `GG2_Master_Table_20260907.xlsx` | `8bb01cade8f28c4d` |
| `GE2_white_Master_Table_20260914.xlsx` | `627c4e3bf3f6c764` |
| `GEC_master_scrna.csv` | `c6ae3fd6ab8f20ce` |
| `GG1_white_Master_Table_20260910.xlsx` | `fc9708223a66d048` |

## Dataset ids

Genome and scRNA cohorts are referred to by code, not by accession: `GG`/`GE`
for genome gastric/oesophageal, `SG`/`SE`/`SM` for scRNA gastric/oesophageal/
mixed, numbered within each group (GG1, GG2, GE1 ...). The code-to-accession
table and the raw-source registry live outside the repository, in
`config/dataset_id_map.csv` and `config/sources.yaml` on the analyst's
machine, so nothing here identifies an unpublished cohort.
