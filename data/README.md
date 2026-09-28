# Data

Source files are not tracked in git (size and redistribution). Download them into
`data/raw/<DATASET>/` with the names below.

| Dataset | Study | Source | Files |
|---|---|---|---|
| PG4 | Gastric cancer, NCC Korea | PDC000214 / PDC000215 (CDAP) | `pg4_proteomics_data.xlsx`, `pg4_phosphoproteomics_data.xlsx` (phosphosite iTRAQ table) |
| PE2 | ESCC | paper supplementary | `PE2.xlsx` |
| PE5 | ESCC, Firmiana | paper supplementary | `PE5_profile.xlsx`, `PE5_proteome_data.xlsx`, `PE5_phosphoproteome_data.xlsx` |
| PG1-PG3, PE1, PE3 | gastric / esophageal proteomes | paper supplementary | see `01_proteomics/README.md` |

Notes
- For CDAP phosphoproteome data use the **phosphosite** iTRAQ table, not the `summary` file, which
  holds spectral counts rather than ratios.
- `data/processed/` holds intermediate files written by the notebooks and is not tracked either.
