# Proteomics

Differential analysis of each proteome dataset, then the cross-dataset DEP table.

```
per_dataset/   one notebook per dataset
integration/   make_dep_summary.py -> gene x dataset table of Up/Down calls
```

`make_dep_summary.py` takes the master table (one sheet per dataset) and writes a summary
workbook:

```bash
python integration/make_dep_summary.py ../data/processed/GEC_tissue_proteomics_DEP_Master_table.xlsx \
    -o ../results/proteomics/GEC_DEP_summary.xlsx --p 0.05 --lfc 1
```

This part does not share code with `02_phosphoproteomics`.
