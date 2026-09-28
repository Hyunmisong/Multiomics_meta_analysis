# Decisions log

- **Per-sample subtraction over MSstatsPTM-style combination.** Both are computed; A is primary
  because proteome and phosphoproteome come from the same lysate here, so the covariance term is
  real and ignoring it inflates the SE (the SE ratio grows with the phospho-protein correlation,
  which has median r about 0.4-0.5).
- **Common-feature normalization.** Median and plain quantile normalization left a median protein
  log2FC of about +1.5 for PE5 T2 vs T0, driven by detection depth, not biology.
- **Technical replicates averaged (PG4).** 50 of 80 patients were measured twice; treating the runs
  as separate patients would inflate n. Replicate correlation is about 0.94-0.97.
- **Patient-level averaging within a group (PE5).** T1 patients contribute 1-6 lesions.
- **Sequence windows for cross-dataset matching.** RefSeq positions shift between isoforms and
  database versions; +/-7 aa windows are more stable.
- **Both DEP lists reported.** Sites whose protein was not quantified exist only in the unadjusted
  list (644 of them in PE5 T2+T3 vs T0).
- **rpy2 output handled by name.** The automatic conversion of R lists changed between rpy2
  versions; `rbridge.r_list_to_dict` reads the list by element name so both behave the same.
