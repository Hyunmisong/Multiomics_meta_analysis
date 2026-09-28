# Methods

## Normalization
Log2 values are normalized per sample. Because tumor samples detect far more (mostly low-abundance)
features than normal samples, aligning all observed values (per-sample median, plain quantile
normalization) shifts tumor samples upward. The default `quantile_common` therefore quantile-
normalizes the features detected in at least 70% of samples and applies the same per-sample offset
to the rest. Each notebook prints the median log2FC of the common features per method so the choice
can be checked rather than assumed.

## Differential analysis
limma, with the design matching the data:
- ratio tables (one log2 tumor/normal value per patient): intercept-only model;
- intensity tables with separate samples: two-group model.

Samples of one patient within one group are averaged first, since several lesions of a patient are
not independent replicates of that group. Patients present in both arms are handled by
`duplicateCorrelation(block = patient)`, which keeps the patients that appear in only one arm.
Features are filtered to a minimum number of values per arm before `eBayes`, so the variance prior
is learned from the features that are actually reported.

DEP call: p < 0.05 and |log2FC| > 1 by default (BH-adjusted p optionally).

## Protein adjustment
Method A subtracts the protein log2 value of the same sample from the phosphosite value and tests
the differences; the within-sample covariance is therefore included. Method B tests the layers
separately and combines the standard errors as if independent (MSstatsPTM-style), which is
conservative under positive correlation. Method A is the primary result and Method B is reported as
a sensitivity check.

## Site categories
Per comparison, each adjustable site is labelled phospho-specific, both (same or opposite
direction), explained by protein, protein only, or not significant.

## Missing values
Default is no imputation: each feature is tested on the samples where it was observed. Imputation
options (downshift, half-min, kNN) make features that are detected almost only in tumors testable,
at the price of inventing part of the signal.
