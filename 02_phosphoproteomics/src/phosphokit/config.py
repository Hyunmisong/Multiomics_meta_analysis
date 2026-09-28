"""Analysis settings shared by all datasets."""
from dataclasses import dataclass


@dataclass
class Settings:
    # --- normalization ---
    # Tumor samples usually detect many more (mostly low-abundance) features than normal samples.
    # Aligning ALL observed values (per-sample median or plain quantile normalization) therefore
    # shifts tumor samples upward: their observed distribution contains the extra low-intensity
    # features. Restricting the normalization to features seen in most samples compares like with
    # like. In PE5 this moved the median protein log2FC of T2 vs T0 from about +1.5 to +0.4.
    normalization: str = "quantile_common"
    norm_common_frac: float = 0.7

    # --- missing values ---
    # "none" keeps NaN and tests each feature on the samples where it was observed. Imputation
    # makes peptides that are detected almost only in tumors testable, at the price of inventing
    # part of the signal.
    missing_method: str = "none"
    impute_min_valid_frac: float = 0.5
    downshift: float = 1.8
    downshift_width: float = 0.3
    knn_k: int = 10

    # --- testing ---
    min_per_group: int = 3      # minimum samples (or patients) with a value in each arm
    p_cut: float = 0.05
    lfc_cut: float = 1.0
    use_adj_p: bool = False     # True -> BH-adjusted p-values decide the DEP call

    seed: int = 0


DEFAULTS = Settings()
