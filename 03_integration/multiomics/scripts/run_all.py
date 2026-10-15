"""Run the whole prototype end to end."""
import _bootstrap  # noqa: F401

import runpy
import sys
from pathlib import Path

STEPS = ["00_build_master", "01_harmonize", "02_enrich", "03_aggregate",
         "04_compare_norm", "05_cellmap", "06_figures", "07_gene_lists", "08_explorer", "09_report_schemes", "10_tissue_split", "11_pathway_components", "12_multiomic", "13_pathway_network",
         "14_layer_reports", "15_cellmap_multiomic",
         ("15_cellmap_multiomic", "--no-phospho")]

here = Path(__file__).parent
for step in STEPS:
    print(f"\n{'=' * 70}\n{step}\n{'=' * 70}")
    sys.argv = [str(here / f"{step}.py")]
    runpy.run_path(str(here / f"{step}.py"), run_name="__main__")
