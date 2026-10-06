#!/bin/bash
# End-to-end reproduction of the analysis.
#
# Inputs: $HCC_HF_RAW holds the L1000 release files, the discovery series
# matrices, the DE tables, the platform probe->gene maps and the Broad Hub file
# (see README). 01_download_geo.sh fetches the validation cohorts into
# $HCC_HF_VAL. Every stochastic step uses SEED = 20260925, so re-running
# reproduces the reported numbers.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "== 1. differential expression =="
python3 dea.py                          # HCC (GSE14520) + heart failure (GSE57345)
python3 dea_cvd.py                      # CAD (GSE113079)

echo "== 2. signatures, enrichment, connectivity =="
python3 prepare_signatures.py           # disease_signatures.json
python3 run_enrichment.py               # enrichment_results.json
python3 connectivity_all3.py            # connectivity_all3.csv, drug_scores_all3.csv

echo "== 3. score reproduction and limma comparison =="
python3 A1b_reconcile_observed.py
python3 A2_export_expression.py
Rscript A2_limma_run.R
python3 A2_compare.py

echo "== 4. independent-cohort validation =="
bash 01_download_geo.sh
python3 A4_export_cohorts.py
Rscript A4_limma_run.R
python3 A4_concordance.py

echo "== 5. calibration and robustness =="
python3 A1_permutation_null.py 1000
python3 A3_size_matched_sensitivity.py
python3 A5_mechanism_enrichment.py

echo "== 6. chain =="
python3 A12_revised_chain.py

echo "== 7. figures and tables =="
python3 A13_workflow_figure.py
python3 S12_figure1.py
python3 S14_figure2.py
python3 S13_figure3.py
python3 A14_validation_figure.py
python3 S15_figureS1.py
python3 A15_table1.py
python3 S29_supp_excel.py

echo "done. Outputs in output/, figures/, tables/ and supplementary/"
