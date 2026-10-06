# Dual-signature connectivity mapping — analysis pipeline

Code for the analysis in *Dual-signature connectivity mapping identifies cardiac
glycosides as dual-effect candidates for hepatocellular carcinoma and heart
failure*. It takes two disease transcriptomes, screens the L1000 perturbation
library for compounds that reverse **both** signatures, and prioritises the
candidates with a versioned annotation source.

Running the scripts in the order given by `run_all.sh` reproduces every number,
figure and table in the paper.

## Run

```bash
pip install -r requirements.txt            # R + limma required separately
export HCC_HF_RAW=/path/to/large_inputs    # inputs, see below
bash run_all.sh
```

Large inputs are not committed. `01_download_geo.sh` fetches the validation
cohorts into `$HCC_HF_VAL` (default `./downloads`). Everything else is placed by
hand in `$HCC_HF_RAW` (default `./large_inputs`):

```
GEO series matrices — GSE14520-GPL3921, GSE57345-GPL11532, GSE113079
L1000 release files — GSE92742 level-5 COMPZ.MODZ .gctx, plus the
                      gene_info / sig_info / pert_info tables
platform maps       — GPL3921_probe2gene.csv, GPL11532_probe2gene.csv,
                      GPL20115_probe2gene.csv
annotation          — repurposing_drugs_20200324.txt (place it in data/, see data/README.md)
```

Every stochastic step uses a fixed seed (`SEED` in `config.py`), so re-running
reproduces the reported numbers.

## Stages

| Stage | Script | What it does |
|---|---|---|
| 1 | `dea.py`, `dea_cvd.py` | Welch's t-test on each discovery contrast (HCC, heart failure, CAD) |
| 2 | `prepare_signatures.py` | up/down signature gene lists |
| 2 | `run_enrichment.py` | GO biological-process enrichment per signature, against the gene sets frozen under `data/genesets/` |
| 2 | `connectivity_all3.py` | reversal score of every L1000 profile, aggregated to compounds |
| 3 | `A1b_reconcile_observed.py` | reproduction of the published compound scores |
| 3 | `A2_export_expression.py`, `A2_limma_run.R`, `A2_compare.py` | limma vs Welch's t-test on the discovery cohorts |
| 4 | `A4_export_cohorts.py`, `A4_limma_run.R`, `A4_concordance.py` | signature reproduction in held-out cohorts |
| 5 | `A1_permutation_null.py` | permutation null for the dual-reversal criterion |
| 5 | `A3_size_matched_sensitivity.py` | signature-size sensitivity |
| 5 | `A5_mechanism_enrichment.py` | mechanism-class enrichment (Fisher, BH) |
| 6 | `A12_revised_chain.py` | the filtering chain and the prioritised candidate set |
| 6 | `A14_validation_figure.py`, `A15_table1.py` | Figure 4 and Table 1 |
| 7 | `A13_workflow_figure.py`, `S12_figure1.py`, `S14_figure2.py`, `S13_figure3.py`, `S15_figureS1.py` | Figures 1-4 and Figure S1 |
| 7 | `S29_supp_excel.py` | supplementary tables S1-S5 (one workbook) |

Outputs land in `output/`, `figures/`, `tables/` and `supplementary/`.

## Inputs

* Discovery — GSE14520 (HCC), GSE57345 (heart failure)
* Comorbidity — GSE113079 (coronary artery disease)
* Validation — GSE76427 (HCC), GSE141910 (heart failure)
* Screening — LINCS L1000 phase-I, GSE92742 level-5 COMPZ.MODZ matrix
* Annotation — Broad Drug Repurposing Hub, release 2020-03-24

## Notes

* `config.py` derives every path from its own location; nothing is hard-coded.
  Where each file lives is fixed, and every script follows it:

  | directory | holds | written by |
  |---|---|---|
  | `$HCC_HF_RAW` | the downloaded inputs **and** the pipeline intermediates (`DEG_*.csv`, `disease_signatures.json`, `enrichment_results.json`, `drug_scores_all3.csv`, `candidate_drugs_*.csv`) | `dea*`, `prepare_signatures`, `run_enrichment`, `connectivity_all3` |
  | `output/` | the analysis outputs (`A1_*`, `A2_*`, `A3_*`, `A4_*`, `A5_*`, `A12_*`) | the `A*` scripts |
  | `results/` | the archived copies of the outputs the paper reports | committed |
  | `data/genesets/` | the frozen GO gene sets | committed |

  So the figure scripts read intermediates from `$HCC_HF_RAW` and analysis
  outputs from `output/`; pointing `$HCC_HF_RAW` at `results/` is not supported.
* Signature sizes: HCC 150 up / 150 down, CAD 150 up / 150 down, heart failure
  25 up / 23 down (|log2 fold change| > 1 and adjusted P < 0.05; the heart-failure
  contrast yields fewer significant genes, so its whole set is used).
* `run_enrichment.py` reads the gene sets frozen under `data/genesets/`, so its output
  does not depend on the Enrichr service as served at run time.
* Two enrichment results appear in the paper and they are not interchangeable.
  `A12_revised_chain.py` tests mechanism-class over-representation against the full set of
  compound groups (19,811); that is the test the manuscript and Supplementary Table S3(B)
  quote. `A5_mechanism_enrichment.py` instead tests against the annotated compounds
  (2,947) and is reported separately. Reading only the A5 output gives a different, and
  weaker, picture of the same classes.
* The candidate chain is a prioritisation filter: it restricts the candidate set
  to compounds with a documented target and a documented indication, and it
  provides no independent evidence for a dual effect. The statistical claim
  rests on the permutation calibration and the mechanism-class enrichment.

## Processed outputs

`results/` holds the outputs the paper reports, produced by the run above:

| File | What it holds |
|---|---|
| `DEG_HCC.csv`, `DEG_CVD.csv`, `DEG_CAD.csv` | differential expression for each discovery contrast |
| `disease_signatures.json` | the up/down signature gene lists |
| `enrichment_results.json` | GO enrichment per signature |
| `drug_scores_all3.csv` | reversal score of every L1000 compound |
| `candidate_drugs_*.csv` | dual-reversal candidates per disease pair |
| `A1_summary.json`, `A1_*.npy`, `A1_candidate_perm_pvalues.csv` | permutation calibration |
| `A1b_recomputed_drug_scores.csv` | reproduction of the published scores |
| `A2_comparison.json`, `A2_limma_signature_*.json` | limma vs Welch's t-test |
| `A3_size_matched.json` | signature-size sensitivity |
| `A4_concordance.json` | signature reproduction in held-out cohorts |
| `A5_enrichment_*.csv` | mechanism-class enrichment |
| `A12_revised_chain.json`, `A12_*_*.csv` | the filtering chain and candidate sets |

The large per-cohort expression matrices and the full per-profile connectivity
table are not shipped; `run_all.sh` regenerates them.

## Licence

Apache License 2.0 — see `LICENSE` and `NOTICE`.

The licence covers the code in this repository. The analysis reads the Broad
Drug Repurposing Hub, the GEO series matrices and the LINCS L1000 release; none
of them is redistributed here — the Hub is provided for non-commercial use only
(see `data/README.md`) and the rest are downloaded from their public sources.
