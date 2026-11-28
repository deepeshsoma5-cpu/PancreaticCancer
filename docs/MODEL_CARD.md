# Model card: 8-gene prognostic score for pancreatic ductal adenocarcinoma (research prototype)

## What it is

A weighted sum of the within-cohort z-scored expression of 8 genes, with weights from a Cox-LASSO model fit on TCGA-PAAD (145 patients, 87 deaths). Higher score means shorter predicted overall survival. The signature was locked in commit 521df7d before any external cohort was scored.

## Genes and weights

| gene | coefficient |
|---|---|
| ARNTL2 | 0.1788 |
| CYP27A1 | -0.1055 |
| MET | 0.0499 |
| SIPA1L2 | -0.0263 |
| DPYD | 0.0192 |
| SLC16A2 | -0.0161 |
| MICAL1 | -0.0148 |
| CAV2 | 0.0117 |

## How to compute it

Take a cohort of patients with log-scale expression for the 8 genes. Z-score each gene across the patients of that cohort, multiply by the weights above, and sum. Rscript code/score_new_cohort.R expression.csv scores.csv does this and also assigns a high or low group by the cohort median.

Important: the score is defined relative to the cohort being scored. It cannot be computed for a single patient, and the high/low split depends on who else is in the cohort.

## Performance (each external cohort scored once)

| cohort | n | events | C_index_95CI | HR_per_SD_95CI | meets_replication_criterion |
|---|---|---|---|---|---|
| GSE62452 (viewed before) | 65 | 49 | 0.592 (0.506-0.678) | 1.40 (1.07-1.82) | TRUE |
| GSE21501 (fresh) | 102 | 66 | 0.565 (0.472-0.654) | 1.22 (0.94-1.58) | FALSE |
| GSE57495 (fresh, after screen) | 59 | 38 | 0.558 (0.449-0.664) | 1.31 (0.90-1.90) | FALSE |

The replication criterion (lower bound of the bootstrap 95 percent CI of the C-index above 0.5 and HR per SD above 1) was written down after GSE62452 and GSE21501 had been scored and before GSE57495 was scored. See docs/INTEGRITY_REPORT.md.

## Intended use

Research and teaching only. This is a student science-fair prototype. It must not be used for any clinical or patient decision.

## Known limitations

- Discrimination is weak: external C-index about 0.56 to 0.59 with wide confidence intervals.
- Small validation cohorts (59 to 102 patients).
- Training and validation patients were resected pancreatic cancers, mostly early stage in TCGA, so it says nothing about unresectable disease.
- RNA-seq and microarray data were harmonized only by within-cohort z-scoring.
- Associations only. Nothing here shows the genes cause the outcome.
- See results/EXTENDED_RESULTS.md for how it compares with proliferation scores, single genes and random gene sets.

