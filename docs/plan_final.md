# Final plan: close validation and finish pathways (written 2026-11-27T19:59Z, before code/30 and code/31 are run)

## Validation is closed
- The locked 8 gene signature was scored once on three external cohorts (GSE62452, GSE21501, GSE57495 after the duplicate screen). GSE79668 was skipped by the inclusion rule.
- E-MTAB-6134, ICGC and the LASSO + RSF ensemble idea are NOT run. This decision was made after seeing the scores, so it is stated here. Stopping at three cohorts cannot inflate the result.
- Headline summary: the post hoc fixed effect pooled HR per SD, shown next to the per cohort table.

## Pathway analysis (primary)
- Unit of inference is the patient. Genes are z-scored inside each cohort. Pathway score = mean z of the member genes of each MSigDB Hallmark set, with the 8 signature genes removed. Sets need at least 10 genes in every cohort.
- Spearman correlation between pathway score and the locked risk score in each cohort, BH FDR across sets within a cohort.
- Discovery cohort: TCGA-PAAD. External cohorts: GSE62452, GSE21501, GSE57495 (after screen), combined by fixed effect meta analysis of Fisher z correlations (weights n minus 3), BH FDR across sets.
- A pathway is called replicated only if TCGA FDR < 0.05, external meta FDR < 0.05, and the sign is the same in TCGA and in all three external cohorts.
- Calibration: the risk scores are shuffled within each cohort 300 times and the number of replicated sets under the null is reported.
- Also reported: gene level Spearman correlation with the risk score, compared between cohorts.
- Earlier pathway analyses (cameraPR, camera, adjusted score tests) stay as exploratory or sensitivity analyses and are not used as primary evidence. Results are associations, not causal. Tumor purity is an unresolved confounder and is listed as a limitation.
