# Extended analyses and documentation (written 2026-11-28T01:44Z, before code/41 to code/43 are run)

Status: every analysis below is secondary or post hoc. The locked 8 gene signature (results/10_locked_signature.csv, commit 521df7d) is never refit or changed. Nothing here replaces the per cohort results in results/25_fresh_scores.csv. All results are reported whatever they show. No multiplicity correction is applied.

## A. Integrity and duplicates
1. code/40 checks commit order and that the locked signature file is unchanged since commit 521df7d.
2. code/41 is an all pairs duplicate audit across TCGA, GSE62452, GSE21501, GSE57495 (all 63) and GSE71729, using the thresholds in docs/plan_fresh_cohort.md. For a sensitivity rescore, a scored cohort patient is dropped if flagged against TCGA, or against a scored cohort earlier in the order GSE62452, GSE21501, GSE57495.

## B. Pooled individual patient analysis of the three scored cohorts (post hoc)
3. Cox model of the locked risk score (z scored within cohort) stratified by cohort: HR per SD.
4. Pooled Kaplan-Meier of within cohort median split groups with a cohort stratified log-rank test and stratified HR (high vs low).
5. Restricted mean survival time difference (high minus low) at 24 months, per cohort with a bootstrap SE, pooled by inverse variance. The horizon of 24 months is fixed here in advance.

## C. Is it just proliferation, or one gene? (post hoc)
6. Comparators scored the same way: the top weighted gene alone, a unit weight version of the 8 genes, the Hallmark G2M checkpoint score and the Hallmark E2F targets score (signature genes removed from the sets).
7. The signature adjusted for the G2M score and, separately, for the E2F score (pooled stratified Cox).
Interpretation rule: the signature is said to carry information beyond proliferation only if its adjusted pooled HR has a 95 percent CI above 1 in BOTH adjusted models. Otherwise the report says it does not show independent information beyond proliferation in these cohorts.

## D. Robustness of the locked signature (no refitting)
8. Leave one gene out scores with the locked weights, per gene HR per SD in each cohort and pooled, and rank based (inverse normal) scaling instead of z scoring.

## E. Random signature benchmarks (5000 draws each)
9. Benchmark A: random 8 gene sets from the 2000 gene candidate pool used for model building (signature genes excluded), with the 8 locked weights randomly permuted over the genes.
10. Benchmark B: random 8 gene sets from the genes that passed the TCGA univariate screen (p < 0.01), unit weights with the sign of the TCGA univariate coefficient, compared with the unit weight version of the locked signature.
Statistics: the pooled stratified HR per SD, and the mean C-index over the three cohorts. Empirical one sided p = (1 + number of random sets at least as large) / (draws + 1). The signature is said to be distinguishable from random sets only if p < 0.05 on the pooled HR. Otherwise the report says it is not distinguishable from random sets drawn from that pool.

## F. Documentation and tools
11. Integrity report, model card, reproduction guide, session info, input checksums, and a scoring script for new cohorts that is tested against the pipeline.
