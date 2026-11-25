# Plan: fresh-cohort test of the locked signature (written 2026-11-25, before any new cohort is scored)
- Signature: results/10_locked_signature.csv (8 genes, TCGA-only, locked in commit 521df7d). No refitting, no changes to the model.
- Candidates: GSE57495, GSE78229, GSE79668, GSE28735 (GEO), E-MTAB-6134 (ArrayExpress), ICGC PACA-CA/AU if downloadable.
- Inclusion (fixed now): primary tumors, overall survival available, at least 50 patients and 30 events, all 8 signature genes measured.
- Duplicate screen first: compare each new cohort with TCGA, GSE62452, GSE21501, GSE71729 and with each other (best expression match + identical time and status). Flagged patients are dropped from the new cohort.
- Scoring: z-score genes inside each cohort, weight by the locked coefficients, score once per cohort. Report C-index with bootstrap CI and HR per SD.
- Reported whatever the result. No cohort is dropped after scoring. A failed download or failed inclusion rule is documented and skipped.
- Pooling with GSE62452 and GSE21501 (fixed-effect HR per SD) is a post hoc extension and is labeled that way.
