# Plan: fresh-cohort test of the locked signature (written 2026-11-25, before any new cohort is scored)
- Signature: results/10_locked_signature.csv (8 genes, TCGA-only, locked in commit 521df7d). No refitting, no changes to the model.
- Candidates: GSE57495, GSE78229, GSE79668, GSE28735 (GEO), E-MTAB-6134 (ArrayExpress), ICGC PACA-CA/AU if downloadable.
- Inclusion (fixed now): primary tumors, overall survival available, at least 50 patients and 30 events, all 8 signature genes measured.
- Duplicate screen first: compare each new cohort with TCGA, GSE62452, GSE21501, GSE71729 and with each other (best expression match + identical time and status). Flagged patients are dropped from the new cohort.
- Scoring: z-score genes inside each cohort, weight by the locked coefficients, score once per cohort. Report C-index with bootstrap CI and HR per SD.
- Reported whatever the result. No cohort is dropped after scoring. A failed download or failed inclusion rule is documented and skipped.
- Pooling with GSE62452 and GSE21501 (fixed-effect HR per SD) is a post hoc extension and is labeled that way.

## Additions (written 2026-11-25T19:14Z, before any new cohort is scored)
- Excluded: GSE78229 (same samples as GSE62452) and GSE28735 (only 45 tumors, and its samples are also part of GSE62452 per the GSE62452 series description).
- ICGC PACA-AU/CA: listed under controlled access (DACO). Skipped unless a download route is confirmed. Documented, not silently dropped.
- Definition of "replicates": the 95% bootstrap CI of the C-index excludes 0.5 AND the HR per SD of the locked risk score is above 1. Stated now, applied to every cohort the same way.
- A cohort that fails the duplicate screen, or loses patients to it so that it drops below 50 patients / 30 events, is skipped and documented.

## Duplicate screen rules (written 2026-11-25T19:18Z, before the screen is run on GSE57495 / GSE79668)
- Each gene is z-scored inside its own cohort; the 3000 most variable shared genes (mean MAD rank) are correlated between every new-cohort sample and every sample of each reference cohort (TCGA, GSE62452, GSE21501, GSE71729 all 123, and the other new cohort).
- A new-cohort patient is flagged and dropped if best r >= 0.8, OR if a match has r >= 0.5 AND survival time within 0.5 months AND the same status.
- After dropping, a cohort must still have at least 50 patients, at least 30 events and all 8 signature genes, or it is skipped and documented.
- Scoring happens only after this screen is finished. Nothing is scored in the build or screen scripts.

## Screen outcome (written 2026-11-25T19:25Z, after the screen, before any scoring of the new cohorts)
- GSE57495: 63 patients, 4 flagged and dropped by the pre-specified rule (screen thresholds committed in 383af9f), leaving 59 patients and 38 events with all 8 signature genes. Passes the inclusion rule. Scored once, with an all-63 sensitivity run reported separately.
- GSE79668: 51 patients, 2 flagged, leaving 49 (under 50) and SLC16A2 is missing from its counts. Fails the inclusion rule, skipped, not scored, no rescue attempted.
