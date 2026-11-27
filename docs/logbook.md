# Research logbook

Dates below are UTC and reconstructed from the repository's timestamped commits
and saved analysis outputs. The commit hashes are included so each entry can be
audited. Exploratory and post hoc analyses are labeled as such; they did not
change the locked TCGA-only signature.

## November 23, 2026

- **21:22:45 — Pooled-model pilot (`d3c65c9`).** Began a TCGA + GSE71729
  model-selection attempt. This was exploratory and not used as the primary
  signature.
- **21:23:28 — Pilot scoring (`4a00f8a`).** Scored the pilot signatures on
  GSE21501 and GSE62452. These were preliminary results; the later overlap
  investigation changed the training eligibility decision.

## November 24, 2026

- **22:25:38 — Overlap screen (`78ce2b1`).** Identical survival time and event
  status flagged 86 of 123 GSE71729 samples as overlapping GSE21501 despite
  low expression correlations. Removed those patients from pooled training;
  37 patients and 26 events remained. A filename bug in the 14a workflow that
  would have retained all 123 was caught before pooled-model test scoring.
- **22:34:13 — Pooled model locked (`efe0ffc`).** Locked the pooled model
  before test scoring: TCGA (145) plus non-overlapping GSE71729 (37), 182
  patients, 113 events, seven genes. This was an attempted alternative, not a
  replacement for the TCGA-only primary model.
- **22:34:59 — Test scoring (`ac526e2`).** Compared pooled and TCGA-only
  signatures on GSE21501 and GSE62452. The pooled model did not improve the
  fresh-cohort result and was not promoted to primary.
- **22:36:45 — HR per SD (`dc10315`, post hoc).** Estimated continuous-score
  hazard ratios in GSE62452 and GSE21501 and a fixed-effect pooled estimate.
- **23:02:37 — Immune estimates (`d5fa5ca`, exploratory).** Ran MCP-counter
  estimates, checkpoint-gene comparisons, cytolytic-score comparisons, and
  immune/stromal proxy adjustment. The separately attempted
  `immunedeconv` installation was unavailable; MCP-counter was the immune
  deconvolution method used.
- **23:08:45 — Subtype and mutation analyses (`31b9006`, post hoc).** Tested
  Moffitt basal-like/classical subtype associations and TCGA mutation
  associations. These are exploratory; Moffitt overlaps GSE21501.
- **23:11:41 — Mutation adjustment (`47c5aa9`, post hoc).** Adjusted KRAS,
  TP53, and mutation-count associations for the immune proxy.

## November 25, 2026

- **18:45:14 — Pathway robustness (`31af8cf`, exploratory).** Re-ran Hallmark
  enrichment with `camera()` and estimated inter-gene correlation, then
  examined an immune/stromal-adjusted version and GSE21501 as an additional
  exploratory cohort. This did not support the first `cameraPR` pathway
  result under the estimated-correlation method.
- **18:48:42 — Complementary pathway scores (`fe4bd66`, exploratory).**
  Tested mean-expression Hallmark scores against continuous risk. The
  unadjusted score test found concordant associations; the adjusted test was
  treated cautiously because of over-adjustment concerns.
- **18:59:04 — Prospective cohort plan (`de61008`).** Committed the plan to
  test the unchanged TCGA-only signature once in fresh eligible cohorts,
  after a duplicate screen.
- **18:59:17 — Candidate inspection (`51fb346`).** Inspected candidate
  datasets before scoring. GSE78229 was later excluded as the same cohort as
  GSE62452; GSE28735 failed the minimum-size rule.
- **19:14:51 — Plan additions (`70e8749`).** Added candidate exclusions,
  ICGC access limits, and the operational replication rule before scoring;
  inspected GSE57495 and GSE79668 metadata.
- **19:18:32 — Duplicate thresholds fixed (`383af9f`).** Committed the
  expression/survival matching thresholds before screening the two candidate
  cohorts.
- **19:18:45 — Cohort build and screen (`1d219b5`).** Built GSE57495 and
  GSE79668 and ran the pre-specified duplicate screen, without scoring.
- **19:21:05 — GSE79668 label/matching correction (`1ad3706`).** Corrected
  count-matrix gene labels and sample matching; reran the duplicate screen,
  still without scoring.
- **19:25:15 — Screen outcome recorded (`02ba724`).** GSE57495 passed after
  screening (59 patients, 38 events, all signature genes). GSE79668 failed
  (49 after screen and SLC16A2 absent); it was skipped without rescue.
- **19:25:31 — One-time scoring (`e9f8aea`).** Applied the locked signature to
  eligible cohorts under the pre-committed rule. GSE57495 did not meet the
  replication definition. The GSE62452 result is a previously viewed
  reference, not a fresh validation.
- **19:31 — Final report and figure finish.** Re-ran the summary-figure script
  against the saved CV and locked-signature score CSVs. Added an AI-use
  disclosure and explicit E-MTAB-6134/ICGC status to the report. The final
  interpretation remains non-replicating in both fresh cohorts; no test cohort
  was used to alter the locked signature.

## Interpretation guardrails

- The 8-gene coefficients were selected in TCGA and were not refit on external
  cohorts.
- The pathway enrichment result based on `cameraPR` was superseded by the
  estimated-correlation `camera()` analysis; pathway claims must describe the
  robust result and the complementary score test separately.
- HR pooling across cohorts, immune adjustment, subtype, mutation, pathway
  score, and purity-related analyses are exploratory/post hoc unless a
  contemporaneous commit explicitly pre-specified the test.
- Public human datasets were analyzed without participant contact. The
  applicable fair's SRC/IRB determination remains to be confirmed before
  submission.

## Late entries (written 2026-11-27T20:12Z for events on Nov 25 to 27; labeled as late so the log stays honest)

- Nov 25, retracted claim: the first pathway analysis (limma cameraPR, code/11b) reported 14 replicated Hallmark pathways. This was retracted because cameraPR assumes an inter-gene correlation of 0.01, while the data give about 0.04 to 0.08, which inflated the p-values. A rerun with the correlation estimated from the data (31af8cf) found 0 replicated pathways in TCGA and GSE62452 (best FDR 0.061 for G2M checkpoint and E2F targets, directions consistent). A patient-level score test (fe4bd66) found 17 replicated unadjusted. The two kinds of test answer different questions, and the discrepancy is reported in the paper.
- Nov 26 02:46Z: final pathway analysis specified in writing (docs/plan_final.md, ea823fb) and run in the same session. The results commit (59a6bfb) came 18 seconds later, so this is "committed before the run", not a long advance registration. Result: 23 of 50 Hallmark sets replicated across four cohorts, and 0 of 300 shuffled-risk controls replicated anything.
- Nov 26: the decision to stop validation at three cohorts (GSE62452, GSE21501, GSE57495) was made after seeing the scores, and is recorded as such in docs/plan_final.md. GSE79668 was skipped by the inclusion rule (49 patients after the screen, SLC16A2 missing). E-MTAB-6134, ICGC and a LASSO + RSF ensemble were not run.
- Nov 27 19:59Z, mistake and fix: a repeated run of the plan-writing step overwrote docs/plan_final.md with a new timestamp (5ec9d27). It was reverted in 27f2c2f and the original header (2026-11-26T02:46Z) is intact. No analysis results changed.
- Nov 27 20:09Z: removed the "original 14" comparison lines from the exploratory pathway scripts (7ca2010). All values in the columns shared with the committed CSVs were verified identical, and only the was14 column was dropped. Committed the mutation-file download list as docs/gdc_manifest_mutation_files.txt (it had overwritten MANIFEST.txt).
