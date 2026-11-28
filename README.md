# PancreaticCancer

Science-fair project evaluating a locked eight-gene expression score for
overall survival in pancreatic ductal adenocarcinoma.

- [Research report draft](docs/report_draft.md) summarizes methods, results,
  limitations, references, and figure links.
- [Research logbook](docs/logbook.md) records dated decisions and links them
  to commits and saved outputs.
- [Fresh-cohort test plan](docs/plan_fresh_cohort.md) records the pre-specified
  screening, eligibility, and replication rules.

Scripts are numbered in `code/`; saved analysis outputs and figures are in
`results/`. Fresh-cohort eligibility and screening results are recorded before
scoring in the plan and the corresponding result files.

- [Model card](docs/MODEL_CARD.md), [integrity report](docs/INTEGRITY_REPORT.md), [extended results](results/EXTENDED_RESULTS.md) and [how to reproduce](docs/REPRODUCE.md).
- Score a new cohort with the locked signature: Rscript code/score_new_cohort.R expression.csv scores.csv (requirements are in the model card).
