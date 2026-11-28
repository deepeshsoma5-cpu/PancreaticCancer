# builds docs/MODEL_CARD.md and tests code/score_new_cohort.R against the saved pipeline scores
f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
mdt <- function(df) {
  df[] <- lapply(df, as.character)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
}
sig <- read.csv("results/10_locked_signature.csv")
sig <- sig[order(-abs(sig$coef)), ]
v <- read.csv("results/25_fresh_scores.csv")
v <- v[v$meta %in% TRUE, ]
perf <- data.frame(cohort = v$cohort, n = v$n, events = v$events,
                   C_index_95CI = paste0(f3(v$cindex), " (", f3(v$ci_lo), "-", f3(v$ci_hi), ")"),
                   HR_per_SD_95CI = paste0(f2(v$hr_sd), " (", f2(v$hr_lo), "-", f2(v$hr_hi), ")"),
                   meets_replication_criterion = v$replicates, stringsAsFactors = FALSE)
card <- c(
  "# Model card: 8-gene prognostic score for pancreatic ductal adenocarcinoma (research prototype)", "",
  "## What it is", "",
  "A weighted sum of the within-cohort z-scored expression of 8 genes, with weights from a Cox-LASSO model fit on TCGA-PAAD (145 patients, 87 deaths). Higher score means shorter predicted overall survival. The signature was locked in commit 521df7d before any external cohort was scored.", "",
  "## Genes and weights", "",
  mdt(data.frame(gene = sig$gene, coefficient = formatC(sig$coef, format = "f", digits = 4), stringsAsFactors = FALSE)),
  "## How to compute it", "",
  "Take a cohort of patients with log-scale expression for the 8 genes. Z-score each gene across the patients of that cohort, multiply by the weights above, and sum. Rscript code/score_new_cohort.R expression.csv scores.csv does this and also assigns a high or low group by the cohort median.", "",
  "Important: the score is defined relative to the cohort being scored. It cannot be computed for a single patient, and the high/low split depends on who else is in the cohort.", "",
  "## Performance (each external cohort scored once)", "",
  mdt(perf),
  "The replication criterion (lower bound of the bootstrap 95 percent CI of the C-index above 0.5 and HR per SD above 1) was written down after GSE62452 and GSE21501 had been scored and before GSE57495 was scored. See docs/INTEGRITY_REPORT.md.", "",
  "## Intended use", "",
  "Research and teaching only. This is a student science-fair prototype. It must not be used for any clinical or patient decision.", "",
  "## Known limitations", "",
  "- Discrimination is weak: external C-index about 0.56 to 0.59 with wide confidence intervals.",
  "- Small validation cohorts (59 to 102 patients).",
  "- Training and validation patients were resected pancreatic cancers, mostly early stage in TCGA, so it says nothing about unresectable disease.",
  "- RNA-seq and microarray data were harmonized only by within-cohort z-scoring.",
  "- Associations only. Nothing here shows the genes cause the outcome.",
  "- See results/EXTENDED_RESULTS.md for how it compares with proliferation scores, single genes and random gene sets.", "")
writeLines(card, "docs/MODEL_CARD.md")
cat("wrote docs/MODEL_CARD.md\n")

# test the scoring tool against the scores saved by the original pipeline (GSE62452)
ex <- readRDS("data/gse62452_tumors.rds")
inp <- data.frame(gene = sig$gene, ex$expr[sig$gene, , drop = FALSE], check.names = FALSE)
write.csv(inp, "results/example_input_gse62452_8genes.csv", row.names = FALSE)
out <- system2("Rscript", c("code/score_new_cohort.R", "results/example_input_gse62452_8genes.csv",
                            "results/example_scores_gse62452.csv"), stdout = TRUE, stderr = TRUE)
cat(out, sep = "\n")
now <- read.csv("results/example_scores_gse62452.csv", stringsAsFactors = FALSE)
old <- read.csv("results/10_risk_scores_external.csv", stringsAsFactors = FALSE)
mm <- merge(old, now, by = "sample")
ok <- nrow(mm) == nrow(old) && max(abs(mm$risk - mm$risk_score)) < 1e-8 && all(mm$grp == (mm$group == "high"))
cat("SCORING TOOL MATCHES THE ORIGINAL PIPELINE SCORES:", ok, "\n")
