# builds results/FINAL_RESULTS.md from the saved result files so no number is typed by hand
f <- function(x, d = 3) formatC(x, format = "f", digits = d)
md <- function(df) {
  df[] <- lapply(df, as.character)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
}
out <- c("# Final results (generated from results/*.csv, do not edit by hand)", "")

# 1. locked signature
s <- read.csv("results/10_locked_signature.csv")
s <- s[order(-abs(s$coef)), ]
out <- c(out, "## Locked signature (TCGA only, 8 genes, positive coef = higher expression, higher risk)", "",
         md(data.frame(gene = s$gene, coef = f(s$coef, 4))))

# 2. model selection
cvs <- function(file, label) {
  r <- read.csv(file)
  a <- aggregate(cindex ~ model, r, function(v) c(mean(v), sd(v)))
  data.frame(run = label, model = a$model, mean_cindex = f(a$cindex[, 1]), sd = f(a$cindex[, 2]))
}
ms <- rbind(cvs("results/first_run/07_cv_results.csv", "first run (variance filter only)"),
            cvs("results/07_cv_results.csv", "final run (median expression filter added)"))
nn <- read.csv("results/09_nn_cv_results.csv")
ms <- rbind(ms, data.frame(run = "final run (median expression filter added)", model = "NeuralNetwork",
                           mean_cindex = f(mean(nn$cindex)), sd = f(sd(nn$cindex))))
out <- c(out, "## TCGA cross-validated C-index (screening inside each fold)", "", md(ms))

# 3. external validation
v <- read.csv("results/25_fresh_scores.csv")
tabv <- data.frame(cohort = v$cohort, n = v$n, events = v$events,
                   cindex_95ci = paste0(f(v$cindex), " (", f(v$ci_lo), "-", f(v$ci_hi), ")"),
                   hr_per_sd_95ci = paste0(f(v$hr_sd, 2), " (", f(v$hr_lo, 2), "-", f(v$hr_hi, 2), ")"),
                   p = f(v$p), meets_rule = v$replicates)
out <- c(out, "## External validation of the locked signature (each cohort scored once)", "", md(tabv))

# 4. pooled hr per sd, post hoc
pool <- function(d) {
  w <- 1 / d$se^2
  est <- sum(w * d$loghr) / sum(w)
  se <- sqrt(1 / sum(w))
  q <- sum(w * (d$loghr - est)^2)
  paste0("HR per SD ", f(exp(est), 2), " (", f(exp(est - 1.96 * se), 2), "-", f(exp(est + 1.96 * se), 2),
         "), p = ", f(2 * pnorm(-abs(est / se)), 4),
         ", heterogeneity p = ", if (nrow(d) > 1) f(1 - pchisq(q, nrow(d) - 1), 2) else "NA",
         ", patients ", sum(d$n), ", events ", sum(d$events))
}
out <- c(out, "## Pooled hazard ratio (post hoc, fixed effect)", "",
         paste0("- All three external cohorts: ", pool(v[v$meta, ])),
         paste0("- Fresh cohorts only: ", pool(v[v$meta & v$fresh, ])), "")

# 5. pathways
p <- read.csv("results/30_pathways_final.csv")
rep <- p[p$replicated, ]
out <- c(out, "## Hallmark pathways associated with the risk score (patient level)", "",
         paste0("Replicated pathways (TCGA FDR < 0.05, external meta FDR < 0.05, same sign in all four cohorts): ", nrow(rep), " of ", nrow(p), "."), "")
if (nrow(rep) > 0) {
  rep <- rep[order(-abs(rep$rho_meta)), ]
  out <- c(out, md(data.frame(pathway = rep$set, direction = ifelse(rep$rho_meta > 0, "higher in high risk", "lower in high risk"),
                              rho_tcga = f(rep$rho_tcga, 2), rho_external_meta = f(rep$rho_meta, 2),
                              fdr_tcga = formatC(rep$fdr_tcga, format = "e", digits = 1),
                              fdr_meta = formatC(rep$fdr_meta, format = "e", digits = 1))))
}

# 6. what the results do and do not support
out <- c(out, "## What these results support", "",
         "- The locked signature points the same way in every external cohort. Only GSE62452 met the pre-registered rule individually.",
         "- The effect is modest (C-index about 0.56 to 0.59). It is not a clinically usable predictor.",
         "- Pathway results are associations at the patient level, replicated across cohorts. They are not causal and tumor purity is an unresolved confounder.", "",
         "## What they do not support", "",
         "- Any C-index near 0.7, any claim that pooled training or the neural network helped, any claim that the signature beats stage, any use of the earlier cameraPR p-values as evidence, any causal wording.")
writeLines(out, "results/FINAL_RESULTS.md")
cat(out, sep = "\n")
