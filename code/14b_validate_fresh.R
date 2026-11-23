# score the locked signatures on the test cohorts, one time each
source("code/helpers.R")
set.seed(42)
suppressMessages(library(timeROC))

new <- read.csv("results/14_locked_signature_pooled.csv")
old <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]
tests <- list("GSE21501 (fresh, primary test)" = readRDS("data/gse21501_tumors.rds"),
              "GSE62452 (viewed before, secondary)" = readRDS("data/gse62452_tumors.rds"))

# z score the signature genes inside the test cohort, then weight by the locked coefs
score_it <- function(expr, sig) {
  have <- intersect(sig$gene, rownames(expr))
  if (length(have) < nrow(sig)) cat("  missing genes:", setdiff(sig$gene, have), "\n")
  z <- scale(t(expr[have, , drop = FALSE]))
  z[!is.finite(z)] <- 0
  as.numeric(z %*% sig$coef[match(have, sig$gene)])
}

eval_it <- function(label, risk, cl) {
  tm <- cl$time_months
  s <- cl$status
  obs <- cidx(risk, tm, s)
  bs <- replicate(2000, {
    i <- sample(length(tm), replace = TRUE)
    tryCatch(cidx(risk[i], tm[i], s[i]), error = function(e) NA)
  })
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  perm <- replicate(5000, {
    j <- sample(length(tm))
    cidx(risk, tm[j], s[j])
  })
  p <- (1 + sum(perm >= obs)) / (length(perm) + 1)
  df <- data.frame(time_months = tm, status = s, grp = risk > median(risk))
  hr <- hr_ci(df)
  cat(sprintf("%-26s c-index %.3f (95%% ci %.3f-%.3f) | shuffle p %.4f | log-rank p %.4f | HR %.2f (%.2f-%.2f)\n",
              label, obs, ci[1], ci[2], p, logrank_p(df), hr[1], hr[2], hr[3]))
  rc <- tryCatch(timeROC(T = tm, delta = s, marker = risk, cause = 1, times = c(12, 24, 36), iid = FALSE),
                 error = function(e) NULL)
  if (!is.null(rc)) cat("   time dependent auc at 12/24/36 months:", paste(round(rc$AUC, 3), collapse = " / "), "\n")
}

for (nm in names(tests)) {
  ct <- tests[[nm]]
  cl <- ct$clin[colnames(ct$expr), ]
  cat("\n#####", nm, "| patients", nrow(cl), "| events", sum(cl$status), "\n")
  rn <- score_it(ct$expr, new)
  eval_it("new pooled signature", rn, cl)
  ro <- score_it(ct$expr, old)
  eval_it("old tcga-only signature", ro, cl)
  if (grepl("^GSE21501", nm)) {
    df <- data.frame(time_months = cl$time_months, status = cl$status, grp = rn > median(rn))
    tryCatch(save_km(df, "results/14_km_gse21501.png"), error = function(e) cat("km plot failed\n"))
  }
}
