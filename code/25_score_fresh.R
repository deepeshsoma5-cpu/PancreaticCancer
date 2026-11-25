# score the locked 8 gene signature once per cohort and apply the rule from docs/plan_fresh_cohort.md
# replicates = bootstrap 95% ci of the c-index excludes 0.5 AND hr per sd > 1
source("code/helpers.R")
set.seed(42)

sig <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]

get <- function(f) {
  x <- readRDS(f)
  list(expr = x$expr, clin = x$clin[colnames(x$expr), ])
}
cohorts <- list(
  list(name = "GSE62452 (viewed before)", d = get("data/gse62452_tumors.rds"), meta = TRUE, fresh = FALSE),
  list(name = "GSE21501 (fresh)", d = get("data/gse21501_tumors.rds"), meta = TRUE, fresh = TRUE),
  list(name = "GSE57495 (fresh, after screen)", d = get("data/fresh_gse57495.rds"), meta = TRUE, fresh = TRUE),
  list(name = "GSE57495 (all 63, sensitivity)", d = get("data/gse57495_tumors.rds"), meta = FALSE, fresh = TRUE))

# z score the 8 genes inside the cohort, weight by the locked coefs
score_it <- function(expr) {
  stopifnot(all(sig$gene %in% rownames(expr)))
  z <- scale(t(expr[sig$gene, , drop = FALSE]))
  z[!is.finite(z)] <- 0
  as.numeric(z %*% sig$coef)
}

out <- data.frame()
for (ct in cohorts) {
  cl <- ct$d$clin
  risk <- score_it(ct$d$expr)
  z <- as.numeric(scale(risk))
  obs <- cidx(risk, cl$time_months, cl$status)
  bs <- replicate(2000, {
    i <- sample(nrow(cl), replace = TRUE)
    tryCatch(cidx(risk[i], cl$time_months[i], cl$status[i]), error = function(e) NA)
  })
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  s <- summary(coxph(Surv(cl$time_months, cl$status) ~ z))
  out <- rbind(out, data.frame(cohort = ct$name, n = nrow(cl), events = sum(cl$status),
                               cindex = obs, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
                               hr_sd = s$conf.int[1, 1], hr_lo = s$conf.int[1, 3], hr_hi = s$conf.int[1, 4],
                               p = s$coefficients[1, 5], loghr = s$coefficients[1, 1], se = s$coefficients[1, 3],
                               replicates = ci[1] > 0.5 && s$conf.int[1, 1] > 1,
                               meta = ct$meta, fresh = ct$fresh))
}
write.csv(out, "results/25_fresh_scores.csv", row.names = FALSE)

cat("\n===== LOCKED 8 GENE SIGNATURE, SCORED ONCE PER COHORT =====\n")
for (i in seq_len(nrow(out))) {
  o <- out[i, ]
  cat(sprintf("%-34s n %3d events %3d | c-index %.3f (95%% ci %.3f-%.3f) | HR per SD %.2f (%.2f-%.2f) p %.3f | replicates: %s\n",
              o$cohort, o$n, o$events, o$cindex, o$ci_lo, o$ci_hi, o$hr_sd, o$hr_lo, o$hr_hi, o$p, o$replicates))
}

# fixed effect pooling of the hr per sd (post hoc, label it that way in the paper)
meta <- function(d) {
  w <- 1 / d$se^2
  est <- sum(w * d$loghr) / sum(w)
  se <- sqrt(1 / sum(w))
  q <- sum(w * (d$loghr - est)^2)
  k <- nrow(d)
  cat(sprintf("  pooled HR per SD %.2f (95%% ci %.2f-%.2f), p %.4f | heterogeneity p %s\n",
              exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se), 2 * pnorm(-abs(est / se)),
              if (k > 1) sprintf("%.2f", 1 - pchisq(q, k - 1)) else "NA"))
}
cat("\n===== POOLED HR PER SD, POST HOC =====\n")
cat("all three cohorts (GSE62452, GSE21501, GSE57495 after screen):\n")
meta(out[out$meta, ])
cat("fresh cohorts only (GSE21501, GSE57495 after screen):\n")
meta(out[out$meta & out$fresh, ])
cat("\nskipped: GSE79668 (49 patients after the screen and SLC16A2 missing)\n")
