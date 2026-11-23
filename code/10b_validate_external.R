# one shot test on gse62452 using the signature locked in 10a
source("code/helpers.R")
set.seed(42)
suppressMessages(library(timeROC))

d <- readRDS("data/model_input.rds")
xe <- d$Xe
ce <- d$clin_e
sig <- read.csv("results/10_locked_signature.csv")
ce <- ce[match(rownames(xe), rownames(ce)), ]
stopifnot(all(rownames(xe) == rownames(ce)))

risk_e <- as.numeric(xe[, sig$safe_name] %*% sig$coef)
df_e <- data.frame(sample = rownames(xe), time_months = ce$time_months, status = ce$status,
                   risk = risk_e, grp = risk_e > median(risk_e),
                   stage_num = stage_to_num(ce$stage))
cat("external patients:", nrow(df_e), " events:", sum(df_e$status), "\n")

# --- c-index with bootstrap ci and a shuffle test ---
obs <- cidx(df_e$risk, df_e$time_months, df_e$status)
bs <- replicate(2000, {
  i <- sample(nrow(df_e), replace = TRUE)
  tryCatch(cidx(df_e$risk[i], df_e$time_months[i], df_e$status[i]), error = function(e) NA)
})
ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)

# shuffle who lived and died, see how often luck matches the real c-index
perm <- replicate(5000, {
  j <- sample(nrow(df_e))
  cidx(df_e$risk, df_e$time_months[j], df_e$status[j])
})
p_perm <- (1 + sum(perm >= obs)) / (length(perm) + 1)

cat("\n===== EXTERNAL VALIDATION (GSE62452), LOCKED SIGNATURE =====\n")
cat(sprintf("c-index %.3f | 95%% bootstrap ci %.3f to %.3f | shuffle test p = %.4f\n",
            obs, ci[1], ci[2], p_perm))

# --- risk groups, km, log-rank, hazard ratio (median of external risk scores) ---
cat("\nlog-rank p (median split):", signif(logrank_p(df_e), 3), "\n")
cat("hazard ratio high vs low (95% ci):", hr_ci(df_e), "\n")
tryCatch(save_km(df_e, "results/10_km_external.png"),
         error = function(e) cat("km plot failed:", conditionMessage(e), "\n"))

# same thing but cutting at the tcga median instead
tr <- read.csv("results/10_risk_scores_tcga.csv")
df_e2 <- df_e
df_e2$grp <- df_e$risk > median(tr$risk)
cat("log-rank p using tcga median cutoff:", signif(logrank_p(df_e2), 3),
    " (high-risk n =", sum(df_e2$grp), ")\n")

# --- time dependent roc ---
roc <- tryCatch(timeROC(T = df_e$time_months, delta = df_e$status, marker = df_e$risk,
                        cause = 1, times = c(12, 24, 36), iid = FALSE),
                error = function(e) { cat("timeROC failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(roc)) {
  cat("\ntime dependent auc at 12, 24, 36 months:\n")
  print(round(roc$AUC, 3))
}

# --- integrated brier score (cox on tcga risk, scored on external) ---
tryCatch({
  suppressMessages(library(pec))
  m <- coxph(Surv(time_months, status) ~ risk, data = tr, x = TRUE, y = TRUE)
  pe <- pec(list(signature = m), formula = Surv(time_months, status) ~ 1, data = df_e,
            times = seq(3, 36, by = 3), cens.model = "marginal", reference = TRUE,
            verbose = FALSE)
  cat("\nintegrated brier score up to 36 months (lower is better):\n")
  print(round(crps(pe, times = 36), 4))
}, error = function(e) cat("brier failed:", conditionMessage(e), "\n"))

# --- does it add anything beyond stage? ---
mv <- df_e[, c("time_months", "status", "stage_num")]
mv$risk_z <- as.numeric(scale(df_e$risk))
mv <- mv[complete.cases(mv), ]
m2 <- coxph(Surv(time_months, status) ~ risk_z + stage_num, data = mv)
s2 <- summary(m2)
cat("\n===== MULTIVARIATE COX (external, n =", nrow(mv), ") =====\n")
print(round(cbind(s2$conf.int[, c(1, 3, 4)], p = s2$coefficients[, 5]), 3))
cat("\nproportional hazards check:\n")
print(cox.zph(m2))

# --- sensitivity: different signature sizes (secondary, the primary is already locked) ---
xt <- d$Xt
ct <- d$clin_t
p <- uni_p(xt, ct$time_months, ct$status)
sel <- names(p)[p < 0.01]
if (length(sel) < 20) sel <- names(sort(p))[1:20]
cat("\n===== SENSITIVITY: signature size =====\n")
for (rng in list(c(5, 6), c(7, 8), c(9, 10))) {
  s3 <- fit_sig(xt[, sel], Surv(ct$time_months, ct$status), rng[1], rng[2])
  r3 <- as.numeric(xe[, s3$genes] %*% s3$beta)
  g3 <- data.frame(time_months = df_e$time_months, status = df_e$status, grp = r3 > median(r3))
  cat(sprintf("genes %d | c-index %.3f | log-rank p %.3f\n", length(s3$genes),
              cidx(r3, df_e$time_months, df_e$status), logrank_p(g3)))
}

write.csv(df_e, "results/10_risk_scores_external.csv", row.names = FALSE)
