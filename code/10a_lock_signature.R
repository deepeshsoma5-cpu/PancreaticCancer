# lock the signature using tcga only
# the external dataset is never loaded in this script
source("code/helpers.R")
set.seed(42)

d <- readRDS("data/model_input.rds")
xt <- d$Xt
ct <- d$clin_t
gmap <- d$gene_map
time <- ct$time_months
status <- ct$status

# gene filter on all 145 patients, then lasso
p <- uni_p(xt, time, status)
sel <- names(p)[p < 0.01]
if (length(sel) < 20) sel <- names(sort(p))[1:20]
cat("genes passing p < 0.01:", length(sel), "\n")

sig <- fit_sig(xt[, sel], Surv(time, status), 5, 10)
out <- data.frame(gene = gmap[sig$genes], safe_name = sig$genes, coef = sig$beta)
write.csv(out, "results/10_locked_signature.csv", row.names = FALSE)

cat("\n===== LOCKED SIGNATURE =====\n")
print(out[order(-abs(out$coef)), ], row.names = FALSE)
cat("(positive coef = higher expression, higher risk)\n")

# risk score in tcga. this is the training score so its optimistic
risk_t <- as.numeric(xt[, sig$genes] %*% sig$beta)
grp_t <- risk_t > median(risk_t)
df_t <- data.frame(patient = rownames(xt), time_months = time, status = status,
                   risk = risk_t, grp = grp_t)
write.csv(df_t, "results/10_risk_scores_tcga.csv", row.names = FALSE)

cat("\n===== TRAINING (TCGA) PERFORMANCE, OPTIMISTIC =====\n")
cat("apparent c-index:", round(cidx(risk_t, time, status), 3), "\n")
cat("log-rank p:", signif(logrank_p(df_t), 3), "\n")
cat("hazard ratio high vs low (95% ci):", hr_ci(df_t), "\n")
tryCatch(save_km(df_t, "results/10_km_training.png"),
         error = function(e) cat("km plot failed:", conditionMessage(e), "\n"))

# does the score add anything beyond stage and age?
mv <- data.frame(time_months = time, status = status,
                 risk_z = as.numeric(scale(risk_t)),
                 stage_num = stage_to_num(ct$ajcc_pathologic_stage),
                 age = ct$age_at_diagnosis / 365.25)
mv <- mv[complete.cases(mv), ]
m <- coxph(Surv(time_months, status) ~ risk_z + stage_num + age, data = mv)
s <- summary(m)
cat("\n===== MULTIVARIATE COX (TCGA, n =", nrow(mv), ") =====\n")
print(round(cbind(s$conf.int[, c(1, 3, 4)], p = s$coefficients[, 5]), 3))
cat("\nproportional hazards check:\n")
print(cox.zph(m))
