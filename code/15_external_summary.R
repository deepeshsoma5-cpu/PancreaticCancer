# exploratory: locked tcga-only signature, hr per 1 sd of risk in each external cohort, then combined
# no tuning, nothing here changes the locked signature
suppressMessages(library(survival))
sig <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]
tests <- list(GSE62452 = readRDS("data/gse62452_tumors.rds"), GSE21501 = readRDS("data/gse21501_tumors.rds"))

out <- data.frame()
for (nm in names(tests)) {
  ct <- tests[[nm]]
  cl <- ct$clin[colnames(ct$expr), ]
  have <- intersect(sig$gene, rownames(ct$expr))
  z <- scale(t(ct$expr[have, , drop = FALSE]))
  z[!is.finite(z)] <- 0
  r <- as.numeric(z %*% sig$coef[match(have, sig$gene)])
  r <- as.numeric(scale(r))
  f <- coxph(Surv(cl$time_months, cl$status) ~ r)
  s <- summary(f)$coefficients
  out <- rbind(out, data.frame(cohort = nm, n = nrow(cl), events = sum(cl$status),
                               b = s[1, 1], se = s[1, 3], p = s[1, 5]))
}
w <- 1 / out$se^2
bc <- sum(w * out$b) / sum(w)
sc <- sqrt(1 / sum(w))
q <- sum(w * (out$b - bc)^2)
out$hr <- exp(out$b); out$lo <- exp(out$b - 1.96 * out$se); out$hi <- exp(out$b + 1.96 * out$se)
print(out[, c("cohort", "n", "events", "hr", "lo", "hi", "p")], digits = 3, row.names = FALSE)
cat(sprintf("\ncombined hr per sd: %.2f (%.2f-%.2f) | p %.4f | heterogeneity p %.3f\n",
            exp(bc), exp(bc - 1.96 * sc), exp(bc + 1.96 * sc), 2 * pnorm(-abs(bc / sc)), pchisq(q, 1, lower.tail = FALSE)))
