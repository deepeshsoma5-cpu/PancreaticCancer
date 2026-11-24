# moffitt subtype vs the locked tcga-only risk score (exploratory, nothing here changes the signature)
# subtype codes are 1/2 in the file, so label them from marker genes first
suppressMessages(library(survival))
sig <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]
mo <- readRDS("data/gse71729_tumors.rds")
cl <- mo$clin[colnames(mo$expr), ]
x <- mo$expr
dup <- readLines("results/overlap_drop_gse71729.txt")

have <- intersect(sig$gene, rownames(x))
cat("signature genes present:", length(have), "of", nrow(sig), "\n")
z <- scale(t(x[have, , drop = FALSE])); z[!is.finite(z)] <- 0
risk <- as.numeric(z %*% sig$coef[match(have, sig$gene)])
rz <- as.numeric(scale(risk))

# which code is basal-like? use marker genes, no survival
mk <- function(g) { g <- intersect(g, rownames(x)); rowMeans(scale(t(x[g, , drop = FALSE]))) }
bas <- mk(c("KRT5", "KRT6A", "S100A2")); cla <- mk(c("GATA6", "TFF1", "AGR2", "VSIG1"))
chk <- tapply(bas - cla, cl$subtype, mean)
cat("\nmean (basal markers - classical markers) by subtype code:\n"); print(round(chk, 2))
bcode <- names(which.max(chk))
cat("labeling code", bcode, "as basal-like\n")
isb <- cl$subtype == bcode

run <- function(i, lab) {
  cat("\n=====", lab, "| n", sum(i), "| basal-like", sum(isb[i]), "| events", sum(cl$status[i]), "\n")
  r <- rz[i]; b <- isb[i]; t <- cl$time_months[i]; s <- cl$status[i]
  rk <- rank(r); n1 <- sum(b); n0 <- sum(!b)
  auc <- (sum(rk[b]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  cat(sprintf("risk score in basal vs classical: wilcoxon p %.3g | auc for basal %.2f\n", wilcox.test(r ~ b)$p.value, auc))
  f <- summary(coxph(Surv(t, s) ~ b))$coefficients
  cat(sprintf("basal-like vs classical survival: hr %.2f (p %.3f)\n", exp(f[1, 1]), f[1, 5]))
  f0 <- summary(coxph(Surv(t, s) ~ r))$coefficients
  f1 <- summary(coxph(Surv(t, s) ~ r + b))$coefficients
  cat(sprintf("risk hr per sd: alone %.2f (p %.3f) | adjusted for subtype %.2f (p %.3f)\n",
              exp(f0[1, 1]), f0[1, 5], exp(f1[1, 1]), f1[1, 5]))
}
run(rep(TRUE, nrow(cl)), "all 123 moffitt")
run(!rownames(cl) %in% dup, "37 moffitt not duplicated in gse21501")
