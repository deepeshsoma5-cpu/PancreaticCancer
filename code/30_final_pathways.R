# final pathway analysis, rules are in docs/plan_final.md (committed before this ran)
# unit of inference is the patient, pathway score = mean z of the member genes
library(msigdbr)
set.seed(42)

sig <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]
tc <- readRDS("data/tcga_paad_pdac.rds")
get <- function(f) list(expr = readRDS(f)$expr)
cohorts <- list(
  tcga = list(expr = tc$logcpm),
  gse62452 = get("data/gse62452_tumors.rds"),
  gse21501 = get("data/gse21501_tumors.rds"),
  gse57495 = get("data/fresh_gse57495.rds"))
ext <- c("gse62452", "gse21501", "gse57495")

# z score every gene inside its own cohort, patients x genes
zs <- function(m) {
  z <- scale(t(m))
  z[!is.finite(z)] <- 0
  z
}

# hallmark sets without the 8 signature genes
m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
m <- as.data.frame(m)
sets <- split(m$gene_symbol, m$gs_name)
sets <- lapply(sets, function(g) setdiff(unique(g), sig$gene))

# locked risk score in each cohort, same recipe as the validation
objs <- lapply(cohorts, function(co) {
  z <- zs(co$expr)
  stopifnot(all(sig$gene %in% colnames(z)))
  list(z = z, risk = as.numeric(z[, sig$gene, drop = FALSE] %*% sig$coef), n = nrow(z))
})
cat("patients per cohort:", paste(names(objs), sapply(objs, function(o) o$n), collapse = ", "), "\n")

ok_set <- sapply(names(sets), function(nm) all(sapply(objs, function(o) sum(sets[[nm]] %in% colnames(o$z)) >= 10)))
keep <- names(sets)[ok_set]
cat("hallmark sets kept (at least 10 genes in every cohort):", length(keep), "of", length(sets), "\n")

# pathway score per patient
for (nm in names(objs)) {
  o <- objs[[nm]]
  o$ps <- sapply(keep, function(s) rowMeans(o$z[, intersect(sets[[s]], colnames(o$z)), drop = FALSE]))
  objs[[nm]] <- o
}

# spearman rho and a t approximation p value
rho_p <- function(ps, risk, n) {
  r <- as.numeric(cor(ps, risk, method = "spearman"))
  t <- r * sqrt((n - 2) / (1 - r^2))
  list(rho = r, p = 2 * pt(-abs(t), n - 2))
}

analyze <- function(rk) {
  r <- lapply(names(objs), function(nm) {
    x <- rho_p(objs[[nm]]$ps, rk[[nm]], objs[[nm]]$n)
    x$fdr <- p.adjust(x$p, "BH")
    x
  })
  names(r) <- names(objs)
  rho_ext <- sapply(ext, function(e) r[[e]]$rho)
  colnames(rho_ext) <- paste0("rho_", ext)
  w <- sapply(ext, function(e) objs[[e]]$n - 3)
  zm <- as.numeric(atanh(rho_ext) %*% w) / sum(w)
  p_meta <- 2 * pnorm(-abs(zm * sqrt(sum(w))))
  fdr_meta <- p.adjust(p_meta, "BH")
  same <- apply(sign(cbind(r$tcga$rho, rho_ext)), 1, function(s) all(s == s[1]))
  rep <- r$tcga$fdr < 0.05 & fdr_meta < 0.05 & same
  rep[is.na(rep)] <- FALSE
  list(r = r, rho_ext = rho_ext, rho_meta = tanh(zm), fdr_meta = fdr_meta, same = same, rep = rep)
}

real <- analyze(lapply(objs, function(o) o$risk))
tab <- data.frame(set = sub("HALLMARK_", "", keep),
                  genes_in_tcga = sapply(keep, function(s) sum(sets[[s]] %in% colnames(objs$tcga$z))),
                  rho_tcga = real$r$tcga$rho, fdr_tcga = real$r$tcga$fdr,
                  real$rho_ext, rho_meta = real$rho_meta, fdr_meta = real$fdr_meta,
                  same_sign_all = real$same, replicated = real$rep, row.names = NULL)
write.csv(tab, "results/30_pathways_final.csv", row.names = FALSE)

cat("\n===== HALLMARK PATHWAYS, patient level (rho > 0 = higher in high risk) =====\n")
show <- tab[order(tab$fdr_tcga), c("set", "rho_tcga", paste0("rho_", ext), "rho_meta", "fdr_tcga", "fdr_meta", "replicated")]
num <- sapply(show, is.numeric)
show[num] <- lapply(show[num], function(x) signif(x, 2))
print(head(show, 25), row.names = FALSE)
cat("\nreplicated pathways (", sum(tab$replicated), "):\n", sep = "")
print(tab$set[tab$replicated])

# calibration: shuffle the risk scores inside each cohort and count false replications
B <- 300
nul <- replicate(B, sum(analyze(lapply(objs, function(o) sample(o$risk)))$rep))
cat(sprintf("\nnull check, %d shuffles: mean replicated %.2f | 95th percentile %d | max %d | share of shuffles with any: %.3f\n",
            B, mean(nul), as.integer(quantile(nul, 0.95)), max(nul), mean(nul > 0)))

# gene level: do genes track the risk score the same way in different cohorts?
gene_rho <- function(o) {
  g <- setdiff(colnames(o$z), sig$gene)
  r <- as.numeric(cor(o$z[, g], o$risk, method = "spearman"))
  names(r) <- g
  r
}
gr <- lapply(objs, gene_rho)
cat("\n===== GENE LEVEL AGREEMENT (spearman of gene-wise rho with the risk score) =====\n")
for (e in ext) {
  sh <- intersect(names(gr$tcga), names(gr[[e]]))
  cat(sprintf("tcga vs %-9s %.2f over %d genes\n", e, cor(gr$tcga[sh], gr[[e]][sh], method = "spearman"), length(sh)))
}

# figure
top <- head(tab[order(tab$fdr_tcga), ], 20)
rr <- as.matrix(top[, c("rho_tcga", paste0("rho_", ext))])
y <- rev(seq_len(nrow(rr)))
png("results/30_pathways_final.png", width = 1900, height = 1400, res = 200)
par(mar = c(4, 14, 2, 1))
plot(NULL, xlim = range(rr, na.rm = TRUE) * 1.1, ylim = c(1, nrow(rr)), yaxt = "n",
     xlab = "spearman rho between pathway score and locked risk score", ylab = "")
axis(2, at = y, labels = top$set, las = 1, cex.axis = 0.7)
abline(v = 0, lty = 2)
for (j in 1:4) points(rr[, j], y, pch = c(16, 1, 2, 5)[j], cex = 1.1)
legend("bottomright", legend = c("tcga", ext), pch = c(16, 1, 2, 5), cex = 0.8, bty = "n")
dev.off()
