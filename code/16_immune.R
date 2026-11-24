# immune cell estimates (mcp-counter) by risk group, in tcga and gse62452
# also a purity check: does risk just track stromal/immune content, and does the hr survive adjusting for it
# uses the locked tcga-only risk scores, nothing here changes the signature
suppressMessages({ library(MCPcounter); library(survival) })
set.seed(42)

sig <- read.csv("results/10_locked_signature.csv")
tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
rt <- read.csv("results/10_risk_scores_tcga.csv")
re <- read.csv("results/10_risk_scores_external.csv")
stopifnot(all(rt$patient %in% colnames(tc$logcpm)), all(re$sample %in% colnames(ex$expr)))

cohorts <- list(TCGA = list(x = tc$logcpm[, rt$patient], r = rt),
                GSE62452 = list(x = ex$expr[, re$sample], r = re))
imm <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage",
         "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils")
stro <- c("Endothelial cells", "Fibroblasts")
chk <- c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2")
fdr <- function(p) p.adjust(p, "BH")

allres <- list()
for (nm in names(cohorts)) {
  x <- cohorts[[nm]]$x
  r <- cohorts[[nm]]$r
  x <- x[!rownames(x) %in% sig$gene, ]
  cat("\n##########", nm, "| patients", ncol(x), "| high risk", sum(r$grp), "\n")

  mcp <- MCPcounter.estimate(x, featuresType = "HUGO_symbols")
  cat("populations estimated:", nrow(mcp), "\n")
  stopifnot(all(c(imm, stro) %in% rownames(mcp)))
  zp <- t(scale(t(mcp)))
  rz <- as.numeric(scale(r$risk))

  # per population: correlation with risk + high vs low
  o <- data.frame(pop = rownames(zp), rho = NA, p_cor = NA, p_wilcox = NA)
  for (i in seq_len(nrow(zp))) {
    ct <- suppressWarnings(cor.test(rz, zp[i, ], method = "spearman"))
    o$rho[i] <- unname(ct$estimate)
    o$p_cor[i] <- ct$p.value
    o$p_wilcox[i] <- wilcox.test(zp[i, ] ~ r$grp)$p.value
  }
  o$fdr_cor <- fdr(o$p_cor)
  o$fdr_wilcox <- fdr(o$p_wilcox)
  cat("\n-- mcp-counter vs risk (rho > 0 = more in high risk) --\n")
  print(o, digits = 3, row.names = FALSE)

  # purity proxies
  st <- as.numeric(scale(colMeans(zp[stro, ])))
  im <- as.numeric(scale(colMeans(zp[imm, ])))
  cat(sprintf("\nrisk vs stromal proxy: rho %.2f (p %.2g) | risk vs immune proxy: rho %.2f (p %.2g)\n",
              cor(rz, st, method = "spearman"), cor.test(rz, st, method = "spearman")$p.value,
              cor(rz, im, method = "spearman"), cor.test(rz, im, method = "spearman")$p.value))

  # does risk stay prognostic after adjusting for the proxies
  f0 <- summary(coxph(Surv(r$time_months, r$status) ~ rz))$coefficients
  f1 <- summary(coxph(Surv(r$time_months, r$status) ~ rz + st + im))$coefficients
  cat(sprintf("risk hr per sd: unadjusted %.2f (p %.3f) | adjusted for stromal+immune proxies %.2f (p %.3f)\n",
              exp(f0[1, 1]), f0[1, 5], exp(f1[1, 1]), f1[1, 5]))

  # immune populations again, adjusting for the stromal proxy (skip the stromal ones, circular)
  o$p_adj_stromal <- NA
  for (i in which(o$pop %in% imm)) {
    o$p_adj_stromal[i] <- summary(lm(zp[i, ] ~ rz + st))$coefficients["rz", 4]
  }
  o$fdr_adj_stromal <- fdr(o$p_adj_stromal)

  # checkpoint genes + cytolytic score
  hv <- intersect(chk, rownames(x))
  cyt <- intersect(c("GZMA", "PRF1"), rownames(x))
  g <- data.frame(gene = hv, diff_high_minus_low = NA, p = NA)
  for (i in seq_along(hv)) {
    g$diff_high_minus_low[i] <- mean(x[hv[i], r$grp]) - mean(x[hv[i], !r$grp])
    g$p[i] <- wilcox.test(x[hv[i], ] ~ r$grp)$p.value
  }
  if (length(cyt) == 2) {
    cs <- colMeans(x[cyt, ])
    g <- rbind(g, data.frame(gene = "cytolytic (GZMA+PRF1)", diff_high_minus_low = mean(cs[r$grp]) - mean(cs[!r$grp]),
                             p = wilcox.test(cs ~ r$grp)$p.value))
  }
  g$fdr <- fdr(g$p)
  cat("\n-- checkpoint genes and cytolytic score (log scale difference) --\n")
  print(g, digits = 3, row.names = FALSE)

  write.csv(o, paste0("results/16_mcp_", nm, ".csv"), row.names = FALSE)
  write.csv(g, paste0("results/16_checkpoint_", nm, ".csv"), row.names = FALSE)
  allres[[nm]] <- o
}

# replication: same direction in both cohorts, fdr < 0.05 in both
a <- allres$TCGA
b <- allres$GSE62452
m <- merge(a[, c("pop", "rho", "fdr_cor")], b[, c("pop", "rho", "fdr_cor")], by = "pop", suffixes = c("_tcga", "_ext"))
m$same_dir <- sign(m$rho_tcga) == sign(m$rho_ext)
m$replicated <- m$same_dir & m$fdr_cor_tcga < 0.05 & m$fdr_cor_ext < 0.05
cat("\n########## replication across cohorts ##########\n")
print(m, digits = 3, row.names = FALSE)
write.csv(m, "results/16_mcp_replication.csv", row.names = FALSE)

png("results/16_mcp_rho.png", width = 1000, height = 600, res = 130)
par(mar = c(4, 11, 2, 1))
mm <- rbind(TCGA = m$rho_tcga, GSE62452 = m$rho_ext)
colnames(mm) <- m$pop
barplot(mm, beside = TRUE, horiz = TRUE, las = 1, xlab = "spearman rho with risk score",
        legend.text = TRUE, args.legend = list(x = "topright", cex = 0.8))
abline(v = 0)
dev.off()
