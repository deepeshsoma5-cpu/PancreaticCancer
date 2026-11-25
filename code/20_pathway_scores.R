# pathway activity scores (mean z of set genes) vs continuous risk score, unadjusted and adjusted for stromal+immune proxies
# rule set before running: all 50 hallmark sets, bh fdr per cohort, replicated = same direction + fdr < 0.05 in tcga and gse62452
# camera (script 19) stays the primary result, this is a complementary self-contained test
suppressMessages({ library(msigdbr); library(MCPcounter) })
sig <- read.csv("results/10_locked_signature.csv")
tc <- readRDS("data/tcga_paad_pdac.rds"); ex <- readRDS("data/gse62452_tumors.rds"); s21 <- readRDS("data/gse21501_tumors.rds")
rt <- read.csv("results/10_risk_scores_tcga.csv"); re <- read.csv("results/10_risk_scores_external.csv")
have <- intersect(sig$gene, rownames(s21$expr))
z <- scale(t(s21$expr[have, , drop = FALSE])); z[!is.finite(z)] <- 0
r21 <- as.numeric(z %*% sig$coef[match(have, sig$gene)])
dropsig <- function(x) x[!rownames(x) %in% sig$gene, , drop = FALSE]
imm <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils")
m <- as.data.frame(tryCatch(msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) msigdbr(species = "Homo sapiens", category = "H")))
sets <- split(m$gene_symbol, m$gs_name)
coh <- list(tcga = list(x = dropsig(tc$logcpm)[, rt$patient], r = rt$risk),
            ext = list(x = dropsig(ex$expr)[, re$sample], r = re$risk),
            g21 = list(x = dropsig(s21$expr), r = r21))
res <- list()
for (nm in names(coh)) {
  x <- coh[[nm]]$x; rz <- as.numeric(scale(coh[[nm]]$r))
  zx <- t(scale(t(x))); zx[!is.finite(zx)] <- 0
  zp <- t(scale(t(MCPcounter.estimate(x, featuresType = "HUGO_symbols"))))
  st <- as.numeric(scale(colMeans(zp[c("Endothelial cells", "Fibroblasts"), ]))); im <- as.numeric(scale(colMeans(zp[imm, ])))
  o <- data.frame()
  for (s in names(sets)) {
    g <- intersect(sets[[s]], rownames(zx)); if (length(g) < 10) next
    sc <- as.numeric(scale(colMeans(zx[g, , drop = FALSE])))
    ct <- cor.test(rz, sc, method = "spearman", exact = FALSE)
    f <- summary(lm(sc ~ rz + st + im))$coefficients
    o <- rbind(o, data.frame(set = s, rho = unname(ct$estimate), p = ct$p.value, beta_adj = f["rz", 1], p_adj = f["rz", 4]))
  }
  o$fdr <- p.adjust(o$p, "BH"); o$fdr_adj <- p.adjust(o$p_adj, "BH")
  res[[nm]] <- o
}
mg <- Reduce(function(a, b) merge(a, b, by = "set"), lapply(names(res), function(n) { d <- res[[n]]; names(d)[-1] <- paste0(names(d)[-1], "_", n); d }))
old14 <- { o <- read.csv("results/11_hallmark_camera.csv"); o$set[o$replicated] }
mg$was14 <- mg$set %in% old14
mg$rep <- sign(mg$rho_tcga) == sign(mg$rho_ext) & mg$fdr_tcga < 0.05 & mg$fdr_ext < 0.05
mg$rep_adj <- sign(mg$beta_adj_tcga) == sign(mg$beta_adj_ext) & mg$fdr_adj_tcga < 0.05 & mg$fdr_adj_ext < 0.05
mg$set <- sub("HALLMARK_", "", mg$set); mg <- mg[order(mg$fdr_tcga), ]
cat("replicated (tcga + gse62452), unadjusted:", sum(mg$rep), "| of original 14:", sum(mg$rep & mg$was14), "\n")
cat("replicated (tcga + gse62452), adjusted for stromal+immune:", sum(mg$rep_adj), "\n")
print(mg[mg$rep | mg$rep_adj | mg$was14, c("set", "rho_tcga", "fdr_tcga", "rho_ext", "fdr_ext", "fdr_g21", "fdr_adj_tcga", "fdr_adj_ext", "rep", "rep_adj")], digits = 2, row.names = FALSE)
write.csv(mg, "results/20_pathway_scores.csv", row.names = FALSE)
