# what differs between high and low risk patients? done in both cohorts
# de with limma (voom for tcga counts), then hallmark sets with cameraPR
# only pathways that go the same way in both cohorts count
# the 8 signature genes are dropped because the groups were built from them
library(limma)
library(edgeR)
library(msigdbr)
set.seed(42)

tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
sig <- read.csv("results/10_locked_signature.csv")
rt <- read.csv("results/10_risk_scores_tcga.csv")
re <- read.csv("results/10_risk_scores_external.csv")
stopifnot(all(rt$patient %in% colnames(tc$counts)))
stopifnot(all(re$sample %in% colnames(ex$expr)))

# tcga
cn <- tc$counts[, rt$patient]
y <- calcNormFactors(DGEList(cn))
grp <- factor(ifelse(rt$grp, "high", "low"), levels = c("low", "high"))
des <- model.matrix(~ grp)
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
dt <- topTable(fit, coef = "grphigh", number = Inf, sort.by = "none")
dt <- dt[!rownames(dt) %in% sig$gene, ]

# external, microarray so no voom
xe <- ex$expr[, re$sample]
grp_e <- factor(ifelse(re$grp, "high", "low"), levels = c("low", "high"))
des_e <- model.matrix(~ grp_e)
fit_e <- eBayes(lmFit(xe, des_e))
de <- topTable(fit_e, coef = "grp_ehigh", number = Inf, sort.by = "none")
de <- de[!rownames(de) %in% sig$gene, ]

write.csv(dt, "results/11_de_tcga.csv")
write.csv(de, "results/11_de_external.csv")

both <- intersect(rownames(dt), rownames(de))
cat("\ngenes tested in both:", length(both), "\n")
cat("spearman correlation of t stats between cohorts:",
    round(cor(dt[both, "t"], de[both, "t"], method = "spearman"), 3), "\n")
sig_t <- rownames(dt)[dt$adj.P.Val < 0.05]
cat("genes de in tcga (fdr < 0.05):", length(sig_t), "\n")
g <- intersect(sig_t, both)
g <- g[de[g, "P.Value"] < 0.05 & sign(de[g, "logFC"]) == sign(dt[g, "logFC"])]
cat("of those, same direction and p < 0.05 in external:", length(g), "\n")

# hallmark sets
m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
m <- as.data.frame(m)
sets <- split(m$gene_symbol, m$gs_name)

run_sets <- function(tt) {
  stat <- setNames(tt$t, rownames(tt))
  idx <- lapply(sets, function(s) which(names(stat) %in% s))
  idx <- idx[sapply(idx, length) >= 10]
  r <- cameraPR(stat, idx, use.ranks = FALSE)
  r$set <- rownames(r)
  r
}
gt <- run_sets(dt)
ge <- run_sets(de)
mg <- merge(gt[, c("set", "NGenes", "Direction", "PValue", "FDR")],
            ge[, c("set", "Direction", "PValue", "FDR")],
            by = "set", suffixes = c("_tcga", "_ext"))
mg$replicated <- mg$Direction_tcga == mg$Direction_ext &
  mg$FDR_tcga < 0.05 & mg$FDR_ext < 0.05
mg <- mg[order(mg$FDR_tcga), ]
write.csv(mg, "results/11_hallmark_camera.csv", row.names = FALSE)

cat("\n===== HALLMARK PATHWAYS (Up = higher in high risk) =====\n")
show <- mg
show$set <- sub("HALLMARK_", "", show$set)
show$PValue_tcga <- signif(show$PValue_tcga, 2)
show$FDR_tcga <- signif(show$FDR_tcga, 2)
show$PValue_ext <- signif(show$PValue_ext, 2)
show$FDR_ext <- signif(show$FDR_ext, 2)
print(head(show, 25), row.names = FALSE)
cat("\npathways significant in both cohorts, same direction:", sum(mg$replicated), "\n")

top <- head(mg[mg$FDR_tcga < 0.05, ], 15)
if (nrow(top) > 0) {
  sc <- function(d, f) ifelse(d == "Up", 1, -1) * -log10(f)
  mat <- rbind(tcga = sc(top$Direction_tcga, top$FDR_tcga),
               external = sc(top$Direction_ext, top$FDR_ext))
  png("results/11_hallmark_pathways.png", width = 1800, height = 1200, res = 200)
  par(mar = c(4, 16, 2, 1))
  barplot(mat, beside = TRUE, horiz = TRUE, names.arg = sub("HALLMARK_", "", top$set),
          las = 1, cex.names = 0.7, legend.text = c("tcga", "external"),
          xlab = "signed -log10(fdr), positive = up in high risk")
  dev.off()
}
