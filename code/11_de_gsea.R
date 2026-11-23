# what is different between high and low risk patients?
# differential expression + hallmark gsea, done in both cohorts
# only stuff that shows up in both the same direction counts as a finding
# the 8 signature genes are left out because the groups were built from them
library(limma)
library(edgeR)
library(clusterProfiler)
library(msigdbr)
set.seed(42)

tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
sig <- read.csv("results/10_locked_signature.csv")
rt <- read.csv("results/10_risk_scores_tcga.csv")
re <- read.csv("results/10_risk_scores_external.csv")
stopifnot(all(rt$patient %in% colnames(tc$counts)))
stopifnot(all(re$sample %in% colnames(ex$expr)))

# tcga, voom on the raw counts
cn <- tc$counts[, rt$patient]
y <- calcNormFactors(DGEList(cn))
grp <- factor(ifelse(rt$grp, "high", "low"), levels = c("low", "high"))
des <- model.matrix(~ grp)
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
dt <- topTable(fit, coef = "grphigh", number = Inf, sort.by = "none")
dt <- dt[!rownames(dt) %in% sig$gene, ]

# external, plain limma on the log2 microarray values
xe <- ex$expr[, re$sample]
grp_e <- factor(ifelse(re$grp, "high", "low"), levels = c("low", "high"))
fit_e <- eBayes(lmFit(xe, model.matrix(~ grp_e)))
de <- topTable(fit_e, coef = "grp_ehigh", number = Inf, sort.by = "none")
de <- de[!rownames(de) %in% sig$gene, ]

write.csv(dt, "results/11_de_tcga.csv")
write.csv(de, "results/11_de_external.csv")

# do the two cohorts agree?
both <- intersect(rownames(dt), rownames(de))
cat("genes tested in both:", length(both), "\n")
cr <- cor(dt[both, "t"], de[both, "t"], method = "spearman")
cat("spearman correlation of t stats between cohorts:", round(cr, 3), "\n")
sig_t <- rownames(dt)[dt$adj.P.Val < 0.05]
cat("genes de in tcga (fdr < 0.05):", length(sig_t), "\n")
g <- intersect(sig_t, both)
g <- g[de[g, "P.Value"] < 0.05 & sign(de[g, "logFC"]) == sign(dt[g, "logFC"])]
cat("of those, same direction and p < 0.05 in external:", length(g), "\n")

# hallmark gsea, ranked by t stat. positive nes = more active in high risk
m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
t2g <- unique(as.data.frame(m)[, c("gs_name", "gene_symbol")])

run_gsea <- function(tt) {
  gl <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  r <- GSEA(gl, TERM2GENE = t2g, minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
  as.data.frame(r)[, c("ID", "NES", "pvalue", "p.adjust")]
}
gt <- run_gsea(dt)
ge <- run_gsea(de)
mg <- merge(gt, ge, by = "ID", suffixes = c("_tcga", "_ext"))
mg$replicated <- sign(mg$NES_tcga) == sign(mg$NES_ext) &
  mg$p.adjust_tcga < 0.05 & mg$p.adjust_ext < 0.05
mg <- mg[order(mg$p.adjust_tcga), ]
write.csv(mg, "results/11_hallmark_gsea.csv", row.names = FALSE)

cat("\n===== HALLMARK PATHWAYS (positive NES = up in high risk) =====\n")
show <- mg
show$ID <- sub("HALLMARK_", "", show$ID)
print(head(show, 20), digits = 2, row.names = FALSE)
cat("\npathways significant in both cohorts, same direction:", sum(mg$replicated), "\n")

top <- head(mg[mg$p.adjust_tcga < 0.05, ], 15)
if (nrow(top) > 0) {
  png("results/11_hallmark_nes.png", width = 1800, height = 1200, res = 200)
  par(mar = c(4, 16, 2, 1))
  barplot(t(as.matrix(top[, c("NES_tcga", "NES_ext")])), beside = TRUE, horiz = TRUE,
          names.arg = sub("HALLMARK_", "", top$ID), las = 1, cex.names = 0.7,
          legend.text = c("tcga", "external"), xlab = "NES (positive = up in high risk)")
  dev.off()
}
