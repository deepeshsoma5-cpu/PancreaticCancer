# do the kras/tp53/mutation count differences survive adjusting for the immune+stromal proxy? (exploratory)
suppressMessages({ library(MCPcounter) })
sig <- read.csv("results/10_locked_signature.csv")
tc <- readRDS("data/tcga_paad_pdac.rds")
maf <- readRDS("data/tcga_paad_maf.rds")
rt <- read.csv("results/10_risk_scores_tcga.csv")
x <- tc$logcpm[!rownames(tc$logcpm) %in% sig$gene, rt$patient]
zp <- t(scale(t(MCPcounter.estimate(x, featuresType = "HUGO_symbols"))))
imm <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils")
rt$st <- as.numeric(scale(colMeans(zp[c("Endothelial cells", "Fibroblasts"), ])))
rt$im <- as.numeric(scale(colMeans(zp[imm, ])))
rt$rz <- as.numeric(scale(rt$risk))
rt$id <- substr(rt$patient, 1, 12)
maf$id <- substr(maf$Tumor_Sample_Barcode, 1, 12)
rt <- rt[rt$id %in% maf$id, ]
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "IGR", "Splice_Region")
ns <- maf[!maf$Variant_Classification %in% silent & maf$id %in% rt$id, ]
rt$tmb <- as.numeric(table(factor(ns$id, levels = rt$id)))
cat("n:", nrow(rt), "\n")
for (g in c("KRAS", "TP53")) {
  rt$m <- as.numeric(rt$id %in% ns$id[ns$Hugo_Symbol == g])
  a <- summary(glm(m ~ rz, binomial, rt))$coefficients["rz", ]
  b <- summary(glm(m ~ rz + st + im, binomial, rt))$coefficients
  cat(sprintf("%s: or per sd risk alone %.2f (p %.3g) | adjusted %.2f (p %.3g) | immune proxy or %.2f (p %.3g)\n",
              g, exp(a[1]), a[4], exp(b["rz", 1]), b["rz", 4], exp(b["im", 1]), b["im", 4]))
}
a <- summary(lm(log1p(tmb) ~ rz, rt))$coefficients["rz", ]
b <- summary(lm(log1p(tmb) ~ rz + st + im, rt))$coefficients
cat(sprintf("log(1+mutation count): risk alone p %.3g | adjusted p %.3g | immune proxy p %.3g\n", a[4], b["rz", 4], b["im", 4]))
cat("kras mutated by immune-proxy tertile:\n")
rt$m <- as.numeric(rt$id %in% ns$id[ns$Hugo_Symbol == "KRAS"])
print(round(tapply(rt$m, cut(rt$im, quantile(rt$im, 0:3 / 3), include.lowest = TRUE, labels = c("low", "mid", "high")), mean), 2))
