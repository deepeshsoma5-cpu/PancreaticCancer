# tcga-paad mutations by risk group (locked tcga-only risk scores), exploratory
suppressMessages({ library(TCGAbiolinks); library(survival) })
f <- "data/tcga_paad_maf.rds"
if (!file.exists(f)) {
  q <- GDCquery(project = "TCGA-PAAD", data.category = "Simple Nucleotide Variation",
                data.type = "Masked Somatic Mutation", access = "open")
  GDCdownload(q)
  maf <- GDCprepare(q)
  saveRDS(as.data.frame(maf), f)
}
maf <- readRDS(f)
rt <- read.csv("results/10_risk_scores_tcga.csv")
rt$id <- substr(rt$patient, 1, 12)
maf$id <- substr(maf$Tumor_Sample_Barcode, 1, 12)
cat("maf rows:", nrow(maf), "| patients with mutation data:", length(unique(maf$id)), "\n")
rt <- rt[rt$id %in% maf$id, ]
cat("risk-scored patients with mutation data:", nrow(rt), "| high risk:", sum(rt$grp), "\n")

silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "IGR", "Splice_Region")
ns <- maf[!maf$Variant_Classification %in% silent & maf$id %in% rt$id, ]

out <- data.frame()
for (g in c("KRAS", "TP53", "SMAD4", "CDKN2A")) {
  mut <- rt$id %in% ns$id[ns$Hugo_Symbol == g]
  tb <- table(high = rt$grp, mutated = mut)
  out <- rbind(out, data.frame(gene = g, pct_mut_high = round(100 * mean(mut[rt$grp]), 1),
                               pct_mut_low = round(100 * mean(mut[!rt$grp]), 1), p = fisher.test(tb)$p.value))
}
out$fdr <- p.adjust(out$p, "BH")
print(out, digits = 3, row.names = FALSE)

tmb <- as.numeric(table(factor(ns$id, levels = rt$id)))
cat(sprintf("\nnonsilent mutations per patient: median high %.0f, low %.0f | wilcoxon p %.3g\n",
            median(tmb[rt$grp]), median(tmb[!rt$grp]), wilcox.test(tmb ~ rt$grp)$p.value))
cat("(this is a mutation count from the masked MAF, not a calibrated per-megabase TMB)\n")
write.csv(out, "results/18_mutations.csv", row.names = FALSE)
