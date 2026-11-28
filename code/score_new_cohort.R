#!/usr/bin/env Rscript
# scores a cohort with the locked 8 gene signature (results/10_locked_signature.csv)
# usage: Rscript code/score_new_cohort.R expression.csv scores.csv
# expression.csv: first column = gene symbols, other columns = samples, values on a log scale
# needs a whole cohort (ideally 20 or more patients): each gene is z-scored across the cohort
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript code/score_new_cohort.R expression.csv scores.csv")
sig <- read.csv("results/10_locked_signature.csv")[, c("gene", "coef")]
x <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
genes <- as.character(x[[1]])
m <- as.matrix(x[, -1, drop = FALSE])
storage.mode(m) <- "double"
rownames(m) <- genes
m <- m[!duplicated(genes), , drop = FALSE]
miss <- setdiff(sig$gene, rownames(m))
if (length(miss)) stop("missing signature genes: ", paste(miss, collapse = ", "))
if (ncol(m) < 20) warning("fewer than 20 samples: the within-cohort z-scoring is unreliable with so few patients")
z <- scale(t(m[sig$gene, , drop = FALSE]))
if (any(!is.finite(z))) {
  warning("a gene had no variation or missing values, treated as 0")
  z[!is.finite(z)] <- 0
}
risk <- as.numeric(z %*% sig$coef)
out <- data.frame(sample = colnames(m), risk_score = risk, risk_z = as.numeric(scale(risk)),
                  group = ifelse(risk > median(risk), "high", "low"), stringsAsFactors = FALSE)
write.csv(out, args[2], row.names = FALSE)
cat("scored", nrow(out), "samples, wrote", args[2], "\n")
