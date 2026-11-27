# pathway robustness: camera() with estimated inter-gene correlation, then adjusted for stromal+immune proxies
# groups come from the locked tcga-only signature, nothing here changes it. gse21501 is an extra exploratory cohort
# rule (set before running): replicated = same direction and fdr < 0.05 in tcga and gse62452
suppressMessages({ library(limma); library(edgeR); library(msigdbr); library(MCPcounter) })
set.seed(42)
sig <- read.csv("results/10_locked_signature.csv")
tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
s21 <- readRDS("data/gse21501_tumors.rds")
rt <- read.csv("results/10_risk_scores_tcga.csv")
re <- read.csv("results/10_risk_scores_external.csv")

# gse21501 risk groups from the locked signature (same scoring as 14b)
have <- intersect(sig$gene, rownames(s21$expr))
z <- scale(t(s21$expr[have, , drop = FALSE])); z[!is.finite(z)] <- 0
r21 <- as.numeric(z %*% sig$coef[match(have, sig$gene)])
g21 <- r21 > median(r21)

imm <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils")
proxy <- function(x) {
  zp <- t(scale(t(MCPcounter.estimate(x, featuresType = "HUGO_symbols"))))
  data.frame(st = as.numeric(scale(colMeans(zp[c("Endothelial cells", "Fibroblasts"), ]))),
             im = as.numeric(scale(colMeans(zp[imm, ]))))
}
m <- as.data.frame(tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                            error = function(e) msigdbr(species = "Homo sapiens", category = "H")))
sets <- split(m$gene_symbol, m$gs_name)
mk_idx <- function(g) { i <- lapply(sets, function(s) which(g %in% s)); i[sapply(i, length) >= 10] }
camres <- function(y, des, idx) {
  r <- camera(y, idx, des, contrast = 2, inter.gene.cor = NA)
  r$set <- rownames(r)
  r[, c("set", "NGenes", "Correlation", "Direction", "PValue", "FDR")]
}
dropsig <- function(x) x[!rownames(x) %in% sig$gene, , drop = FALSE]
mkdes <- function(grp, a, cv) {
  g <- factor(ifelse(grp, "high", "low"), levels = c("low", "high"))
  if (a == "est_cor") model.matrix(~ g) else model.matrix(~ g + cv$st + cv$im)
}

out <- list()
for (a in c("est_cor", "adjusted")) {
  # tcga (voom on counts)
  cn <- dropsig(tc$counts)[, rt$patient]
  y <- calcNormFactors(DGEList(cn))
  cv <- if (a == "adjusted") proxy(dropsig(tc$logcpm)[, rt$patient]) else NULL
  des <- mkdes(rt$grp, a, cv)
  out[[paste("TCGA", a)]] <- camres(voom(y, des), des, mk_idx(rownames(cn)))
  # gse62452
  xe <- dropsig(ex$expr)[, re$sample]
  cv <- if (a == "adjusted") proxy(xe) else NULL
  des <- mkdes(re$grp, a, cv)
  out[[paste("GSE62452", a)]] <- camres(xe, des, mk_idx(rownames(xe)))
  # gse21501 (exploratory; reference-ratio data, so the proxy is less reliable here)
  x21 <- dropsig(s21$expr)
  cv <- if (a == "adjusted") tryCatch(proxy(x21), error = function(e) NULL) else NULL
  if (a == "est_cor" || !is.null(cv)) {
    des <- mkdes(g21, a, cv)
    out[[paste("GSE21501", a)]] <- camres(x21, des, mk_idx(rownames(x21)))
  }
}

w <- function(k, tag) { d <- out[[k]]; if (is.null(d)) return(NULL); names(d) <- c("set", paste0(c("n", "cor", "dir", "p", "fdr"), "_", tag)); d }
for (a in c("est_cor", "adjusted")) {
  parts <- Filter(Negate(is.null), list(w(paste("TCGA", a), "tcga"), w(paste("GSE62452", a), "ext"), w(paste("GSE21501", a), "g21")))
  d <- Reduce(function(x, y) merge(x, y, by = "set"), parts)
  d$rep2 <- d$dir_tcga == d$dir_ext & d$fdr_tcga < 0.05 & d$fdr_ext < 0.05
  if ("dir_g21" %in% names(d)) d$rep3 <- d$rep2 & d$dir_g21 == d$dir_tcga & d$fdr_g21 < 0.05
  d$set <- sub("HALLMARK_", "", d$set)
  d <- d[order(d$fdr_tcga), ]
  cat("\n##########", a, "##########\n")
  cat("median estimated inter-gene correlation: tcga", round(median(d$cor_tcga), 3), "| ext", round(median(d$cor_ext), 3), "(cameraPR assumed 0.01)\n")
  cat("replicated in tcga + gse62452:", sum(d$rep2))
  if ("rep3" %in% names(d)) cat(" | also in gse21501:", sum(d$rep3))
  cat("\n")
  cols <- intersect(c("set", "dir_tcga", "fdr_tcga", "dir_ext", "fdr_ext", "dir_g21", "fdr_g21", "rep2"), names(d))
  print(d[d$rep2, cols, drop = FALSE], digits = 2, row.names = FALSE)
  write.csv(d, paste0("results/19_pathways_", a, ".csv"), row.names = FALSE)
}
