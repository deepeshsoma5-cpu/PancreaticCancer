# pooled training on tcga + gse71729, cox lasso only
# the test cohorts (gse21501, gse62452) are NOT loaded in this script
source("code/helpers.R")
set.seed(42)
nrep <- as.integer(Sys.getenv("N_REPEATS", "10"))
nvar <- 2000
pcut <- 0.01

tc <- readRDS("data/tcga_paad_pdac.rds")
mo <- readRDS("data/gse71729_tumors.rds")
drop <- if (file.exists("results/14_dropped_moffitt_samples.txt")) readLines("results/14_dropped_moffitt_samples.txt") else character(0)
drop <- drop[drop != ""]
mo$expr <- mo$expr[, !colnames(mo$expr) %in% drop, drop = FALSE]
mo$clin <- mo$clin[colnames(mo$expr), ]

xa <- tc$logcpm
ca <- tc$clin[colnames(xa), ]
xb <- mo$expr
cb <- mo$clin[colnames(xb), ]
genes <- intersect(rownames(xa), rownames(xb))
xa <- xa[genes, ]
xb <- xb[genes, ]
cat("tcga patients:", ncol(xa), " gse71729 patients:", ncol(xb), " shared genes:", length(genes), "\n")

# keep genes above the median expression in both cohorts, then the 2000 most variable
# (no survival info used here)
ma <- rowMeans(xa)
mb <- rowMeans(xb)
ok <- ma > median(ma) & mb > median(mb)
mad_a <- apply(xa[ok, ], 1, mad)
mad_b <- apply(xb[ok, ], 1, mad)
rk <- (rank(-mad_a) + rank(-mad_b)) / 2
top <- names(sort(rk))[1:nvar]

# z score each gene inside its own cohort
za <- scale(t(xa[top, ]))
zb <- scale(t(xb[top, ]))
good <- apply(za, 2, function(v) all(is.finite(v))) & apply(zb, 2, function(v) all(is.finite(v)))
za <- za[, good]
zb <- zb[, good]

X <- rbind(za, zb)
cohort <- c(rep("tcga", nrow(za)), rep("moffitt", nrow(zb)))
time <- c(ca$time_months, cb$time_months)
status <- c(ca$status, cb$status)
cat("pooled patients:", nrow(X), " events:", sum(status), " candidate genes:", ncol(X), "\n")

# survival test per gene, with each cohort getting its own baseline
uni_p2 <- function(x, t, s, co) {
  apply(x, 2, function(g) {
    f <- tryCatch(suppressWarnings(coxph(Surv(t, s) ~ g + strata(co))), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[1, 5]
  })
}
fit_pool <- function(x, t, s, co) {
  p <- uni_p2(x, t, s, co)
  sel <- names(p)[p < pcut]
  if (length(sel) < 20) sel <- names(sort(p))[1:20]
  y <- stratifySurv(Surv(t, s), co)
  fit_sig(x[, sel, drop = FALSE], y, 5, 10)
}
score <- function(x, sg) as.numeric(x[, sg$genes, drop = FALSE] %*% sg$beta)

# train on one cohort, test on the other
cat("\n===== TRAIN ON ONE COHORT, TEST ON THE OTHER =====\n")
for (d in list(c("tcga", "moffitt"), c("moffitt", "tcga"))) {
  tr <- which(cohort == d[1])
  te <- which(cohort == d[2])
  sg <- fit_pool(X[tr, ], time[tr], status[tr], cohort[tr])
  r <- score(X[te, ], sg)
  cat(sprintf("train %-8s -> test %-8s | c-index %.3f | genes %d\n",
              d[1], d[2], cidx(r, time[te], status[te]), length(sg$genes)))
}

# repeated 5 fold cv on the pooled set, c-index measured inside each cohort
folds <- function(s, co, k) {
  f <- integer(length(s))
  key <- paste(co, s)
  for (v in unique(key)) {
    i <- which(key == v)
    f[i] <- sample(rep_len(1:k, length(i)))
  }
  f
}
res <- data.frame()
for (r in 1:nrep) {
  fd <- folds(status, cohort, 5)
  for (k in 1:5) {
    tr <- which(fd != k)
    te <- which(fd == k)
    sg <- fit_pool(X[tr, ], time[tr], status[tr], cohort[tr])
    risk <- score(X[te, ], sg)
    for (co in unique(cohort[te])) {
      i <- cohort[te] == co
      res <- rbind(res, data.frame(rep = r, fold = k, cohort = co,
                                   cindex = cidx(risk[i], time[te][i], status[te][i]),
                                   n_genes = length(sg$genes)))
    }
  }
  cat("done with repeat", r, "of", nrep, "\n")
}
write.csv(res, "results/14_pooled_cv.csv", row.names = FALSE)
cat("\n===== POOLED CV (c-index inside each cohort) =====\n")
for (co in unique(res$cohort)) {
  v <- res$cindex[res$cohort == co]
  cat(sprintf("%-8s mean %.3f | sd %.3f\n", co, mean(v), sd(v)))
}
cat("median genes per signature:", median(res$n_genes), "\n")

# lock it
sg <- fit_pool(X, time, status, cohort)
out <- data.frame(gene = sg$genes, coef = sg$beta)
write.csv(out, "results/14_locked_signature_pooled.csv", row.names = FALSE)
cat("\n===== LOCKED POOLED SIGNATURE =====\n")
print(out[order(-abs(out$coef)), ], row.names = FALSE)

old <- read.csv("results/10_locked_signature.csv")
cat("\ngenes shared with the old tcga-only signature:", paste(intersect(out$gene, old$gene), collapse = ", "), "\n")

# apparent performance (optimistic, same patients)
risk <- score(X, sg)
for (co in unique(cohort)) {
  i <- cohort == co
  df <- data.frame(time_months = time[i], status = status[i], grp = risk[i] > median(risk[i]))
  cat(sprintf("apparent %-8s c-index %.3f | log-rank p %.2g\n", co, cidx(risk[i], time[i], status[i]), logrank_p(df)))
}
