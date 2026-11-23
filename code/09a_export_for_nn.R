# export data, folds and per-fold gene lists so the python net can use them
library(survival)
set.seed(7)
nrep <- 10

d <- readRDS("data/model_input.rds")
xt <- d$Xt
clin <- d$clin_t
time <- clin$time_months
status <- clin$status

uni_p <- function(x, t, s) {
  apply(x, 2, function(g) {
    f <- tryCatch(suppressWarnings(coxph(Surv(t, s) ~ g)), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[1, 5]
  })
}

make_folds <- function(s, nfold) {
  f <- integer(length(s))
  for (v in unique(s)) {
    i <- which(s == v)
    f[i] <- sample(rep_len(1:nfold, length(i)))
  }
  f
}

write.csv(xt, "data/nn_X.csv")
meta <- data.frame(patient = rownames(xt), time_months = time, status = status)
write.csv(meta, "data/nn_meta.csv", row.names = FALSE)

fold_rows <- list()
sel_rows <- list()
n <- 0
for (r in 1:nrep) {
  fold <- make_folds(status, 5)
  fold_rows[[r]] <- data.frame(rep = r, patient = rownames(xt), fold = fold)
  for (k in 1:5) {
    tr <- which(fold != k)
    p <- uni_p(xt[tr, ], time[tr], status[tr])
    sel <- names(p)[p < 0.01]
    if (length(sel) < 20) sel <- names(sort(p))[1:20]
    sel <- head(names(sort(p[sel])), 100)
    n <- n + 1
    sel_rows[[n]] <- data.frame(rep = r, fold = k, gene = sel)
  }
  cat("repeat", r, "done\n")
}
write.csv(do.call(rbind, fold_rows), "data/nn_folds.csv", row.names = FALSE)
write.csv(do.call(rbind, sel_rows), "data/nn_selected.csv", row.names = FALSE)
