# model selection: cox lasso vs random survival forest
# only using tcga here, the external data stays locked until the end
library(survival)
library(glmnet)
library(randomForestSRC)
set.seed(42)

nrep <- as.integer(Sys.getenv("N_REPEATS", "3"))
nfold <- 5
nvar <- 2000
pcut <- 0.01

tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
shared <- intersect(rownames(tc$logcpm), rownames(ex$expr))

# keep the 2000 most variable genes, this doesnt use survival so nothing leaks
avg <- rowMeans(tc$logcpm[shared, ])
well <- shared[avg > median(avg)]
m <- apply(tc$logcpm[well, ], 1, mad)
top <- names(sort(m, decreasing = TRUE))[1:nvar]

# scale each gene inside its own dataset
xt <- scale(t(tc$logcpm[top, ]))
xe <- scale(t(ex$expr[top, ]))
good <- apply(xt, 2, function(v) all(is.finite(v))) & apply(xe, 2, function(v) all(is.finite(v)))
xt <- xt[, good]
xe <- xe[, good]

# gene names with dashes break the model formulas so make them safe
# and keep a lookup back to the real names
orig <- colnames(xt)
safe <- make.names(orig, unique = TRUE)
colnames(xt) <- safe
colnames(xe) <- safe
gene_map <- setNames(orig, safe)

clin <- tc$clin[match(rownames(xt), rownames(tc$clin)), ]
time <- clin$time_months
status <- clin$status
saveRDS(list(Xt = xt, Xe = xe, clin_t = clin, clin_e = ex$clin, gene_map = gene_map),
        "data/model_input.rds")
cat("patients:", nrow(xt), " genes:", ncol(xt), "\n")

# c-index, higher risk should mean dies sooner
cidx <- function(risk, t, s) {
  unname(concordance(Surv(t, s) ~ risk, reverse = TRUE)$concordance)
}

# p value for each gene on its own
uni_p <- function(x, t, s) {
  apply(x, 2, function(g) {
    f <- tryCatch(suppressWarnings(coxph(Surv(t, s) ~ g)), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[1, 5]
  })
}

# split into folds, keeping deaths and censored balanced
make_folds <- function(s, nfold) {
  f <- integer(length(s))
  for (v in unique(s)) {
    i <- which(s == v)
    f[i] <- sample(rep_len(1:nfold, length(i)))
  }
  f
}

res <- list()
sigs <- list()
n <- 0

for (r in 1:nrep) {
  fold <- make_folds(status, nfold)

  for (k in 1:nfold) {
    tr <- which(fold != k)
    te <- which(fold == k)
    xtr <- xt[tr, ]
    xte <- xt[te, ]
    ytr <- Surv(time[tr], status[tr])

    # gene filter using the training patients only
    p <- uni_p(xtr, time[tr], status[tr])
    sel <- names(p)[p < pcut]
    if (length(sel) < 20) sel <- names(sort(p))[1:20]

    # cox lasso, try to land on 5-10 genes
    fit <- glmnet(xtr[, sel], ytr, family = "cox")
    ix <- which(fit$df >= 5 & fit$df <= 10)
    if (!length(ix)) ix <- which(fit$df >= 5)
    li <- if (length(ix)) ix[1] else length(fit$df)
    lam <- fit$lambda[li]
    risk_l <- as.numeric(predict(fit, newx = xte[, sel], s = lam, type = "link"))
    cf <- coef(fit, s = lam)
    g_l <- rownames(cf)[as.numeric(cf) != 0]

    n <- n + 1
    res[[n]] <- data.frame(rep = r, fold = k, model = "Cox-LASSO",
                           cindex = cidx(risk_l, time[te], status[te]),
                           n_genes = length(g_l))
    sigs[[n]] <- data.frame(model = "Cox-LASSO", gene = g_l)

    # random survival forest
    # rank the top 100 genes by importance then refit on the best 10
    sel_rf <- head(names(sort(p[sel])), 100)
    d0 <- data.frame(time = time[tr], status = status[tr], xtr[, sel_rf, drop = FALSE])
    rf0 <- rfsrc(Surv(time, status) ~ ., data = d0, ntree = 300,
                 importance = "permute", seed = -42)
    g_r <- names(sort(rf0$importance, decreasing = TRUE))[1:10]

    d1 <- data.frame(time = time[tr], status = status[tr], xtr[, g_r, drop = FALSE])
    rf1 <- rfsrc(Surv(time, status) ~ ., data = d1, ntree = 300, seed = -42)
    risk_r <- predict(rf1, newdata = data.frame(xte[, g_r, drop = FALSE]))$predicted

    n <- n + 1
    res[[n]] <- data.frame(rep = r, fold = k, model = "RandomSurvivalForest",
                           cindex = cidx(risk_r, time[te], status[te]),
                           n_genes = length(g_r))
    sigs[[n]] <- data.frame(model = "RandomSurvivalForest", gene = g_r)
  }
  cat("done with repeat", r, "of", nrep, "\n")
}

R <- do.call(rbind, res)
S <- do.call(rbind, sigs)
write.csv(R, "results/07_cv_results.csv", row.names = FALSE)
S$gene <- gene_map[S$gene]
write.csv(S, "results/07_cv_signature_genes.csv", row.names = FALSE)

cat("\n===== CROSS-VALIDATED C-INDEX (TCGA only) =====\n")
for (mod in unique(R$model)) {
  v <- R$cindex[R$model == mod]
  cat(sprintf("%-22s mean %.3f | sd %.3f | median genes %d\n",
              mod, mean(v), sd(v), median(R$n_genes[R$model == mod])))
}

w <- reshape(R[, c("rep", "fold", "model", "cindex")],
             idvar = c("rep", "fold"), timevar = "model", direction = "wide")
d <- w[["cindex.Cox-LASSO"]] - w[["cindex.RandomSurvivalForest"]]
cat(sprintf("\nmean difference (lasso - forest): %.3f | lasso won %d of %d folds\n",
            mean(d), sum(d > 0), length(d)))

cat("\nmost often selected genes:\n")
tab <- sort(table(S$gene), decreasing = TRUE)
print(head(tab, 15))
